---
name: stan-mixing-failures
description: Diagnose and recover from Stan model fits that completed sampling but didn't converge (high Rhat, low n_eff, treedepth saturation). Especially relevant to the standard_model_2 family build where init_from a degenerate source fit parks chains at a parameter boundary. Use whenever a fit finishes with Rhat > 1.1, n_eff < 100, large between-chain time spread, or "100% of transitions hit max treedepth" warnings.
---

# Stan mixing failures: the cold-start / boundary-parking pattern

## Symptoms

A fit that sampled all the way through but came back with:

- **Rhat 2.5–3.0** on key parameters (target: <1.05)
- **n_eff ≈ 5** (target: >>100)
- **100% of transitions hit max treedepth=10** in the Stan summary
- **All chains had E-BFMI < 0.3** (energy diagnostic failure)
- **Large between-chain wall-clock spread** during sampling — fastest chain finishes 25-50% earlier than slowest (e.g., 10h vs 14h)

These usually co-occur. If you see any of them, the others are likely true too.

## The most common root cause: warm-start from a degenerate prior fit

In `standard_model_2`'s family build, each variant adds one parameter to the previous one. Warm-starting the next variant from the previous one's posterior summary parks chains where the new parameter's prior is degenerate.

**Concrete example we hit (2026-05-18):**

`long_demo_pure` pins σ_α = 0. `long_demo_alpha` frees σ_α but is initialized from `demo_pure`'s posterior, where σ_α ≈ 0. All 4 chains then start at the σ_α = 0 boundary. The geometry there is highly correlated (σ_α and σ_r trade off ambiguously), and the sampler can't escape — every transition saturates max_treedepth trying. Result: Rhat ~2.7, n_eff ~5.

This is **not** fixed by:
- Bumping `max_treedepth` (the chains still can't move, they just spend more time per trajectory — wall-clock goes 6–8× slower)
- Bumping `adapt_delta` (same — smaller step, more leapfrog, same stuck region)
- More warmup iters (warmup is doing the right thing, the init is wrong)

## The fix: cold-start the variant

Drop `init_from` for any variant whose source fit pinned a parameter the target frees. In the family build that's:

- `long_demo_alpha`: cold-start (source `demo_pure` pins σ_α)
- `long_demo_kappa_pop`: cold-start (source `demo_alpha` pins κ_pop)
- `long_no_freq_slopes`: **keep** warm-start (source `demo_kappa_pop` has freed κ_pop)
- `long_no_freq_slopes_si_signed`: **keep** warm-start (source `slopes` is non-degenerate)

`gcp/run_family.sh` already encodes this policy via the `COLD_START_VARIANTS` list.

## Don't reach for adapt_delta / max_treedepth defensively

`adapt_delta=0.97` and `max_treedepth=12` were tried as belt-and-suspenders after the failure. They cost **~6–8× per-iteration wall-clock** on healthy chains (more leapfrog steps allowed, smaller step size) without fixing the underlying boundary issue. They should only be used **after** confirming a fit fails at defaults with cold-start.

Defaults are `adapt_delta=0.95, max_treedepth=10`. Leave them alone unless a cold-start fit also fails to converge.

## Spread between chain finish times is the early warning

If chain 1 finishes sampling at 10h and chain 3 is still at 70% sampling 13h in, you're likely already in the boundary-parking regime. Don't wait for the final summary — abort and switch to cold-start.

`ps -eo etime,cmd | grep log_irt` and the iter-progression lines in the log show the spread clearly.

## Recovering a stalled family

If the family driver crashed mid-run (e.g., `set -e` killed it after a bad fit), the standard recovery is:

1. Delete the bad fit's `.rds`: `rm fits/long_demo_alpha.rds` (or whatever variant failed)
2. Re-run `gcp/run_family.sh <lang>` — its skip-if-exists logic will skip variants whose `.rds` is still on disk and pick up from the failed one.

If the LOO extract was the only thing that crashed (sampling succeeded, `.rds` saved), use:
```
Rscript sherlock/extract_loo_thinned.R <tag>
```
which uses a chunked PSIS implementation that survives N≥2M obs. See `psis-loo-large-datasets` skill.

## Don't scp a script that's currently running

If you `gcloud compute scp gcp/run_fit.sh` to a VM where `run_fit.sh` is currently running, bash will re-read the file mid-execution and the running shell can throw a syntax error at the next line boundary. This crashed EN's family driver between sampling and extracts on 2026-05-18.

If you must patch a script during a run, patch on a different VM, or wait for the current `run_fit.sh` invocation to finish.

## A note on Rhat-with-3-admins-per-kid

In English (median 3 admins/kid), per-child slopes `ζ_i` are under-identified. The population term `δ` ends up pooling what Norwegian's data (median 8 admins/kid) lets fall into per-child variance. So the `δ` vs `ζ_i` decomposition shifts across languages even though the aggregate `κ_i` distribution is the same. Don't read too much into a population-level term shift between languages with very different admin densities — read the aggregate.
