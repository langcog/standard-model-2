# Active refit plan (2026-05-15)

Self-contained snapshot of what we're doing and where we're going, so
either we can pick up after a context break or hand off to a fresh
session.

## Current state of the world

**Notation conventions (settled):**
- `δ_j` = per-word difficulty (IRT-standard); Stan code uses `delta_j[j]`.
- `δ` = population scaling-exponent shift; Stan code uses `delta`.
- `κ_pop = 1 + δ` = population scaling exponent.
- `κ_i = κ_pop + ζ_i` = per-kid scaling.
- `σ_ζ` = SD of per-kid scaling deviation.
- `s` = population trajectory anchor (now pinned at 6 mo, see below).
- `s_i` = per-kid trajectory phase offset (signed, sum-to-zero centered).
- `σ_s` = SD of `s_i`.
- `α_i` = per-kid efficiency; `σ_α` = SD of `α_i`.
- `ξ_i` = per-kid intercept = `log r_i + log α_i`.
- `π_α = σ_α² / (σ_α² + σ_r²)` = efficiency share of intercept variance.

**s pinning convention (settled):** `s_prior_mean = 6, s_prior_sd = 0.5`
in `model/R/config.R` `DEFAULT_PRIORS`. Set this way (versus the old
s ≈ 0.5) so signed `s_i` deviations around `s` yield positive
effective onsets (~6 ± 2 mo) — interpretively cleaner. **Cost:** the
(s, δ) ridge from §20 means κ_pop will drop substantially. Estimated
new κ_pop is in the 5–7 range (was 10.6 with s = 0.5). Qualitative
claim (κ_pop > 1, super-linear) is invariant.

**Active Stan files:**
- `log_irt_long.stan` — base longitudinal IRT accumulator.
- `log_irt_long_si_signed.stan` — adds signed-`s_i` + (σ_total, p_zeta)
  reparam. Headline variant.
- `log_irt_long_si_corr.stan` — trivariate LKJ (sidebar finding).
- `log_irt_long_lmm.stan`, `log_irt_long_proc.stan`,
  `log_irt_long_io_accel.stan`, `log_irt_long_io_comp.stan` — IO/proc.
- `log_irt.stan`, `log_irt_io.stan` — cross-sectional base.

**Active variants (5-panel build):**
1. `long_demo_pure` — pure accumulator (σ_α=σ_ζ=σ_s=0, δ=0)
2. `long_demo_alpha` — + σ_α free (δ still 0)
3. `long_demo_kappa_pop` — + δ free (σ_ζ still 0) — NEW VARIANT
4. `long_no_freq_slopes` — + σ_ζ free (was "M_best")
5. `long_no_freq_slopes_si_signed` — + σ_s free; signed-`s_i` — HEADLINE

## 1. Optimizations under test

**Threading benchmark** (job 25132788 on Sherlock, currently PENDING):
- `long_no_freq_slopes_si_signed` refit at `--cpus-per-task=32`,
  `STAN_THREADS_PER_CHAIN=8`.
- Baseline: 2h53m at 4 threads/chain on the current I=200, J=198 bundle.
- Target: ≥1.7× speedup, ideally 2×.

**Decision tree from benchmark result:**
- Speedup ≥2.5×: use 8 threads/chain for family refit. I=800 × J=650
  projects to ~12-14h. Within slurm 18h window. ✓
- Speedup 1.7–2.5×: use 8 threads/chain. I=600 × J=650 projects to
  ~12h. ✓
- Speedup <1.7×: fall back to 4 threads/chain, smaller bundle (I=400,
  J=650).

**Other speedups available (deferred unless we need them):**
- Higher adapt_delta only when needed for divergent transitions
  (currently fine at 0.95).
- Init values from a prior fit (skip cold-start exploration).
- 16 threads/chain on the 2 newest Sherlock nodes (huge queue wait).

**Hardware available (`normal` partition):**
- 106 × 24-core nodes
- 98 × 32-core nodes  ← target for benchmark
- 62 × 20-core nodes
- 2 × 64-core nodes (long queue, only use if compute is critical)

## 2. 5-panel additive build for English

**Once benchmark lands**, execute:

### Step 1: Build scaled bundle

```bash
Rscript model/scripts/prepare_longitudinal_data.R "English (American)" 800 650
```

Writes `fits/long_subset_data.rds` (overwrites current I=200, J=198).
Stratified-by-median-admin-age sampling for kids (4 bins);
stratified-by-(lexical_class × log_p_quartile) sampling for items.

Expected bundle dimensions: I=800, A≈3000 admins, J=~650, N≈1.8M obs.

### Step 2: Submit 5 fits on Sherlock

Each variant submits via `sherlock/long_fit.slurm <variant> english`.
Required env vars:
- `STAN_ITER=4000`, `STAN_WARMUP=2500`
- `STAN_ADAPT_DELTA=0.95`
- `--cpus-per-task=32`, `STAN_THREADS_PER_CHAIN=8` (assuming benchmark passes)
- `--time=18:00:00` (max slurm window)

Variants to submit (in this order; can be parallel if fairshare allows):
1. `long_demo_pure`
2. `long_demo_alpha`
3. `long_demo_kappa_pop`
4. `long_no_freq_slopes`
5. `long_no_freq_slopes_si_signed`  (most parameters; slowest)

### Step 3: Extract summaries + LOO + delta_j

For each completed fit:
```bash
srun --time=01:30:00 --mem=48G --cpus-per-task=1 --ntasks=1 \
  bash -c "ml R && export STANDARD_MODEL_FITS_DIR=\$SCRATCH/standard_model_2/fits && \
    Rscript sherlock/extract_summary_table_only.R <tag> && \
    Rscript sherlock/extract_scalar_draws.R <tag> && \
    Rscript sherlock/extract_psi_slim.R <tag> && \
    Rscript sherlock/extract_loo_thinned.R <tag>"
```

Then rsync to `fits/summaries/`:
```bash
rsync sherlock:/scratch/users/mcfrank/standard_model_2/summaries/<tag>.{summary,draws}.rds \
      sherlock:/scratch/users/mcfrank/standard_model_2/summaries/<tag>_psi.csv \
      sherlock:/scratch/users/mcfrank/standard_model_2/summaries/<tag>.loo.rds \
      fits/summaries/
```

### Step 4: Rebuild figure

`Rscript model/scripts/quantile_demo.R` — script already configured to
pick up the 5 variants. Outputs `outputs/figs/longitudinal/quantile_demo.png`.

### Step 5: Update explainer + experiments.md

The explainer currently quotes posterior values from the I=200 fits
(σ_α=1.56, σ_ζ=3.51, σ_s=1.40, δ=9.62, κ_pop=10.62). These will need
to update for the new I=800 fits with s=6 pinning. Expected:
- σ_α, σ_ζ values likely similar (data-identified, robust to scale).
- σ_s might shift (s_i is partially identified — more kids better).
- δ will shift substantially because s=6 (was s=0.5).
- κ_pop = 1 + δ — new value to determine.

## 3. Followup refit sequence

**B. Norwegian.** Same 5-variant family, applied to the Norwegian
longitudinal bundle. Norwegian has 1,562 kids / 5,702 admins / 722
items — bigger than English. Scale-down options:
- I=800 × J=full ~720 items → ~1.6M rows (similar to English scale).
  Should run in ~12-14h per fit.
- Alternative: I=600 × J=720 if compute is tight.

Bundle prep: `Rscript model/scripts/prepare_longitudinal_norwegian.R
"Norwegian" 800 720`. Writes `fits/long_subset_data_nor.rds`.

Variant submission identical to English but with `norwegian` dataset
key. The §19 cross-sample replication table will be updated with new
Norwegian π_α at s=6.

**C. IO models (BabyView + SEEDLingS).** Smaller samples (BabyView
N=20, SEEDLingS N=44) but specialized models with kid-level input
rate (`log_irt_long_io_accel.stan`) or comp+std channels
(`log_irt_long_io_comp.stan`).

Required variants:
- BabyView: `io_no_freq_slopes_si_signed` (need to add this — does the
  io Stan file support the si_signed extension? Check.)
- SEEDLingS: `io_comp_no_freq_slopes_si_signed` (similar question).

For both: with s=6 pinning, expect π_α to change slightly from the
§19 reported values (0.84 for BabyView, 0.94 for SEEDLingS). The
qualitative claim (input rate ~5-15% of variance) should hold.

**Bundle status:** Already on disk (`babyview_subset_data.rds`,
`seedlings_subset_data.rds`). Probably don't need to rebuild unless we
want to incorporate new data.

**D. Peekbank (processing channel).** `log_irt_long_proc.stan`. Smaller
N (62 kids). Doesn't need to be scaled. With s=6 pinning, π_α may
shift slightly from §19's 0.90.

**Important caveat for B/C/D:** All consumers of σ_α / π_α / κ_pop need
to be re-reported at the new s=6 pinning for the §19 cross-sample
table to be internally consistent. Don't mix-and-match s=0.5 fits with
s=6 fits in the same comparison.

## 4. After all refits land

- **Update §19 cross-sample replication table** in experiments.md with
  all-s=6 values.
- **Update §22 mixing diagnostics** with the s=6 numbers for
  slopes_si_signed.
- **Update model_explainer.tex** "Per-child onset" section with new
  σ_s / δ / κ_pop values.
- **Rebuild slide 21** (4-panel demo → 5-panel) in
  `outputs/slides/standard_model_external.tex`.
- **Update slide 17** (unified model comparison) with the new κ_pop.

## 5. Headline numbers as of 2026-05-15 (will be replaced)

| Variant | σ_α | σ_ζ | σ_s | δ | κ_pop | π_α | Rhat | N |
|---|---|---|---|---|---|---|---|---|
| Old M_best (long_no_freq_slopes, s=0.5) | 1.83 | 3.47 | — | 9.39 | 10.39 | 0.91 | 1.04 | 200 |
| signed_si (s=0.5) | 1.56 | 3.51 | 1.40 | 9.62 | 10.62 | 0.89 | 1.02 | 200 |
| **TARGET (s=6, I=800)** | TBD | TBD | TBD | TBD | **~5-7** | TBD | <1.05 | 800 |

The headline κ_pop dropping from 10.6 → ~5-7 with the s=6 re-pin is
the most consequential change. Be ready to update every reference to
the κ_pop ≈ 10 number across the explainer, slides, and any draft
abstracts.

## 6. Risk register

- **Benchmark queue wait**: 32-cpu request is held back by fairshare;
  may need to wait several hours before benchmark runs.
- **18h slurm window**: I=800 × J=650 projects to 12-14h at 8
  threads/chain, but per-iteration cost depends on chain mixing.
  signed_si is well-mixed in current bundle; should still be well-mixed
  at larger scale, but not guaranteed.
- **Mixing at larger scale**: more kids = more s_i and ζ_i parameters
  to sample. Could expose new geometry issues. Plan B: revert to
  I=500 if Rhat > 1.1.
- **Norwegian bundle**: bigger data; may need smaller I to fit in
  slurm window.

## Sentinel for autonomous handoff

If a fresh session needs to pick this up, the entry points are:
- This file (read first for current state).
- `outputs/experiments.md` §22 for the mixing-fix backstory.
- `outputs/model_explainer.tex` "Per-child onset s_i" section for the
  signed-s_i mechanics.
- `model/R/helpers.R` `variant_hyperpriors()` for active-variant
  definitions.
- `model/scripts/quantile_demo.R` for figure pipeline.
- Latest fits + LOO files are in `fits/summaries/`.
