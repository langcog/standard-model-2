# Imputed-panel refits on clean bundles — runbook (2026-08-15)

**Why.** The EN/NO longitudinal bundles were rebuilt on the acceleration repo's
data standard (clean child key `dataset::study_internal_id`, uni_lemma
cross-form harmonization, crater/jump QC, wordbankr 2.0 / Redivis dataset
v2.0). All six imputed fits must be redone on them. Full context:
`journal/paper_models_provenance.md`, 08-15 currency flags.

**Compute: Sherlock only** (owners partition now available; GCP retired for
this. Local cmdstan is CrowdStrike-killed — never fit locally).

## The six fits — SUBMITTED 2026-08-15 (jobs 39218376–39218385, owners, 48h)

**σ_r anchor decision (MCF approved 08-15):** pins = **{0.35, 0.44, 0.58}**.
0.35 = the direct enct/io_count fits' own σ_r estimate (0.35 [0.32, 0.39],
study-centered → population lower bound); 0.44 = channel-matched population
re-anchor (main dashed pin); 0.58 = literature upper. **0.53 (Sperry CDS)
retired** — no baseline (un-suffixed) refits; `build_fig_io_cache.R` now
reads σ_xi from the 0p44 fits and the band is [0.32, 0.58], so it will
error (not silently go stale) until the new summaries land.

| tag | dataset | σ_r pin |
|---|---|---|
| `long_no_freq_slopes_sigmaR_0p35` | english | 0.35 |
| `long_no_freq_slopes_sigmaR_0p44` | english | 0.44 |
| `long_no_freq_slopes_sigmaR_0p58` | english | 0.58 |
| `long_no_freq_slopes_norwegian_sigmaR_0p35` | norwegian | 0.35 |
| `long_no_freq_slopes_norwegian_sigmaR_0p44` | norwegian | 0.44 |
| `long_no_freq_slopes_norwegian_sigmaR_0p58` | norwegian | 0.58 |

*Confirmed: the rebuilt bundles pin σ_r = 0.534 (Sperry-anchored baseline),
so the anchor set stays {0.44 pin, 0.53 baseline, 0.58 pin} and **all six
fits are needed** (build_fig_io_cache.R plots anchors at those three values
and skips missing summaries gracefully).

**Bundle validation (2026-08-15):** rebuilt EN = I 3105 / A 7329 / J 681 /
N 4,253,076 — exactly the sum of the acceleration repo's clean thal+smith+
marchman bundles (Marchman post-QC 2136 ✓, QC removed 90 admins ✓); rebuilt
NO = I 1630 / A 6446 / N 4,442,722 — exactly their norwegian bundle (QC 50 ✓).
The two repos are on identical data. EN a0 moved 19 → 18 (median admin age
with the WG arm restored); NO a0 = 25 unchanged.

## Launch pattern (already run 08-15; kept for resubmission)

Bundles + scripts were rsynced to `$HOME/standard_model_2` on Sherlock
(the clone's git state is old/dirty by design — nothing was committed;
resync with rsync, not git pull, until the reorg is committed):

```bash
for v in 0.35 0.44 0.58; do for d in english norwegian; do
  sbatch -p owners --requeue -t 48:00:00 \
    --export=ALL,STAN_SIGMA_R_OVERRIDE=$v \
    cluster/sherlock/long_fit.slurm no_freq_slopes $d
done; done
```

Be generous with `-t` (see feedback memory: short walltimes silently kill
fits; normal QOS caps at 48h, `--qos=long` beyond). The clean EN bundle is
~2× the old one (Marchman WG arm restored) — expect the EN fits to run
closer to the old Norwegian times than the old English times.

## After the fits land

1. Extract with the existing `cluster/sherlock/extract_summaries.R` flow →
   `fits/summaries/<tag>.{draws,summary}.rds`, pull back locally.
2. `Rscript paper/build_cache.R` (§8 io_partition; §7 input_rate_table
   unaffected) and `Rscript paper/build_fig_io_cache.R` (panel A curves +
   anchors; panel B unaffected — the joint enct fit is current).
3. Re-render `paper/input_paper.qmd`; the imputed section inline values
   (`io_share()`, π_α) update from the caches.

## Queue state (2026-08-15 evening)

Everything EXCEPT the six imputed refits above is now submitted or done —
the imputed refits are deliberately held pending the σ_r anchor-set
discussion (MCF: current anchors {0.44/0.53/0.58} may be wrong; note the
direct enct/io_count fits estimate σ_r ≈ 0.35, below the band's lower edge).

- Sherlock **39216615**: proc-only fit (`log_irt_long_proc_count.stan`, new,
  syntax-checked; mirror of io_count) → tag `joint_io_proc_lean_proc_d2_enct`.
- Sherlock **39216878**: glmer ladder, 35-cell array (5 clean datasets × 7
  models) on the rebased extraction — extracts match the acceleration
  bundles exactly (marchman 2136/QC90, smith 316, thal 653, norwegian
  1630/QC50; japanese 187 kids new-standard). Both on `-p owners --requeue`.
- io_count (no-proc) summary+draws recovered from Sherlock scratch →
  `fits/summaries/joint_io_proc_lean_io_d0_enct.*` (eff_input_k 0.762
  [0.082, 1.41], matches entry 39).
- `01_extract_one.R` REBASED onto `long_items.rds` (supersedes 01b);
  pull now includes Japanese; `data_check_io_panels.R` repointed to the
  english_count bundle; glmer slurm+submit paths fixed (were pre-reorg).

Remaining local chores after fits land: `paper/build_cache.R` §3 must join
demographics via the saved `ckey_map` in each `data_<slug>.rds` (ckey no
longer equals a Wordbank child_id) + a wordbankr-2.0 pass; wire an io_count/
proc_count cache once the separate-vs-joint figure design settles.

## Still queued behind this (separate work)

- **glmer ladder re-extraction** for the demographics BLUPs: recommended
  design = add Japanese to `pull_longitudinal.R` LANGUAGES and rebase
  `studies/glmer_ladder/01_extract_one.R` onto `long_items.rds` (one
  extraction, one standard; kills the child_id/QC divergence for good).
  Then Sherlock `glmer_ladder.slurm` refits ×5 → rebuild
  `blups_demographics.rds` (its build_cache.R §3 demo join must switch to
  ckey and get a wordbankr-2.0 pass).
- **Separate proc-only fit** (mirror of `log_irt_long_io_count.stan`) +
  re-extract the io_count no-proc summary from Sherlock (entry 39, job
  31764182) — unaffected by the Wordbank cleanup (io/proc datasets).
