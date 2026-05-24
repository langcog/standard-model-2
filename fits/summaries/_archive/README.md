# Retired fit summaries

Archived 2026-05-23. These summary / draws / loo / psi files belong to
fit variants no longer cited in the slide deck or referenced by the
current M_best pipeline.

## Categories

### `s_i` family (post-§22 retire)

- `long_no_freq_slopes_si.*` — half-normal `s_i` variant
- `long_no_freq_slopes_si_corr.*` — LKJ-correlated `(ξ, ζ, s_i)`
- `long_no_freq_slopes_si_reparam.*` — `(σ_total, p_zeta)` reparam
- `long_no_freq_slopes_si_signed.*` — signed-normal `s_i` (the
  headline before excision)
- `long_no_freq_slopes_si_signed_psi.csv`, `long_no_freq_slopes_si_corr_psi.csv`,
  `long_no_freq_slopes_si_reparam_psi.csv`, `long_no_freq_slopes_si_psi.csv`
- `long_no_freq_si_only.*` — α + `s_i` (no per-kid slopes)

See `experiments.md` §23 for the retire rationale.

### Pre-M_best lineage (§14 nested family + baseline ablations)

- `long_baseline.*`, `long_baseline_norwegian.*` — Rasch only
- `long_2pl_slopes.*` etc. — 2PL with discrimination
- `long_m0.*`, `long_m0_norwegian.*` — m0 (time = 0)
- `long_m1.*`, `long_m1_norwegian.*` — m1 (rate only)
- `long_m1_time_only.*` (+ Norwegian) — m1 with time on log scale only
- `long_class_beta_slopes.*` — per-class log-p slope (β_c free)
- `long_no_class_slopes.*` — class structure removed (data override)
- `long_fix_delta_slopes.*` — δ pinned at 0 (no acceleration)
- `long_lmm_slopes_norwegian.*` — linear mixed model comparison
- `long_proc_2pl_slopes.*`, `long_proc_slopes.*` — older proc variants
- `free_s_no_freq_slopes.*` — free-s robustness check (cited in §20 of
  experiments.md but not in the current deck)

### Backup / stale snapshots

- `long_no_freq_slopes_psi.csv.stale_may20` — pre-cleanup snapshot
- `long_no_freq_slopes_norwegian.*.I200_backup` — pre-bundle-upgrade NO I=200 results, kept after we scaled to I=500

## What's NOT archived (still active in `..`)

- `long_no_freq_slopes.{summary,draws}.rds`, `_psi.csv` — EN M_best (slide 23)
- `long_no_freq_slopes_norwegian.{summary,draws}.rds`, `_psi.csv` — NO M_best (slide 23)
- `long_no_freq_slopes_english_I200.{summary,draws}.rds`, `_psi.csv` — I=200 M_best (slides 19, 20)
- `long_demo_pure.*`, `long_demo_alpha.*`, `long_demo_kappa.*`, `long_demo_kappa_pop.*` — 4-panel architecture demo (slides 19, 20)
- `io_no_freq_slopes.*`, `io_no_freq_slopes_seedlings.*`, `io_comp_*.*`, `io_std_*.*`, `io_comp_std_*.*` — IO fits (slides 32–33, 36)
- `long_proc_no_freq_slopes.{summary,draws}.rds`, `.draws_full.rds`, `_psi.csv` — Peekbank/proc (slides 34–36)
- `long_no_freq_slopes_sigmaR_0p*.*` — σ_r sensitivity refits (slide 28 black dots)
- `long_demo_alpha_psi.csv`, `long_demo_kappa_psi.csv`, `long_demo_pure_psi.csv` — for the 4-panel demo

See `outputs/PROVENANCE.md` for the full slide-asset → fit map.
