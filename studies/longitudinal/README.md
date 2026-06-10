# longitudinal — Bayesian log-age IRT (the core child model)  ·  Figs 3,4 + demographics

**Provenance stub.** Code lives in the shared engine (`model/`), indexed here.

- **Model:** `model/stan/log_irt_long.stan` (1PL log-age IRT; θ_i(t)=ξ_i+κ_i·log((t−s)/a0))
  and `log_irt_long_dprime.stan` (D′: input on the slope channel).
- **Fit driver:** `model/scripts/fit_longitudinal.R`; analysis `analyze_longitudinal.R`
- **Datasets:** EN (pooled Thal+Smith+Marchman), NO, JP → `fits/long_*` (e.g. `long_no_freq_slopes*`)
- **Demographics BLUPs:** `model/scripts/analyze_longitudinal.R` → `paper/cache/blups_demographics.rds`
  (feeds the longitudinal arm of the demographics figure; xsec arm = `studies/cross_sectional_demographics/`)
- **Exposures-to-learn (Fig 4):** `model/scripts/exposure_to_learn.R` → `paper/cache/fig4_exposure.rds`
- **Input share (Fig 3 A–C):** `paper/build_input_cache.R` reads the `long_*` summaries.

Narrative: see [`/journal/experiments.md`](../../journal/experiments.md).
