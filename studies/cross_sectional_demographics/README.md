# cross_sectional_demographics

Cross-sectional replication + extension of the paper's longitudinal
demographic analysis (sex & maternal education → efficiency vs. acceleration
of vocabulary growth), across ~30 Wordbank languages.

## Model

Per language, one admin per child (monolingual-TD: bilingual + clinical
datasets excluded; **not** restricted to `is_norming`), item-level Rasch GLMM:

```
produces ~ predictor * log(age / a0) + (1 | item) + (1 | child)
```

`predictor` = sex (Male vs Female) or maternal education (z-scored years).
The predictor **main effect** = effect on **efficiency** (intercept at `a0`);
the `predictor:log_age` **interaction** = effect on **acceleration** (slope).
Raw `log(age/a0)` is used so coefficients match the longitudinal raw-BLUP
regressions and are directly comparable.

## Files

- `00_build.R` — pulls one admin/child per language, scores production, fits
  the two models, meta-analyzes, and writes the caches. Resumable (per-language
  frames + fits cached). Wordbank pulls are the bottleneck.
- `cross-sectional_demographics.qmd` — notebook: scatter data checks (sex,
  maternal ed), cross-sectional forest + meta figure, combined
  cross-sectional + longitudinal figure, anomaly diagnostics.
- `cache/fits.rds` (committed) — per-language fits + meta + longitudinal effects.
- `cache/scatter.rds` (committed) — per-child production proportion for scatters.
- `cache/frames/`, `cache/fits/` (gitignored) — regenerable item-level frames
  and per-language fit summaries.

## Key findings (see notebook / experiments.md #31)

- **Sex → efficiency** is near-universal: significant female efficiency
  advantage in ~26/31 languages, ~zero on acceleration. Matches the
  longitudinal estimate closely (e.g. Norwegian −0.86 x-sec vs −0.89 long).
- **Maternal ed** is heterogeneous and noisier; cross-sectional *acceleration*
  is inflated vs. longitudinal (validation: ~1.5–2.3× per language), so it is
  an upper bound. Some per-language matEd estimates are unstable where the
  education distribution is skewed (e.g. Spanish-European, French).

## Reproduce

```bash
Rscript studies/cross_sectional_demographics/00_build.R         # ~hours first run (pulls); cached after
quarto render studies/cross_sectional_demographics/cross-sectional_demographics.qmd
```
