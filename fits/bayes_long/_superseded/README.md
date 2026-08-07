# Superseded `bayes_long` fits — safe to delete

Archived 2026-07-08. See `journal/experiments.md` §40 / §40.1 for the full story.
Everything here is reproducible from `studies/bayes_long/{00_prepare_bundles.R,01_fit.R}`
+ the current `_a3` bundles, so nothing unique is lost by deleting this directory.

Two categories:

1. **Base 2+-admin generation (no `_a3` suffix).** The first `bayes_long` fits, at
   ≥2 administrations/child. Superseded by the ≥3-admin (`_a3`) filter, which
   resolved the σ_b sparse-design inflation (Smith σ_b 8→5) and let Norwegian M3
   converge. Includes `bundle_<slug>.rds` and `summaries/<slug>_m{0-3}.*` for all
   five datasets (Norwegian base has only m0–m2 — its base m3 was a killed
   stuck-chain funnel).

2. **Pre-QC `_a3` M3 fits: `marchman_a3_m3.*`, `norwegian_a3_m3.*`.** Fit on the
   `_a3` bundles *before* the monotone-vocabulary QC filter (commit `ab8b81c`)
   removed the impossible declining trajectories (Marchman −21, Norwegian −13).
   Their σ_b is inflated by those records; the post-QC reruns replace them.
   (`thal/smith/japanese_a3_m3` had 0 QC drops → unchanged → kept live, not here.)

**Current, kept live** in `fits/bayes_long/`: the five `_a3` bundles, the three
unchanged M3 fits (`{thal,smith,japanese}_a3_m3`), and the post-QC M0–M3 reruns as
they land from Sherlock.
