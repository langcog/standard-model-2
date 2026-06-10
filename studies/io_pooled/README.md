# io_pooled — input-observed study  ·  Fig 3 A–D

**Provenance stub.** Code lives in the shared engine (`model/`), indexed here.

- **Per-dataset prep:** `model/scripts/prepare_io_dataset.R` (am2018/fmw2013),
  `model/scripts/parse_seedlings_cdi.R`, `analyze_babyview.R`
- **Pooled prep:** `model/scripts/prepare_io_pooled.R` → `fits/io_pooled_subset_data.rds`
  (4 datasets: BabyView head-cam, SEEDLingS/AM2018/FMW2013 LENA; global delta_j anchored to EN longitudinal)
- **Model:** `model/stan/log_irt_io_pooled.stan` (+ `_gamma_add/_mult` variants for the slope channel)
- **Main fit:** `model/scripts/fit_io_pooled_widedelta.R` → `fits/io_pooled_widedelta.rds`
- **Cluster:** `cluster/sherlock/io_pooled_*_fit.slurm`, `cluster/gcp/`
- **Headline:** input is a small share of between-child intercept variance (~3%), convergent with Coffey 2026 meta (4–7%).
- **Figure:** `paper/build_input_cache.R` → `paper/cache/fig3_input.rds` (panels A–D); `fig-io-partition` chunk.

Narrative: see [`/journal/experiments.md`](../../journal/experiments.md) (io / io_pooled entries).
