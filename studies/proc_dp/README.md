# proc_dp — processing (LWL reaction time) study  ·  Fig 3E

**Provenance stub.** Code lives in the shared engine (`model/`), indexed here.

- **Data prep:** `model/scripts/prepare_proc_dp_bundle.R`
  → `fits/proc_dp_all_subset_data.rds` (I=226, 3 datasets: AM2018, FM2012, FMW2013;
  observed LENA in AM2018/FMW2013, imputed for FM2012; LWL RT measured by looking-while-listening)
- **Model:** `model/stan/log_irt_long_proc_dp.stan` (regression ladder: RT → efficiency ξ and acceleration κ)
- **Fit driver:** `model/scripts/fit_proc_dp.R` (ladder D′0–D′3)
- **Cluster:** `cluster/sherlock/proc_dp_fit.slurm` (--mem=96G; loo(cores=1))
- **delta_j extractor:** `cluster/sherlock/extract_proc_deltaj.R` → `fits/summaries/proc_dp1_all_psi.csv`
- **Selected rung:** D′1 — processing predicts efficiency (β_ξ); both acceleration rungs null.
- **Figure:** `paper/build_input_cache.R` → `paper/cache/fig3_input.rds` (panel E); rendered in `paper/standard_model.qmd` (`fig-io-partition`).

Narrative: see [`/journal/experiments.md`](../../journal/experiments.md) (proc_dp ladder entries).
