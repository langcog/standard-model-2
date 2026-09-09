# paper/tables — LaTeX tables for the Overleaf manuscript / SI

One script per table, reading committed caches only; the `.tex` output in `out/`
is committed and uploaded to Overleaf (`\input{}` it inside a table environment;
needs `booktabs`).

| table | script | reads | output |
|---|---|---|---|
| Datasets and analysis samples (SI Table S1: Wordbank longitudinal, input/processing, Wordbank cross-sectional) | `datasets_table.R` | `paper/cache/table1_datasets.csv` (built by `paper/build_table1.R`) | `out/datasets_table.tex` |
| σ_r in real units: fold change per SD, 10th–90th percentile range in words per waking hour at the datasets' reference rate, 90th:10th ratio, implied input share EN / NO (band edge, three anchors, one value beyond) | `sigma_r_tables.R` | `paper/cache/fig_io_imputed_proc.rds` (`panelA$input_ref`, curves) | `out/sigma_r_ranges.tex` |
| σ_r evidence: all-adult tokens/h in Western samples — literature rows with a within-sample SD (`studies/input_estimation/validation_set.csv`), this paper's four datasets (raw per-child SD), the measurement model's noise-corrected estimate | `sigma_r_tables.R` | same + the validation set | `out/sigma_r_evidence.tex` |

```bash
Rscript paper/tables/datasets_table.R
Rscript paper/tables/sigma_r_tables.R
```

The three English Wordbank longitudinal sources print as "[year TODO]" until the
citation strings in `build_table1.R` are fixed (the cache carries "20XX").
