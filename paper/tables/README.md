# paper/tables — LaTeX tables for the Overleaf manuscript / SI

One script per table, reading committed caches only; the `.tex` output in `out/`
is committed and uploaded to Overleaf (`\input{}` it inside a table environment;
needs `booktabs`).

| table | script | reads | output |
|---|---|---|---|
| Datasets and analysis samples (SI Table S1: Wordbank longitudinal, input/processing, Wordbank cross-sectional) | `datasets_table.R` | `paper/cache/table1_datasets.csv` (built by `paper/build_table1.R`) | `out/datasets_table.tex` |

```bash
Rscript paper/tables/datasets_table.R
```

The three English Wordbank longitudinal sources print as "[year TODO]" until the
citation strings in `build_table1.R` are fixed (the cache carries "20XX").
