# studies/ — analysis provenance map

One entry per analysis in the paper, mapping **paper element → code → fits/cache → figure**.
The goal is provenance: for any claim, this table points to the scripts that produced it.

Two kinds of study live here:
- **Self-contained dirs** (their own code + cache): `glmer_ladder/`, `input_estimation/`,
  `cross_sectional_demographics/`. Moved here from the repo root in the 2026-06 reorg
  (see [`/MOVES.md`](../MOVES.md)).
- **Provenance stubs** (`proc_dp/`, `io_pooled/`, `longitudinal/`, `llm/`): the code lives
  in the shared engine `model/scripts/` + `model/stan/` and is *indexed* here rather than
  moved, to avoid breaking `source()` paths to `model/R/config.R`. Each stub README lists
  the exact scripts, Stan models, cluster jobs, and fit outputs for that study.

## Map

| Paper element | Study | Code | Fits / cache | Figure cache |
|---|---|---|---|---|
| Fig 2 (model ladder), Table 2 (AIC/BIC) | `glmer_ladder/` | `studies/glmer_ladder/0*.R` | `fits/glmer_ladder/data_*.rds`, `fits/glmer_ladder/sim_cache.rds` | `paper/cache/fig2_glmer_ladder.rds` |
| Demographics figure (sex, mat-ed) | `cross_sectional_demographics/` | `studies/cross_sectional_demographics/00_build.R`, `composite_figure.R` | `…/cache/fits.rds`, `…/cache/frames/*.rds` | (built inline) |
| Demographics — longitudinal arm (BLUPs) | `longitudinal/` (stub) | `model/scripts/analyze_longitudinal.R` | `fits/long_*` | `paper/cache/blups_demographics.rds` |
| Fig 3 A–D (input share, input fan) | `io_pooled/` (stub) | `model/scripts/prepare_io_pooled.R`, `fit_io_pooled_widedelta.R` | `fits/io_pooled*`, `fits/io_{am2018,fmw2013}_*` | `paper/cache/fig3_input.rds` |
| Fig 3 E (processing fan) | `proc_dp/` (stub) | `model/scripts/prepare_proc_dp_bundle.R`, `fit_proc_dp.R` | `fits/proc_dp_*` | `paper/cache/fig3_input.rds` |
| Fig 4 (exposures-to-learn) | `longitudinal/` (stub) | `model/scripts/exposure_to_learn.R` | `fits/long_*` | `paper/cache/fig4_exposure.rds` |
| Fig 5 (LLM acceleration) | `llm/` (stub) | `model/scripts/feng_eval/` | `fits/feng_eval/`, `fits/glmer_mbest_*` | `paper/cache/fig6_llm_slopes.rds` |
| Table 1 (dataset inventory) | (all) | `paper/build_table1.R` | reads all bundles above | `paper/cache/table1_datasets.csv` |
| σ_r literature band | `input_estimation/` | `studies/input_estimation/*.R` | `studies/input_estimation/validation_set.csv` | — |

The narrative history behind each of these lives in [`/journal/experiments.md`](../journal/experiments.md)
(child models + input/processing) and [`/journal/experiments_llm.md`](../journal/experiments_llm.md) (LLM work).
