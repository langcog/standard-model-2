# Repository reorganization — 2026-06-10

The repo was reorganized for legibility and provenance. This file is the
old → new path map. **If a script breaks on a missing path, look it up here.**

Policy: paper-render-critical paths were updated in place; the long tail of
cluster/analysis scripts was **not** rewritten (they often assume a different
cluster-side layout anyway). Fix those lazily — when you next touch a script,
update its paths and delete its row here.

## Directory moves

| Old | New | Notes |
|---|---|---|
| `gcp/` | `cluster/gcp/` | compute helpers |
| `sherlock/` | `cluster/sherlock/` | compute helpers |
| `glmer_ladder/` | `studies/glmer_ladder/` | study code (self-contained) |
| `input_estimation/` | `studies/input_estimation/` | study code (self-contained) |
| `cross_sectional_demographics/` | `studies/cross_sectional_demographics/` | **PENDING** — moved after the live uncap refit finishes |
| `outputs/figs/` | `figs/` | merged into the existing figs/; `io/` merged |
| `outputs/slides/` | `reports/slides/` | |
| `outputs/marlowe/` | `reports/marlowe/` | grant proposal |
| `outputs/{model_explainer,chang_bergen_derivation,model_comparison_unified}.{tex,pdf}` | `reports/` | standalone docs |
| `outputs/feng_eval/` | `fits/feng_eval/` | LLM eval data |
| `outputs/glmer_ladder/{sim_cache.rds,*.csv}` | `fits/glmer_ladder/` | study outputs |
| `outputs/experiments.md` | `journal/experiments.md` | canonical log |
| `outputs/experiments_llms.md` | `journal/experiments_llm.md` | LLM log (renamed) |
| `outputs/PROVENANCE.md` | `journal/PROVENANCE.md` | |
| `outputs/_archive/` | `journal/archive/` | old findings |
| `outputs/*.md` (notes/plans) | `journal/notes/` | |
| `outputs/{param,input_rate,quality_variation}_table.*`, `predictors_alpha_zeta*.csv`, etc. | `journal/results/` | result tables |
| `run_gamma_weekend.sh` | `cluster/run_gamma_weekend.sh` | local-laptop fit launcher |
| `outputs/` | *(deleted)* | fully dissolved |

`fits/`, `data/`, `model/`, `paper/`, `papers/` are unchanged. `model/scripts/`
remains the shared engine (the proc_dp / io_pooled / longitudinal / llm studies
are *indexed* from `studies/<name>/README.md`, not moved).

## Paths already fixed
- `paper/build_cache.R` — `outputs/glmer_ladder/sim_cache.rds` → `fits/…`; `input_estimation/…` → `studies/input_estimation/…`
- `Makefile` — `outputs/figs` → `figs`; `outputs/model_explainer` → `reports/model_explainer`
- `.gitignore` — figs/feng_eval/manifest/quarto-artifact rules repointed

## Known-stale references (fix when next touched)
- `cluster/sherlock/*.slurm`, `cluster/sherlock/glmer_ladder_submit.sh` — invoke
  `glmer_ladder/02_fit_one.R` etc.; these run cluster-side under `$HOME/standard_model_2`
  with a layout that may differ from this repo. Update `glmer_ladder/` → `studies/glmer_ladder/`
  if you re-sync the tree.
- `model/scripts/*` — various scripts write to `outputs/...`; repoint to `figs/`,
  `fits/`, or `journal/results/` as appropriate when you re-run them.
- `paper/_helpers.R` — a comment references `glmer_ladder/04b_plot.R` (now `studies/…`).
