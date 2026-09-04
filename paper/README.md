# Manuscript build

The paper in this folder is the **input paper**:

> *Dissociating efficiency and acceleration in the growth of children's early
> vocabulary* (Frank & Marchman) — `input_paper.qmd`

(The acceleration paper — children vs. language models — lives in its own repo,
`~/Projects/acceleration`. Everything retired at the split is in `old/`; see
`old/README.md`.)

## Render

```bash
quarto render paper/input_paper.qmd
```

Produces `input_paper.pdf` + `input_paper.tex` (kept out of git; `keep-tex:
true` so the .tex is there for debugging). Plain `article` class — no journal
extension needed. Chunks resolve paths with `here()`, so rendering works from
the repo root or from `paper/`.

The SI is `full_supplemental.qmd` (no YAML header of its own; its `figs/`
image paths are relative to this directory).

## Build the figure caches first

Figure chunks load small summary RDS/CSV files from `paper/cache/`, which are
committed so the manuscript builds without cluster access. To regenerate after
a fit changes (requires local fits in `fits/`; see `journal/`):

| script | rebuilds |
|---|---|
| `build_cache.R` | core caches (glmer ladder, BLUPs+demographics, IO summary, …) |
| `build_input_cache.R` | `fig3_input.rds` — input-experiments triptych (Fig 3) |
| `build_fig_io_cache.R` | `fig_io_imputed_proc.rds` — 2-panel IO figure |
| `build_table1.R` | `table1_datasets.csv` — dataset characteristics (Table 1) |
| `build_si_io_data_table.R` | `si_io_data_table.csv` — SI io-proc data table |

## What's committed

- `input_paper.qmd`, `full_supplemental.qmd` — manuscript + SI source
  (author-written prose + agent-written chunks; agents don't touch prose)
- `_helpers.R` — shared palettes + schematic helpers
- `build_*.R` — cache rebuilders (table above)
- `cache/` — the small summaries the chunks load
- `standard_model.bib` — bibliography (entries verified against
  Crossref/arXiv, never written from memory)
- `figs/` — static figures (graphical models, TikZ sources)
- `old/` — retired papers (tracked history preserved; new render artifacts
  in there are gitignored)

## 2026-09: manuscript moved to Overleaf; figures built standalone

The manuscript is now on Overleaf (PNAS two-column, with a collaborator), so
`input_paper.qmd` is frozen as the pre-Overleaf source. Figures are built one
script per figure in **`paper/figures/`** (see its README) and uploaded as
vector PDFs from `paper/figures/out/`:

```bash
bash paper/figures/render_all.sh
```

The cache chain is unchanged (`build_cache.R`, `build_fig_io_cache.R` → `cache/`);
the figure scripts read only those committed caches.
