# paper/figures — standalone figure pipeline (Overleaf era)

The manuscript now lives on Overleaf (PNAS two-column, with a collaborator), so
figures are built here, one script per figure, and uploaded as vector PDFs.
Nothing is numbered: order is decided in the manuscript.

```bash
bash paper/figures/render_all.sh            # all five -> out/*.pdf (+ .png previews)
bash paper/figures/render_all.sh ioproc     # just one (name = script suffix)
```

Every `fig_*.R` sources `_common.R` (PNAS widths, palettes, `save_fig()`) and reads
**committed caches only** — `paper/cache/*.rds` (built by `paper/build_cache.R`,
`paper/build_fig_io_cache.R`) and `studies/cross_sectional_demographics/cache/fits.rds`
— so anyone can rebuild without the fits. The one exception is `fig_data.R`, which
plots raw data and needs `data/` plus the fitted io-proc bundle.

| figure | script | width | status |
|---|---|---|---|
| Graphical model with the three measurement models (vocabulary, input, processing) and the structural links into ξ / κ | `fig_graphical_model.tex` (TikZ) | 1 col | **prototype** — first version with measurement models; the acceleration-paper version lacked them |
| Demographics: sex and maternal-ed effects on efficiency vs acceleration, per language, cross-sectional + longitudinal, meta panel | `fig_demographics.R` → `composite_figure.R` | 2 col, full height | unchanged content; longitudinal arm on the clean BLUPs (experiments.md #44) |
| Data entering the io-proc analysis: CDI trajectories, LWL RT, observed input | `fig_data.R` → `data_check_io_panels.R` | 2 col | unchanged content |
| io-proc: (A) ±2 SD schematic from the joint fit, (B) coefficients — each channel fit alone vs jointly | `fig_ioproc.R` | 2 col | **new panel B** (paired separate/joint; experiments.md #43) |
| Imputed population input share vs the meta-analytic 4–7% range (σ_r sweep + six anchors) | `fig_imputed_meta.R` | 1 col | factored out of the io-proc figure (experiments.md #45) |

PNAS sizes used: single column 8.7 cm, 1.5 column 11.4 cm, double 17.8 cm, max height
22.5 cm; ggplot `base_size` 8 keeps type ≥ 6 pt at print scale.

**Rebuild chain after a fit changes:** `build_cache.R` / `build_fig_io_cache.R` →
`render_all.sh` → upload the changed `out/*.pdf` to Overleaf. `out/*.pdf` are
committed; `out/*.png` previews are gitignored.

The pre-Overleaf Quarto manuscript (`paper/input_paper.qmd`) still contains the
original figure chunks; they are superseded by these scripts.
