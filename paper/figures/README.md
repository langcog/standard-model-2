# paper/figures — standalone figure pipeline (Overleaf era)

The manuscript now lives on Overleaf (PNAS two-column, with a collaborator), so
figures are built here, one script per figure, and uploaded as vector PDFs.
Nothing is numbered: order is decided in the manuscript.

```bash
bash paper/figures/render_all.sh            # all five -> out/*.pdf (+ .png previews)
bash paper/figures/render_all.sh ioproc     # just one (name = script suffix)
```

Every `fig_*.R` sources `_common.R` (PNAS widths, the palette, `theme_fig()`,
`tag_panels()`, `save_fig()`) and reads **committed caches only** — `paper/cache/*.rds`
(built by `paper/build_cache.R`, `paper/build_fig_io_cache.R`) and
`studies/cross_sectional_demographics/cache/fits.rds` — so anyone can rebuild without
the fits. The one exception is `fig_data.R`, which plots raw data and needs `data/`
plus the fitted io-proc bundle. R deps beyond the tidyverse: `patchwork`, `ggh4x`
(tinted facet strips; CRAN). The TikZ figure needs pdflatex (TinyTeX is found).

| figure | script | size (in) | notes |
|---|---|---|---|
| Graphical model: the three measurement models (vocabulary, input, processing) and the structural links into ξ / κ | `fig_graphical_model.tex` (TikZ) | 1 col, 3.43 × 3.39 | nodes tinted with the palette; natural width = column width so type prints at size |
| Demographics: sex and maternal-ed effects on efficiency vs acceleration, per language, cross-sectional + longitudinal, meta panel | `fig_demographics.R` → `composite_figure.R` | 2 col, 7.01 × 5.6 | (A) sex forest tall on the left, (B) maternal ed over (C) meta on the right; (C) has its own x-scale, shared by its two components (the forests' acceleration axis is set by the widest per-language CIs); languages sorted by efficiency effect (a `geom_blank` scale-training bug had left them alphabetical); longitudinal arm on the clean BLUPs (experiments.md #44) |
| Data entering the io-proc analysis: CDI trajectories, LWL RT, observed input | `fig_data.R` → `data_check_io_panels.R` | 2 col, 7.01 × 3.7 | dataset small multiples (rows = streams, columns = datasets, empty cell = not collected), stream colours |
| io-proc: (A) ±2 SD schematic from the joint fit, (B) coefficients — each channel fit alone vs jointly | `fig_ioproc.R` | 2 col, 7.01 × 3.4 | (A) one average child, direct-labelled, θ and words stacked; (B) paired separate/joint (experiments.md #43) under tinted trait strips; each panel has its own legend |
| Imputed population input share vs the meta-analytic 4–7% range (σ_r sweep + six anchors) | `fig_imputed_meta.R` | 1 col, 3.43 × 3.08 | input hue, direct-labelled English / Norwegian (experiments.md #45) |

## Style conventions (every figure)

- **Panels.** Each panel is its own box (light panel border, no axes shared across
  panels) with a bold letter at its top-left — `tag_panels()` is patchwork's
  `tag_levels = "A"`, so letters align across panels by construction. A short bold
  title sits inside the panel; explanation goes in the caption. Single-panel figures
  carry no letter.
- **Type** (printed sizes; figures are saved at their final width): axis text 7 pt,
  axis titles / strips / panel titles 8 pt, panel letters 10 pt bold, in-plot
  annotations 7 pt (`TXT`). TikZ: `\footnotesize` (8 pt) node labels, `\scriptsize`
  (7 pt) small nodes, edge and plate labels; Helvetica for text.
- **Palette — one meaning per hue, everywhere** (`_common.R`; `fig_graphical_model.tex`
  repeats the hex values). Validated with the dataviz palette checker: every pair
  passes the deutan / tritan / normal-vision separation checks and 3:1 contrast on white.

  | hue | hex | meaning | appears in |
  |---|---|---|---|
  | blue | `#0072B2` | efficiency (ξ) | ξ box in the graphical model; efficiency strips + points (demographics, io-proc B) |
  | vermillion | `#D55E00` | acceleration (κ) | κ box; acceleration strips + points |
  | green | `#009E73` | input (d_i, r_iv) | graphical model; data row C; io-proc A lines and B points; imputed curves |
  | purple | `#AA3377` | processing (w_i, q_in) | graphical model; data row B; io-proc A lines and B points |
  | dark grey | `#404040` | vocabulary / the outcome | y_ijt; data row A; the average child in io-proc A |

  Contrasts that are not traits or channels are carried by **shape**, not colour:
  cross-sectional = filled circle vs longitudinal = open diamond (demographics);
  separate = open vs joint = filled circle (io-proc B); ± perturbations = solid vs
  dashed (io-proc A). Curves are direct-labelled where a legend would only repeat
  the colours (io-proc A, imputed).

PNAS sizes used: single column 8.7 cm, 1.5 column 11.4 cm, double 17.8 cm, max height
22.5 cm.

**Rebuild chain after a fit changes:** `build_cache.R` / `build_fig_io_cache.R` →
`render_all.sh` → upload the changed `out/*.pdf` to Overleaf. `out/*.pdf` are
committed; `out/*.png` previews are gitignored.

The pre-Overleaf Quarto manuscript (`paper/input_paper.qmd`) still contains the
original figure chunks; they are superseded by these scripts.
