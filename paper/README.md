# Manuscript build

## Render

```bash
cd paper
quarto render standard_model.qmd
```

Produces `standard_model.pdf` + `standard_model.tex` (we keep the
.tex committed for debugging — set by `keep-tex: true` in the YAML).

## Build the figure caches first

Most figure chunks load small summary RDS files from `paper/cache/`.
To regenerate from scratch (after a fit changes, or on a fresh clone):

```bash
Rscript paper/build_cache.R
```

This produces:

- `cache/table1_datasets.rds` — Table 1 contents
- `cache/fig2_glmer_ladder.rds` — quantile predictions + empirical points
- `cache/blups_demographics.rds` — per-kid BLUPs joined with Wordbank demos
- `cache/fig5_io_summary.rds` — IO-pooled posterior medians for Fig 5

These are <30 KB total and *are* committed (cheap to keep in git).

The full Bayesian fits underneath are NOT in git — they live in
`fits/` (gitignored) or Sherlock scratch. Re-running `build_cache.R`
requires the fits to be present locally (see `outputs/experiments.md`).

## Rendering troubleshooting

If the PDF comes out single-column / abstract bleeds wide / footnote
overruns:

1. The PNAS template (`pnasresearcharticle.sty`) sets the
   `shortarticle` boolean which controls the title block + column
   layout. If it doesn't load, you fall back to a wide single-column
   title page.
2. Confirm `paper/pnasresearcharticle.sty` exists after a render
   (Quarto copies it from `_extensions/christopherkenny/pnas/`
   via `format-resources`).
3. Inspect `standard_model.log` for any `Package ... not found` or
   `File ... not found` errors — TinyTeX should auto-fetch but
   sometimes silently skips.
4. Try `quarto render --log-level info standard_model.qmd 2>&1 | tee render.log`
   for a much more verbose trace.
5. Re-install the extension if all else fails:
   ```bash
   cd paper
   quarto remove christopherkenny/pnas
   quarto add christopherkenny/pnas
   ```

The expected first-page layout: narrow title (left half),
single-column author block + abstract + keywords on the left, then
two-column body starting below.

## What's committed

- `standard_model.qmd` — manuscript source (author-written prose +
  agent-written chunks; see `~/.claude/skills/rmarkdown-manuscript-collab.md`
  for the collaboration protocol)
- `_helpers.R` — shared palettes + Fig 1 schematic helper
- `build_cache.R` — single-entry cache rebuilder
- `cache/*.rds` — the tiny summaries the chunks load
- `_extensions/christopherkenny/pnas/` — PNAS Quarto extension
- `standard_model.bib` — bibliography (entries from `bibs/` only, no
  agent-generated citations)
- `widetext.sty`, `pnas-new.cls`, `pnasresearcharticle.sty`, etc. —
  PNAS support files (kept locally for portability with TinyTeX)
