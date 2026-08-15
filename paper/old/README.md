# Retired papers (archived 2026-08-15)

Nothing in here is built anymore. The live manuscript is `paper/input_paper.qmd`;
this folder preserves the two earlier papers that shared this directory, plus
their format scaffolding.

## `short_paper/` — the acceleration paper

"Children, but not language models, show accelerating returns in word learning"
(`standard_model_short.qmd`, Science format). **This project moved to its own
repo, `~/Projects/acceleration`, which is now canonical** — the copies here were
byte-identical to that repo's at the time of retirement. Includes its SI
(`supplemental.qmd`, the bayes_long M0–M3 one), the Science supplement
(`science/`), the Science submission template, its cache builders
(`build_cache_short.R`, `build_cache_si_gamlss.R`, `build_cache_si_settings.R`),
`_setup_shared.R`, and the Fig 1 graphical-model `mockups/`.

One exception: `si_clippings.qmd` (SI prose/chunk clippings, acceleration
content) was **not** copied to the acceleration repo before the split — this is
the only copy.

## `long_paper/` — the original combined paper

"Children's early vocabulary growth is a process of accelerating accumulation"
(`standard_model.qmd`, PNAS format) — the single long paper from which both the
acceleration paper and the input paper were carved. Includes its rendered
pdf/tex/html, its SI renders (`supplemental.pdf`/`.tex`), and the journal format
scaffolding only it used: the PNAS Quarto extension plus support files
(`_extensions/christopherkenny/pnas`, `pnas-new.*`, `pnasresearcharticle.sty`,
`jabbrv*`, `widetext.sty`) and the unused ACS extension
(`_extensions/quarto-journals/acs`, `achemso.bst`, `acs-standard_model.bib`).

## Loose

`Standard Model 2.0 BUCLD abstract.pdf` — conference abstract for the program.
