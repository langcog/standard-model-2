# paper/figs — committed figure assets

TikZ graphical models (rendered to PDF, included by `supplemental.qmd` → `@fig-graphical-models`):
- `graphical_model_accumulator.{tex,pdf}` — (A) the accelerating accumulator
- `graphical_model_ioproc.{tex,pdf}`      — (B) + input & processing measurement models

Regenerate after editing a `.tex`:
```
cd paper/figs && pdflatex graphical_model_accumulator.tex && pdflatex graphical_model_ioproc.tex && rm -f *.aux *.log
```
