#!/usr/bin/env bash
## Render every figure in paper/figures/ to paper/figures/out/.
## R figures read the committed caches only (paper/cache/ + study caches);
## fig_data.R additionally needs data/ + the fitted io-proc bundle locally.
## The TikZ graphical model compiles with pdflatex (TinyTeX is found even
## when it is not on PATH). Each figure is attempted independently; failures
## are listed at the end and the exit code is non-zero if any failed.
##
##   bash paper/figures/render_all.sh            # everything
##   bash paper/figures/render_all.sh ioproc     # one figure (name = script suffix)
set -uo pipefail
cd "$(dirname "$0")/../.."          # repo root (here() needs it)
ONLY="${1:-}"
FAILED=()
run() {
  local name="$1"; shift
  if [ -n "$ONLY" ] && [ "$ONLY" != "$name" ]; then return 0; fi
  echo "== $name"
  if ! "$@"; then FAILED+=("$name"); echo "   !! $name FAILED"; fi
}
PDFLATEX="$(command -v pdflatex || ls "$HOME"/Library/TinyTeX/bin/*/pdflatex "$HOME"/.TinyTeX/bin/*/pdflatex /Library/TeX/texbin/pdflatex 2>/dev/null | head -1 || true)"

tikz_fig() {   # tikz_fig <tex basename without .tex> <out name>
  local src="$1" out="$2"
  [ -n "$PDFLATEX" ] || { echo "   pdflatex not found (install TinyTeX: quarto install tinytex)"; return 1; }
  ( cd paper/figures && "$PDFLATEX" -interaction=batchmode -halt-on-error "$src.tex" > /dev/null \
      && mv "$src.pdf" "out/$out.pdf" && rm -f "$src.aux" "$src.log" ) || return 1
  if command -v pdftoppm > /dev/null; then
    pdftoppm -png -r 200 -singlefile "paper/figures/out/$out.pdf" "paper/figures/out/$out" || true
  fi
  echo "wrote paper/figures/out/$out.pdf"
}

run graphical_model tikz_fig fig_graphical_model graphical_model
run demographics    Rscript paper/figures/fig_demographics.R
run data            Rscript paper/figures/fig_data.R
run ioproc          Rscript paper/figures/fig_ioproc.R
run imputed_meta    Rscript paper/figures/fig_imputed_meta.R

if [ ${#FAILED[@]} -gt 0 ]; then echo "FAILED: ${FAILED[*]}"; exit 1; fi
echo "done -> paper/figures/out/"
