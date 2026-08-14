# Science submission targets

Three outputs from one source. Nothing is duplicated between them.

| command | output | what it is |
|---|---|---|
| `quarto render paper/standard_model_short.qmd --to docx` | `standard_model_short.docx` | **main text for submission** — figures embedded, SI excluded |
| `quarto render paper/science/supplement.qmd` | `paper/science/supplement.pdf` | **supplementary materials for submission** — separate file, S-numbered |
| `quarto render paper/standard_model_short.qmd --to pdf` | `standard_model_short.pdf` | combined copy for internal review (unchanged) |

## Why it is arranged this way

Science asks for "a single, complete file that includes all figures and tables in Word's
.docx format", and says plainly: *"If you are using LaTeX, please convert your paper into
.docx format."* The supplementary materials go as "a single separate file in .docx or PDF
format."

The SI stays **PDF**, which is permitted, and that is what makes this cheap: the SI carries
thirteen kableExtra tables using `pack_rows`, `add_header_above` and `landscape`, none of
which survive a Word conversion. The main text carries **no tables at all**, so it converts
with nothing lost.

## Things that will bite if they are changed

- **`dev` must follow the format.** The setup chunk sets `cairo_pdf` for LaTeX and
  `ragg_png` otherwise. A docx built with `cairo_pdf` embeds three `.pdf` files that Word
  renders as blank placeholders — the render succeeds and the figures are simply missing.
- **The SI include is wrapped in `.content-visible when-format="pdf"`.** Removing that puts
  the whole supplement inside the submission docx.
- **The S in `fig. S1` lives in the LaTeX counter,** not in `fig-title`. Quarto joins label
  and number with a hardcoded space, and for PDF that join happens inside LaTeX's caption
  machinery where a Lua filter on the Pandoc AST cannot reach it. Hence
  `\renewcommand{\thefigure}{S\arabic{figure}}` in `supplement.qmd`.
- **`paper/_setup_shared.R` is sourced by both documents.** It was extracted from the main
  .qmd's setup chunk because the SI quotes values (`n_excl`, `pct_excl`, `qc_mar`,
  `qc_nor`) that used to exist only inside that chunk — which is why rendering
  `supplemental.qmd` on its own always failed with `object 'n_excl' not found`. Do not
  duplicate it back into either file.

## Still outstanding

Checked against
<https://www.science.org/content/page/instructions-preparing-initial-manuscript>:

- Main text is **3,047 words** against the 3,000 limit for a standard Research Article.
- **84 main-text references** against ~50.
- Required section keywords in order (`Title:`, `Authors:`, `Affiliations:`, `Abstract:`,
  `Main Text:`, `References and Notes`, `Acknowledgments:`, `List of Supplementary
  Materials:`).
- Acknowledgments subsections: Funding, Author contributions, Competing interests,
  Data/code/materials availability.
- Figure captions must open with a short **bold title**, max 200 words.
- **Red and green together is prohibited** ("this creates problems for individuals with
  color-deficient vision"). Fig 1B's percentile palette runs blue → green → red; Fig 2 uses
  Dark2, which pairs green with pink. Both need recolouring.
