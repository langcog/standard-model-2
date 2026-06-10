#!/usr/bin/env python3
"""
Preprocess standard_model.qmd into a PPTX-buildable variant.

Steps:
  1. Compile all 6 TikZ blocks to standalone PNGs (via pdflatex + pdftocairo).
  2. Convert irt.pdf -> irt.png (PPTX can't embed PDFs).
  3. Expand custom LaTeX commands (\\pia, \\sigxi, etc.) inline.
  4. Strip LaTeX layout commands that PPTX ignores (\\vspace, \\Large, etc.).
  5. Convert \\includegraphics{...} to markdown ![](...) form.
  6. Replace the YAML header with a PPTX-targeted version.
  7. Render: quarto render standard_model_pptx.qmd --to pptx

Run from this directory: python3 build_pptx.py
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

SLIDES_DIR = Path(__file__).resolve().parent
SRC        = SLIDES_DIR / "standard_model.qmd"
DST        = SLIDES_DIR / "standard_model_pptx.qmd"
TIKZ_DIR   = SLIDES_DIR / "tikz_pngs"
FIGS_DIR   = SLIDES_DIR.parent / "figs"

# Resource paths the original deck uses via \graphicspath in _preamble.tex.
# Quarto pptx finds images via resource-path; we hand it the same list.
RESOURCE_PATHS = [
    ".",
    "tikz_pngs",
    "../figs",
    "../figs/longitudinal",
    "../figs/io",
    "../figs/proc",
    "../figs/schematic",
]

# Custom LaTeX command -> math expansion (matches _preamble.tex defs).
# Order matters: longer names first so \sigxi doesn't match the \sig prefix
# of \sigxi/\siga/etc.
CMD_REPLACEMENTS = [
    (r"\siglam",  r"\sigma_{\lambda}"),
    (r"\kappap",  r"\kappa_{\text{pop}}"),
    (r"\sigxi",   r"\sigma_{\xi}"),
    (r"\sigz",    r"\sigma_{\zeta}"),
    (r"\sigr",    r"\sigma_{r}"),
    (r"\sigs",    r"\sigma_{s}"),
    (r"\siga",    r"\sigma_{\alpha}"),
    (r"\pia",     r"\pi_{\alpha}"),
]

# Standalone preamble for compiling each TikZ block to its own PDF.
# Mirrors the libraries / colors / custom commands loaded in _preamble.tex.
TIKZ_STANDALONE_PREAMBLE = r"""\documentclass[border=4pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{matrix, positioning, arrows.meta}
\usepackage{xcolor}
\usepackage{amsmath, amssymb}
\definecolor{accent}{RGB}{195,30,55}
\definecolor{cool}{RGB}{31,120,180}
\newcommand{\pia}{\pi_{\alpha}}
\newcommand{\sigxi}{\sigma_{\xi}}
\newcommand{\siga}{\sigma_{\alpha}}
\newcommand{\sigr}{\sigma_{r}}
\newcommand{\sigz}{\sigma_{\zeta}}
\newcommand{\siglam}{\sigma_{\lambda}}
\newcommand{\sigs}{\sigma_{s}}
\newcommand{\kappap}{\kappa_{\text{pop}}}
\begin{document}
"""
TIKZ_STANDALONE_END = "\n\\end{document}\n"


def run(cmd, **kw):
    """Run a subprocess and surface its stderr if it fails."""
    res = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if res.returncode != 0:
        print(f"FAILED: {' '.join(map(str, cmd))}", file=sys.stderr)
        print(res.stderr[-2000:], file=sys.stderr)
        sys.exit(res.returncode)
    return res


def compile_tikz_block(idx, src):
    """Compile one TikZ block to PNG. Returns the PNG path relative to slides dir."""
    tex_path = TIKZ_DIR / f"tikz_{idx}.tex"
    pdf_path = TIKZ_DIR / f"tikz_{idx}.pdf"
    png_stem = TIKZ_DIR / f"tikz_{idx}"   # pdftocairo appends .png

    tex_path.write_text(TIKZ_STANDALONE_PREAMBLE + src + TIKZ_STANDALONE_END)
    run(["pdflatex", "-interaction=nonstopmode",
         f"-output-directory={TIKZ_DIR}", str(tex_path)])
    run(["pdftocairo", "-png", "-r", "240", "-singlefile",
         str(pdf_path), str(png_stem)])
    return f"tikz_pngs/tikz_{idx}.png"


def convert_includegraphics(match):
    """\\includegraphics[width=0.95\\linewidth]{foo.png} -> ![](foo.png){width=95%}"""
    opts = match.group(1) or ""
    filename = match.group(2)
    w = re.search(r"width\s*=\s*([0-9.]+)\s*\\linewidth", opts)
    if w:
        pct = int(round(float(w.group(1)) * 100))
        return f"![]({filename}){{width={pct}%}}"
    return f"![]({filename})"


def preprocess(text):
    # 1. Replace TikZ blocks with image references.
    tikz_pat = re.compile(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", re.DOTALL)
    blocks = list(tikz_pat.finditer(text))
    print(f"Found {len(blocks)} TikZ blocks; compiling to PNG...")
    pieces = []
    last = 0
    for i, m in enumerate(blocks, 1):
        png_rel = compile_tikz_block(i, m.group(0))
        pieces.append(text[last:m.start()])
        pieces.append(f"![]({png_rel}){{width=80%}}")
        last = m.end()
        print(f"  block {i}: {png_rel}")
    pieces.append(text[last:])
    text = "".join(pieces)

    # 2. Expand custom commands. Word-boundary on the right so \sigxi
    #    doesn't get partially matched against \sig.
    for cmd, repl in CMD_REPLACEMENTS:
        # Match \cmd followed by non-letter or end-of-string (LaTeX command boundary).
        pattern = re.escape(cmd) + r"(?![a-zA-Z])"
        text = re.sub(pattern, lambda m, r=repl: r, text)

    # 3. Strip layout commands.
    text = re.sub(r"\\vspace\{[^}]*\}", "", text)
    text = re.sub(r"\\vspace\*\{[^}]*\}", "", text)
    text = re.sub(r"\\hspace\{[^}]*\}", "", text)
    text = re.sub(r"\\centering(?![a-zA-Z])", "", text)
    text = re.sub(r"\\quad(?![a-zA-Z])", " ", text)
    text = re.sub(r"\\qquad(?![a-zA-Z])", "  ", text)
    # Font-size commands -- PPTX uses its own; strip all of them.
    SIZE_CMDS = ("tiny", "scriptsize", "footnotesize", "small", "normalsize",
                 "large", "Large", "LARGE", "huge", "Huge")
    for cmd in SIZE_CMDS:
        text = re.sub(r"\\" + cmd + r"(?![a-zA-Z])", "", text)
    # \color{...} -- strip with its braced argument. Doesn't handle
    # \color{...}{...} (two-arg form) since the deck only uses the
    # one-arg form for inline color shifts that PPTX ignores anyway.
    text = re.sub(r"\\color\{[^}]*\}", "", text)
    # Booktabs rules in raw tabular blocks -- pandoc handles plain
    # \hline but strips booktabs commands as unknown. Replace with
    # \hline so the table borders survive.
    text = re.sub(r"\\(?:toprule|midrule|bottomrule)\b", r"\\hline", text)

    # 4. Strip \centerline{...} - keep the inner content.
    #    Has to allow nested braces (one level deep is enough for our deck).
    text = re.sub(r"\\centerline\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}",
                  r"\1", text)

    # 5. Convert \includegraphics to markdown.
    text = re.sub(r"\\includegraphics(?:\[([^\]]*)\])?\{([^}]+)\}",
                  convert_includegraphics, text)

    # Convert irt.pdf reference to irt.png (we'll have built the PNG above).
    text = text.replace("irt.pdf", "irt.png")

    # 6. Replace YAML header.
    yaml_pat = re.compile(r"^---\n.*?^---\n", re.DOTALL | re.MULTILINE)
    new_yaml = (
        "---\n"
        'title: "Acceleration and variability in vocabulary growth for young '
        'children and language models"\n'
        'author: "Michael C. Frank"\n'
        'institute: "Stanford University"\n'
        'date: "2026-06-01"\n'
        'date-format: "YYYY"\n'
        "format:\n"
        "  pptx:\n"
        "    slide-level: 2\n"
        "resource-path:\n"
    )
    for rp in RESOURCE_PATHS:
        new_yaml += f"  - {rp}\n"
    new_yaml += "---\n"
    text = yaml_pat.sub(new_yaml, text, count=1)

    return text


def ensure_irt_png():
    """Convert outputs/slides/irt.pdf -> outputs/slides/irt.png once."""
    pdf = SLIDES_DIR / "irt.pdf"
    png = SLIDES_DIR / "irt.png"
    if pdf.exists() and not png.exists():
        print("Converting irt.pdf -> irt.png...")
        run(["pdftocairo", "-png", "-r", "240", "-singlefile",
             str(pdf), str(SLIDES_DIR / "irt")])


def main():
    TIKZ_DIR.mkdir(exist_ok=True)
    ensure_irt_png()
    text = SRC.read_text()
    out = preprocess(text)
    DST.write_text(out)
    print(f"Wrote {DST}")
    print("\nRender with:  cd outputs/slides && quarto render standard_model_pptx.qmd --to pptx")


if __name__ == "__main__":
    main()
