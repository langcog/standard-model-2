# standard_model_2

Bayesian accumulator model of early word learning, linking language
input, vocabulary outcomes, and individual differences. Successor to
Kachergis, Marchman, & Frank (2021). See
[`notes/model_explainer.pdf`](notes/model_explainer.pdf) for the model
specification and [`notes/experiments.md`](notes/experiments.md) for a
running log of fits and findings.

## Layout

```
standard_model_2/
├── Makefile                        ← single entry point for all workflows
├── model/
│   ├── stan/                       ← Stan models (cross-sectional + longitudinal)
│   ├── R/                          ← config + helpers + PPC (sourced by scripts)
│   ├── scripts/                    ← driver scripts (numbered in workflow order)
│   ├── fits/                       ← output .rds files (gitignored)
│   └── figs/                       ← output .png figures (gitignored)
├── notes/
│   ├── model_explainer.{tex,pdf}   ← durable model specification
│   ├── experiments.md              ← running log: each fit + backlog
│   ├── model_summaries.md          ← literature review notes
│   └── _archive/                   ← superseded standalone findings files
├── data/                  ← raw external inputs (Sperry, BabyView, etc.)
└── sherlock/                       ← SLURM scripts for remote fits
```

## Local workflow (small fits on your laptop)

```bash
make smoke              # sanity-check everything loads
make recovery           # parameter recovery on simulated data
make data               # build Wordbank subsample Stan data (reads CDI + CHILDES)
make variant NAME=2pl   # fit one cross-sectional variant
make analyze NAME=2pl   # plots + scalar summary
```

Full list of targets: `make` with no argument prints usage.

## Remote workflow (Sherlock — for the bigger fits)

See [`sherlock/README.md`](sherlock/README.md) for step-by-step. One-liner:

```bash
# On Sherlock login node, after one-time setup:
sbatch sherlock/long_fit.slurm long_2pl_slopes_nor
```

Results land in `$SCRATCH/standard_model_2/fits/` and are synced home via
`rsync`. See Sherlock README for the full cycle.

## Getting started (on a fresh clone)

```bash
git clone <repo> standard_model_2
cd standard_model_2

# Local only — install R packages:
Rscript sherlock/setup_R.R          # works both locally and on Sherlock

# Wordbank longitudinal data is pulled by model/scripts/pull_longitudinal.R
# (requires childesr / wordbankr; uses preprocessed bundles when available
# at fits/long_subset_data.rds).
#
# The Sperry / Hart-Risley / Weisleder-Fernald per-recording rate CSV
# lives at data/sperry/hourly_tokens_Sperry_HartRisley.csv.

# Sanity check
make smoke
```

The code auto-detects the project root (by searching for `Makefile` in
the cwd and parents), or respects the env var `STANDARD_MODEL_ROOT`.
Output paths can be redirected with `STANDARD_MODEL_FITS_DIR` and
`STANDARD_MODEL_FIGS_DIR` — used on Sherlock to send outputs to
`$SCRATCH`.
