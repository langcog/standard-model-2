# standard_model_2

Bayesian accumulator model of early word learning, linking language
input, vocabulary outcomes, and individual differences. Successor to
Kachergis, Marchman, & Frank (2021). See
[`reports/model_explainer.pdf`](reports/model_explainer.pdf) for the model
specification and [`journal/experiments.md`](journal/experiments.md) for the
running log of fits and findings (LLM arc: [`journal/experiments_llm.md`](journal/experiments_llm.md)).

## Layout

```
standard_model_2/
├── Makefile                ← local build targets (make with no args prints usage)
├── MOVES.md                ← old→new path map from the 2026-06 reorg
├── paper/                  ← the manuscript (.qmd, build_*.R, cache/)
├── model/                  ← shared engine
│   ├── stan/               ← Stan models
│   ├── R/                  ← config + helpers (sourced by scripts)
│   └── scripts/            ← driver scripts (the studies below index into here)
├── studies/                ← one analysis per subdir + provenance map (README.md)
│   ├── glmer_ladder/             ← Fig 2, Table 2 (model ladder)
│   ├── cross_sectional_demographics/  ← demographics figure
│   ├── input_estimation/         ← σ_r literature band
│   └── {proc_dp,io_pooled,longitudinal,llm}/  ← provenance stubs → model/scripts
├── cluster/                ← compute helpers
│   ├── gcp/                ← Google Cloud VM launchers
│   └── sherlock/           ← Stanford Sherlock SLURM jobs + extractors
├── data/                   ← raw external inputs (Sperry, BabyView, peekbank, …)
├── fits/                   ← model fit outputs (heavy .rds gitignored; summaries/ tracked)
├── figs/                   ← figures (PNGs gitignored) + their source CSVs
├── reports/                ← standalone docs: explainer, derivations, slides, proposal
├── journal/                ← project history (the system of record)
│   ├── experiments.md      ← numbered log of every fit + finding + backlog
│   ├── experiments_llm.md  ← the LLM / GPU arc
│   ├── PROVENANCE.md       ← per-asset provenance for the slide deck
│   └── notes/  results/  archive/
└── papers/                 ← literature PDFs
```

For per-claim provenance (which scripts/fits/figures back each paper element),
see [`studies/README.md`](studies/README.md).

## Local workflow (small fits on your laptop)

```bash
make smoke              # sanity-check everything loads
make recovery           # parameter recovery on simulated data
make data               # build Wordbank subsample Stan data (reads CDI + CHILDES)
make variant NAME=2pl   # fit one cross-sectional variant
make analyze NAME=2pl   # plots + scalar summary
```

Full list of targets: `make` with no argument prints usage.

## Remote workflow (Sherlock / GCP — for the bigger fits)

See [`cluster/sherlock/README.md`](cluster/sherlock/README.md) for step-by-step. One-liner:

```bash
# On Sherlock login node, after one-time setup:
sbatch cluster/sherlock/long_fit.slurm long_2pl_slopes_nor
```

Results land in `$SCRATCH/standard_model_2/fits/` and are synced home via
`rsync`. GCP launchers live in [`cluster/gcp/`](cluster/gcp/).

## Getting started (on a fresh clone)

```bash
git clone <repo> standard_model_2
cd standard_model_2

# Install R packages (works locally and on Sherlock):
Rscript cluster/sherlock/setup_R.R

# Wordbank longitudinal data is pulled by model/scripts/pull_longitudinal.R
# (requires childesr / wordbankr; uses preprocessed bundles when available
# at fits/long_subset_data.rds). The Sperry / Hart-Risley / Weisleder-Fernald
# per-recording rate CSV lives at data/sperry/hourly_tokens_Sperry_HartRisley.csv.

make smoke
```

The code auto-detects the project root (by searching for `Makefile` in
the cwd and parents), or respects the env var `STANDARD_MODEL_ROOT`.
Output paths can be redirected with `STANDARD_MODEL_FITS_DIR` and
`STANDARD_MODEL_FIGS_DIR` — used on Sherlock to send outputs to `$SCRATCH`.
