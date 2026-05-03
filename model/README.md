# `model/` — log-linear IRT accumulator implementation

Reproducible code for the Bayesian accumulator model described in
[`outputs/model_explainer.pdf`](../outputs/model_explainer.pdf).

## Layout

```
model/
├── stan/
│   └── log_irt.stan          # The Stan model (collapsed xi parameterization,
│                             #   logit link, hyperprior-parameterized so the
│                             #   2x2 diagnostic variants share one file).
├── R/
│   ├── config.R              # Paths + default priors + fit config.
│   └── helpers.R             # All reusable functions (data loading, sim,
│                             #   Stan-data construction, fitting, summaries,
│                             #   plots). All drivers source this.
├── scripts/
│   ├── 01_recovery.R         # Parameter-recovery sim on synthetic data.
│   ├── 02_prepare_data.R     # Build stratified Wordbank subsample + Stan data.
│   ├── 03_diagnostics.R      # 2x2 diagnostic sweep + LOO-CV comparison.
│   └── 04_analyze.R          # Post-hoc plots + posterior summaries.
├── fits/                     # Output: .rds fits and cached Stan data.
└── figs/                     # Output: .png figures.
```

## Entry points — use the Makefile

All commands run from the project root. Each target is self-contained and
idempotent (re-running skips existing fits unless you `make clean-fits`).

| Target | What it does | Typical runtime |
|---|---|---|
| `make recovery` | Param-recovery on 250 children × 150 items of synthetic data | ~4 min |
| `make data` | Build Wordbank subsample (default 500 × 200) | ~20 sec |
| `make diagnostics` | Fit the 2×2 sweep (baseline, fix_delta, fix_s, both_fixed) + LOO | ~2–3 hrs |
| `make analyze` | RQ1/2/3/4 plots + tables for baseline | ~30 sec |
| `make analyze-all` | Same for all four variants | ~1 min |
| `make explainer` | Recompile `outputs/model_explainer.pdf` | ~5 sec |
| `make all` | recovery → data → diagnostics → analyze-all | ~3 hrs |
| `make clean-fits` / `clean-figs` / `clean` | Remove outputs | — |

Overriding defaults:
```
make data N_CHILDREN=1000 N_ITEMS=400
```

## Workflow

1. **First time:** `make recovery` (verifies the model can recover known
   parameters on simulated data; this is the green-light before running on
   real data).
2. **Prepare the Wordbank subsample:** `make data`.
3. **Run the diagnostic sweep:** `make diagnostics`.
4. **Inspect:** `make analyze-all` produces plots + tables in `outputs/figs/`.

## Migrating from the pre-refactor flat layout

The first diagnostic run wrote outputs to `model/wordbank_fit_*.rds` (flat
layout). To move them into `fits/` under the new naming convention:

```
make migrate
```

Idempotent — safe to run any number of times.

## Extending

All variants share one Stan file; new variants only need a new entry in
`variant_hyperpriors()` in `R/helpers.R`. The default priors are declared
in `R/config.R:DEFAULT_PRIORS`.
