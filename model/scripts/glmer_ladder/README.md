# glmer_ladder

Self-contained model-comparison pipeline. Tests the four design choices in
M_best (acceleration, log vs linear age, kid intercept, kid slope) on as many
longitudinal CDI datasets as we have. Pure-frequentist `glmer` (nAGQ=0); no
Bayesian machinery and no σ_r decomposition (those belong to the input-share
section of the talk, not this section).

## The model ladder

Cumulative nested sequence, 7 models per language:

| Tag | Formula | What it tests |
|---|---|---|
| `A`     | `produces ~ offset(log_age) + (1 \| item)` | pure unit accumulator (κ=1, no kid RE) |
| `B_log` | `produces ~ 1 + log_age + (1 \| item)` | free κ on log_age |
| `B_lin` | `produces ~ 1 + age_c   + (1 \| item)` | free slope on linear age |
| `C_log` | `+ (1 \| child)` | + random kid intercept |
| `C_lin` | `+ (1 \| child)` | + random kid intercept |
| `D_log` | `(1 + log_age \| child)` | + random kid slope (= M_best logit) |
| `D_lin` | `(1 + age_c   \| child)` | + random kid slope |

Key questions each Δ-AIC answers:

- `A → B_log`: does free κ beat the unit accumulator?
- `B_lin → B_log`: log vs linear age (the power-law-vs-exponential test)
- `B → C`: does per-kid intercept matter?
- `C → D`: does per-kid slope (σ_ζ > 0) matter?

## Candidate languages (≥100 kids with ≥2 admins)

From `outputs/glmer_ladder/00_language_survey.csv`:

| Language | Long kids | Forms |
|---|---|---|
| English (American) | 1874 | WG + WS |
| Norwegian | 1676 | WG + WS |
| Finnish | 236 | WGProd + WS |
| French (Quebecois) | 179 | WG + WS |
| Japanese | 178 | WG + WS |
| Spanish (Mexican) | 119 | WG + WS |
| French (French) | 111 | WG + WS |

Skipped: English (British) (TEDS/Oxford, non-standard instrument); Dutch
(too many form variants to combine cleanly).

Total compute: 7 langs × 7 models = **49 fits**.

## Pipeline

```
00_survey_languages.R   # one-shot: which languages have longitudinal data?
01_extract_one.R        # pull one language's tidy data → fits/glmer_ladder/data_<slug>.rds
01_extract_all.R        # loops 01 over all 7 languages
02_fit_one.R            # fit ONE model on ONE language → summary csv + fit rds
03_aggregate.R          # combine summaries → ΔAIC table + figure
04a_simulate.R          # BLUP-bootstrap predictions (slow) → sim_cache.rds
04b_plot.R              # read cache, build mega-figures (fast; iterate here)
```

The 04 step is split so you can iterate on plot aesthetics without
re-running the 500-kid bootstrap across all 42 fits. Run 04a once after
fits land; then re-run 04b as many times as you like.

## Compute on Sherlock

Each fit runs on 1 core (glmer is single-threaded). 49 SLURM array tasks
each requesting 1 core × 64 GB × 12 hr.

```
# On laptop, locally:
Rscript model/scripts/glmer_ladder/01_extract_all.R       # ~5 min
rsync -avz fits/glmer_ladder/data_*.rds sherlock:standard_model_2/fits/glmer_ladder/

# On Sherlock:
bash sherlock/glmer_ladder_submit.sh    # writes manifest, prints sbatch line
sbatch --array=1-49 sherlock/glmer_ladder.slurm

# After fits complete, locally:
bash sherlock/glmer_ladder_sync.sh                        # pulls summaries
Rscript model/scripts/glmer_ladder/03_aggregate.R         # builds figure
```

## Two pre-flight tests (Mike's spec)

1. **Infrastructure smoke test** — one tiny cell (Japanese B_log, ~30 sec)
   submitted as `sbatch --array=1-1 sherlock/glmer_ladder.slurm` with the
   manifest restricted to that one row. Validates SLURM, paths, sync.
2. **Norwegian D_log timing test** — the biggest single cell. Submitted
   as a standalone sbatch (or `--array=14-14` once the full manifest is
   built and we know which index NO D_log is). Validates 12 hr is enough
   and 64 GB RAM holds.

Only after both pass do we submit the full 49-cell array.
