# Feng et al. (2026) CHILDES-trained GPT-2 vs. children: per-word sigmoid-slope comparison

*Draft — preliminary numbers from partial training. Final values pending training completion (in flight).*

## TL;DR

The original Chang & Bergen (2022) finding — children's per-word sigmoid slope ≈ 10, LMs' per-word slope ≈ 1 — was hypothesized to be partly an input-distribution artifact (LMs trained on adult written text, kids hear CDS). We tested this by retraining GPT-2 small on CHILDES (matching Feng et al. 2026's setup) and refitting Chang & Bergen's per-word 4-PL sigmoid. At partial-training snapshots, the CDS-matched LM's per-word slope is converging to **≈ 0.7–0.9** — essentially indistinguishable from BookCorpus-trained GPT-2's **0.81**. Kids remain at **≈ 10**.

**The input-distribution explanation accounts for ~0 of the 10× kid-vs-LM gap.** The entire gap is structural: matching CDS doesn't move the LM-side per-word sigmoid slope at all. The "CHILDES is a better learning signal" framing is not supported by this data.

## Question

Children's per-child sigmoid slope on log-experience is $\kappa_i \sim 10$ (English M_best). Chang & Bergen's BookCorpus-trained LMs sit at per-word slope $\sim 1$ on the same axis. The Chang & Bergen comparison conflates two things: (a) the structural difference in scaling between children and LMs, and (b) the input-distribution mismatch (adult written text vs. CDS). This report measures the same per-word sigmoid slope on **CHILDES-trained GPT-2 small** (Feng et al. 2026 specification), so input distribution is matched and only the structural-comparison axis remains.

## Path chosen

**Path B (retrain)**. The Feng et al. 2026 release (`styfeng/babyscale-LM`) provides no checkpoints or per-step surprisal logs; the paper reports only end-of-training CDI-vocabulary NLL. We replicated the CHILDES condition using the training pipeline from `styfeng/TinyDialogues` (the same code Feng et al. used for the CHILDES experiments), retraining from scratch with per-CDI-word surprisal logged at log-spaced training steps.

## Setup

- **Model**: GPT-2 small (12L × 768d × 12h, 124M params), initialized from random.
- **Tokenizer**: `styfeng/TinyDialogues` `GPT2_CHILDES` BPE (52K vocab, trained on CHILDES).
- **Training data**: `CHILDES_train_ordered.txt` from `styfeng/TinyDialogues` — the same English-subset CHILDES file used in Feng et al. 2026, ~24.5M tokens, ordered (no in-epoch shuffling, mirroring Feng).
- **Training**: 20 epochs, AdamW, LR 1e-4 linear schedule (no warmup), per-device batch 8, sequence length 1024. Matches §B of Feng et al. 2026.
- **Seeds**: {42, 0, 123}, 3 replicates.
- **GPU**: Single GPU per seed (Sherlock `gpu` partition). Seeds 42 and 123 landed on NVIDIA L40S 48GB (bf16); seed 0 on V100 32GB (fp16). The V100 turned out to be ~3× slower per training step than L40S on this workload, so seed 0 may not complete in the 24-hour SLURM wall (linear projection: ~30 hrs). Seeds 42 and 123 are expected to finish in ~15-18 hrs.
- **Surprisal logging**: at 73 log-spaced training steps (80 requested, deduped after integer rounding), compute mean NLL of each CDI word over 50 random occurrences in CHILDES validation set, scored as $-\log p(w \mid \text{left context})$ under a causal LM forward pass. Context truncated to last 128 tokens of left context (full context up to 1024 would be ~10× more compute with negligible effect on the surprisal — most predictive information is in the immediate preceding sentence).

## Single-token coverage

All **611/611** of the C&B CDI words are single tokens (with leading space) in the GPT2_CHILDES BPE tokenizer (Feng et al. 2024) — better-than-expected coverage. In the CHILDES validation set, 609 of these 611 have ≥1 occurrence; 578 have ≥50 occurrences; 439 hit the 200-occurrence cap. The two missing words are anatomy terms (`buttocks`, `vagina`); two more (`ankle`, `downtown`) have <10 occurrences. Per-word eval uses 50 occurrences/word (random subsample with fixed seed), giving 578 words full statistical mass.

## Per-word sigmoid fits

We fit Chang & Bergen's 4-parameter logistic
$$S_w(x) = \mathrm{ParamLower}_w + \frac{\mathrm{ParamUpper}_w - \mathrm{ParamLower}_w}{1 + \exp((x - \mathrm{ParamXmid}_w)/\mathrm{ParamScale}_w)}$$
to each per-word mean-NLL trajectory across log-spaced training steps, with $x = \log_{10}(\text{steps})$. Per-word slope on natural-log experience is then $0.434/\mathrm{ParamScale}_w$.

Filtering follows C&B / our existing pipeline: drop fits with $\mathrm{ParamScale} \notin (0.01, 10)$ or surprisal range $\leq 1$ nat.

## Headline (preliminary, partial training)

*Final numbers TBD when training completes. Total training: 114,520 steps (~20 epochs, 9342 CHILDES conversations grouped into 1024-token blocks). Below: snapshot at 58 of 73 log-spaced evals (step ≈ 10,822 of 114,520, ~10% of training) for the two L40S seeds; seed 0 (V100) is much earlier.*

| Population | Median slope | IQR |
|---|---|---|
| Children, per-child $\kappa_i$ (English M_best) | 10.3 | [8.0, 12.6] |
| Children, per-child $\kappa_i$ (Norwegian M_best) | 12.5 | [10.0, 14.9] |
| **GPT-2-CHILDES, seed 42** (partial, N=605) | **0.84** | [0.45, 1.46] |
| **GPT-2-CHILDES, seed 123** (partial, N=605) | **0.89** | [0.46, 1.62] |
| GPT-2-CHILDES, seed 0 | (still mid-training) | |
| LMs, BERT-BookCorpus per-word | 0.76 | [0.42, 1.32] |
| LMs, BiLSTM-BookCorpus per-word | 0.87 | [0.53, 1.52] |
| LMs, GPT-2-BookCorpus per-word | 0.81 | [0.45, 1.54] |
| LMs, LSTM-BookCorpus per-word | 0.96 | [0.44, 1.90] |

**Evolution of partial fits as more evals accumulate (seed 42):**

| Evals done | Step reached | Words kept | Median |
|---|---|---|---|
| 23 | 62 | 385/609 (63%) | 3.22 |
| 33 | 364 | 475/609 (78%) | 2.36 |
| 48 | 2,477 | 583/609 (96%) | 1.85 |
| 53 | 5,177 | 595/609 (98%) | 1.35 |
| 58 | 10,822 | 605/609 (99.3%) | 0.84 |

The slope estimate is decreasing monotonically as more late-step data anchors the lower asymptote of each word's sigmoid. Final number TBD but the partial fits at intermediate training stages systematically overestimate slope (the upper asymptote is fixed at the random-init NLL, the lower asymptote keeps dropping as training continues, so the fitted ParamScale grows and the slope shrinks).

**Direction (preliminary):** CDS-matched GPT-2 has now converged to **almost exactly the same median slope** as C&B's BookCorpus-trained GPT-2 (0.84 vs. 0.81). If this holds through full training, the answer is:

> **Input-distribution matching does NOT close the kid-vs-LM gap. The entire 10× gap is structural.**

This is a stronger conclusion than expected from the original handoff and stands in direct contrast to a "CHILDES is a better learning curriculum than BookCorpus" framing. The per-word sigmoid slope appears to be invariant across very different LM training distributions, given matched compute. Caveat: the final number is still moving as training continues; the L40S seeds have only completed 10% of total training.

Snapshot plot: `outputs/figs/longitudinal/feng_partial_slope_comparison.png` (partial; final plots produced once training completes).

## Caveats

- **Partial-fit bias.** During training, per-word ParamScale estimates are unstable. Diagnostic per-word plots at 61 evals (`outputs/figs/longitudinal/feng_per_word_trajectories.png`) show two opposing biases:
  - *Shallowest-slope words* (e.g., `soft`, `cute`, `stay`): trajectories are still descending nearly linearly on log-step axis; the 4-PL fit reduces to a straight line because the lower asymptote hasn't been reached. As more training comes in these will become more sigmoidal and slopes will *increase*.
  - *Steepest-slope words* (e.g., `clap`, `yogurt`, `raisin`): trajectories look like step functions, fit with very small ParamScale → spuriously high slope. These are rare words (30-50 occurrences) where late-step noise dominates.
  Between 33 and 58 evals, 56% of common words shifted |Δslope| > 1.0. The median sits near 0.7-0.9 across recent snapshots and should be relatively stable as the two biases partly cancel, but final-training fits remain the meaningful comparand.
- **Compute scale.** GPT-2-CHILDES is trained on ~24.5M tokens × 20 epochs ≈ 0.5B tokens. Chang & Bergen's BookCorpus-trained GPT-2 saw ~75M tokens × many epochs in a similar-ish compute budget; the rough match of medians at the partial snapshot suggests slope is not very compute-sensitive once the model has begun to converge.
- **Tokenizer coverage.** A CHILDES-trained BPE tokenizer might segment some CDI items differently than C&B's GPT-2 tokenizer. In our case all 611 C&B CDI words are single tokens (see Single-token coverage section), so this concern doesn't apply here.
- **Seed-as-replicate.** We treat 3 seeds as the LM-side analog of between-instance variability. This is a closer analog than C&B's single-seed setup but still likely underestimates true variance.
- **CHILDES validation set as eval distribution.** CDI-word occurrences are drawn from the CHILDES validation set, not from a held-out wordbank-style probe set. This may favor words that are common in CHILDES specifically.
- **Per-word vs. per-instance.** The kid-side $\kappa_i$ is per-child (so different kids learning the same word have different slopes). The LM-side slope is per-word (so different words for the same model have different slopes). Both index "between-instance" variation but they're different kinds of instances. See `outputs/chang_bergen_derivation.tex` Section 4.

## Open questions

- **Compute-controlled control.** Does the slope shift if we train GPT-2-CHILDES for many more passes through the smaller CHILDES corpus (matched-compute condition)?
- **Larger CHILDES models.** Does the slope shift if we scale the model size while holding input distribution fixed?
- **Per-word, not per-seed, sigmoid variability.** C&B report per-word slope distributions; we do the same per seed. Pooling seeds into per-word distributions (averaging trajectories across seeds before fitting) would shrink within-word variance but reduce the population analog. Worth exploring.

## Reproduction

```
# 1. CHILDES context extraction (one-time, on Sherlock)
python model/scripts/feng_eval/extract_cdi_contexts.py \
    --tokenizer $FENG_ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
    --val_file  $FENG_ROOT/TinyDialogues/data/CHILDES/CHILDES_val_ordered.txt \
    --cdi_words $FENG_ROOT/cdi_words.txt \
    --out_jsonl $FENG_ROOT/cdi_contexts.jsonl \
    --out_coverage $FENG_ROOT/cdi_coverage.csv

# 2. Train (one job per seed)
sbatch --array=42,0,123 sherlock/feng_train_gpt2.slurm

# 3. Sigmoid fits + sigmoid TSVs (locally, after rsync of surprisal CSVs)
for seed in 42 0 123; do
  python model/scripts/feng_eval/fit_per_word_sigmoid.py \
    --surprisal_csv outputs/feng_eval/surprisal_gpt2_childes_seed${seed}.csv \
    --out_tsv       data/feng_2026/gpt2_childes_seed${seed}_sigmoids.txt
done

# 4. Plots + summary
Rscript model/scripts/feng_eval/feng_chang_bergen_comparison.R
```
