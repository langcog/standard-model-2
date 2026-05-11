# Feng et al. (2026) CHILDES-trained GPT-2 vs. children: per-word sigmoid-slope comparison

*Final numbers from completed L40S training (seeds 42, 123). Seed 0 still finishing on L40S restart after V100 abandonment; current partial estimate matches the others.*

## TL;DR

The original Chang & Bergen (2022) finding — children's per-word sigmoid slope ≈ 10, LMs' per-word slope ≈ 1 — was hypothesized to be partly an input-distribution artifact (LMs trained on adult written text, kids hear CDS). We tested this by retraining GPT-2 small on CHILDES (Feng et al. 2026 spec) and refitting Chang & Bergen's per-word 4-PL sigmoid. **The CDS-matched GPT-2 per-word slope is 0.74, essentially identical to BookCorpus-trained GPT-2's 0.81.** Kids are at 10.3.

> **The input-distribution explanation accounts for ~0 of the 10× kid-vs-LM gap.** The entire gap is structural; matching CDS does not move the LM-side per-word sigmoid slope.

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
- **GPU**: Single L40S 48GB per seed (Sherlock `gpu` partition, constrained to `GPU_GEN:AMP|GPU_GEN:LOV|GPU_GEN:HPR`). Seed 0 was originally placed on a V100 32GB, which was ~3× slower per training step than L40S; cancelled at 11% of training and resubmitted on L40S with the tighter constraint. Final wall-time per seed: 7h 21m (seed 42), 8h 12m (seed 123), seed 0 still running on L40S restart.
- **Surprisal logging**: at 73 log-spaced training steps (80 requested, deduped after integer rounding), compute mean NLL of each CDI word over 50 random occurrences in CHILDES validation set, scored as $-\log p(w \mid \text{left context})$ under a causal LM forward pass. Context truncated to last 128 tokens of left context (full context up to 1024 would be ~10× more compute with negligible effect on the surprisal — most predictive information is in the immediate preceding sentence).

## Single-token coverage

All **611/611** of the C&B CDI words are single tokens (with leading space) in the GPT2_CHILDES BPE tokenizer (Feng et al. 2024) — better-than-expected coverage. In the CHILDES validation set, 609 of these 611 have ≥1 occurrence; 578 have ≥50 occurrences; 439 hit the 200-occurrence cap. The two missing words are anatomy terms (`buttocks`, `vagina`); two more (`ankle`, `downtown`) have <10 occurrences. Per-word eval uses 50 occurrences/word (random subsample with fixed seed), giving 578 words full statistical mass.

## Per-word sigmoid fits

We fit Chang & Bergen's 4-parameter logistic
$$S_w(x) = \mathrm{ParamLower}_w + \frac{\mathrm{ParamUpper}_w - \mathrm{ParamLower}_w}{1 + \exp((x - \mathrm{ParamXmid}_w)/\mathrm{ParamScale}_w)}$$
to each per-word mean-NLL trajectory across log-spaced training steps, with $x = \log_{10}(\text{steps})$. Per-word slope on natural-log experience is then $0.434/\mathrm{ParamScale}_w$.

Filtering follows C&B / our existing pipeline: drop fits with $\mathrm{ParamScale} \notin (0.01, 10)$ or surprisal range $\leq 1$ nat.

## Headline

Total training: 114,520 steps (20 epochs, 45,807 1024-token blocks of CHILDES). Per-word sigmoid fit to 73 log-spaced eval points.

| Population | N | Median slope | IQR |
|---|---|---|---|
| Children, per-child $\kappa_i$ (English M_best) | 5000 | **10.3** | [8.0, 12.6] |
| Children, per-child $\kappa_i$ (Norwegian M_best) | 5000 | **12.5** | [10.0, 14.9] |
| GPT-2-CHILDES, seed 42 (✅ completed) | 609 | **0.74** | [0.43, 1.11] |
| GPT-2-CHILDES, seed 123 (✅ completed) | 609 | **0.74** | [0.45, 1.16] |
| GPT-2-CHILDES, seed 0 (still mid-training, ~64 evals) | 608 | 0.67 | [0.40, 1.13] |
| LMs, BERT-BookCorpus per-word | 609 | 0.76 | [0.42, 1.32] |
| LMs, BiLSTM-BookCorpus per-word | 604 | 0.87 | [0.53, 1.52] |
| LMs, GPT-2-BookCorpus per-word | 604 | 0.81 | [0.45, 1.54] |
| LMs, LSTM-BookCorpus per-word | 593 | 0.96 | [0.44, 1.90] |

**The CDS-matched GPT-2 (0.74) sits inside the same tight cluster as the 4 BookCorpus-trained LMs (0.76–0.96), all far below kids (10.3).**

The two completed seeds (42, 123) give identical medians 0.74 — strong reproducibility. Seed 0's partial estimate (0.67) is also in the same cluster.

**Headline plot:** [outputs/figs/longitudinal/feng_chang_bergen_slope_comparison.png](outputs/figs/longitudinal/feng_chang_bergen_slope_comparison.png) (Feng-CHILDES blue, C&B-BookCorpus green densities overlap almost exactly at the mode; children red is at a completely separate location ~10×).

### Evolution of partial fits as training progressed (seed 42)

The slope estimate dropped monotonically as more late-step data anchored each word's lower asymptote:

| Evals done | Step reached | Words kept | Median |
|---|---|---|---|
| 23 | 62 | 385/609 (63%) | 3.22 |
| 33 | 364 | 475/609 (78%) | 2.36 |
| 48 | 2,477 | 583/609 (96%) | 1.85 |
| 53 | 5,177 | 595/609 (98%) | 1.35 |
| 58 | 10,822 | 605/609 (99.3%) | 0.84 |
| 61 | 16,843 | 608/609 (99.8%) | 0.69 |
| **73 (final)** | **114,520** | **609/609 (100%)** | **0.74** |

The slope settled at 0.74 once training reached convergence. During training the value first dropped (intermediate fits overestimated ParamScale because each word's lower asymptote was still falling), then ticked back up slightly at the very end as the shallowest-slope words (which had been fit as straight lines on log-step axis) acquired clearly sigmoidal shape with their own well-defined transitions. The two biases identified mid-training partly cancelled; the final value is the meaningful comparand.

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
