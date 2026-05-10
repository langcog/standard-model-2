# Feng et al. (2026) CHILDES-trained GPT-2 vs. children: per-word sigmoid-slope comparison

*Draft. To be filled in after training completes.*

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
- **GPU**: Single GPU (Sherlock `gpu` partition, H100/A100), bf16.
- **Surprisal logging**: at ~80 log-spaced training steps, compute mean NLL of each CDI word over up to 200 occurrences in CHILDES validation set, scored as $-\log p(w \mid \text{left context})$ under a causal LM forward pass.

## Single-token coverage

All **611/611** of the C&B CDI words are single tokens (with leading space) in the GPT2_CHILDES BPE tokenizer (Feng et al. 2024) — better-than-expected coverage. In the CHILDES validation set, 609 of these 611 have ≥1 occurrence; 578 have ≥50 occurrences; 439 hit the 200-occurrence cap. The two missing words are anatomy terms (`buttocks`, `vagina`); two more (`ankle`, `downtown`) have <10 occurrences. Per-word eval uses 50 occurrences/word (random subsample with fixed seed), giving 578 words full statistical mass.

## Per-word sigmoid fits

We fit Chang & Bergen's 4-parameter logistic
$$S_w(x) = \mathrm{ParamLower}_w + \frac{\mathrm{ParamUpper}_w - \mathrm{ParamLower}_w}{1 + \exp((x - \mathrm{ParamXmid}_w)/\mathrm{ParamScale}_w)}$$
to each per-word mean-NLL trajectory across log-spaced training steps, with $x = \log_{10}(\text{steps})$. Per-word slope on natural-log experience is then $0.434/\mathrm{ParamScale}_w$.

Filtering follows C&B / our existing pipeline: drop fits with $\mathrm{ParamScale} \notin (0.01, 10)$ or surprisal range $\leq 1$ nat.

## Headline (preliminary, partial training)

*Final numbers TBD when training completes. Below: snapshot at ~33 of 80 log-spaced evals (step ≈ 364 of 60000) for the two L40S seeds; seed 0 (V100) is much earlier.*

| Population | Median slope | IQR |
|---|---|---|
| Children, per-child $\kappa_i$ (English M_best) | 10.3 | [8.0, 12.6] |
| Children, per-child $\kappa_i$ (Norwegian M_best) | 12.5 | [10.0, 14.9] |
| **GPT-2-CHILDES, seed 42** (partial, N=475) | **2.36** | [1.43, 3.75] |
| **GPT-2-CHILDES, seed 123** (partial, N=499) | **2.34** | [1.41, 3.75] |
| GPT-2-CHILDES, seed 0 | (still on early evals) | |
| LMs, BERT-BookCorpus per-word | 0.76 | [0.42, 1.32] |
| LMs, BiLSTM-BookCorpus per-word | 0.87 | [0.53, 1.52] |
| LMs, GPT-2-BookCorpus per-word | 0.81 | [0.45, 1.54] |
| LMs, LSTM-BookCorpus per-word | 0.96 | [0.44, 1.90] |

**Direction (preliminary):** CDS-matched LM training shifts the per-word sigmoid-slope distribution upward by roughly **2–4×** versus BookCorpus-trained LMs. The gap to kids ($\kappa \approx 10$) narrows but does **not** close — kids remain ~4–5× steeper than CDS-matched GPT-2. The structural-difference hypothesis is winning.

Snapshot plot: `outputs/figs/longitudinal/feng_partial_slope_comparison.png` (partial; final plots produced once training completes).

## Caveats

- **Compute scale.** GPT-2-CHILDES is trained on ~24.5M tokens × 20 epochs ≈ 0.5B tokens. Chang & Bergen's BookCorpus-trained GPT-2 saw ~75M tokens × many epochs in a larger compute budget. We are varying input distribution *and* compute simultaneously here; the cleanest control would be a CHILDES-trained model run to a matched compute budget (deferred — see open questions).
- **Tokenizer coverage.** A CHILDES-trained BPE tokenizer may segment some CDI items differently than C&B's GPT-2 tokenizer; we report single-token coverage in the table above.
- **Seed-as-replicate.** We treat 3 seeds as the LM-side analog of between-instance variability. This is a closer analog than C&B's single-seed setup but still likely underestimates true variance.
- **CHILDES validation set as eval distribution.** CDI-word occurrences are drawn from the CHILDES validation set, not from a held-out wordbank-style probe set.

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
