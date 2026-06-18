# Autoresearch acceleration pilot

A **feasibility pilot** (2026-06-17) for an exploratory idea: repurpose Karpathy's
[`autoresearch`](https://github.com/karpathy/autoresearch) (a 5-min single-GPU
nanochat micro-loop on the ClimbMix web corpus) to measure **per-word acquisition
acceleration** (the Chang & Bergen sigmoid slope κ) instead of hill-climbing
`val_bpb`. Long-term goal: a sandbox for "can we *engineer* an LM that accelerates
like a child" — try architectural/training tricks, see if median κ moves off ~1.

**This pilot only answered the prerequisite question: can we estimate per-word κ
from these tiny runs at all? Answer: yes.**

## Result

| run | steps | tokens | val_bpb | words resolving (≥1 nat) | median κ | mode |
|---|---|---|---|---|---|---|
| 5-min budget | 152 | ~80 M | 1.226 | 349/365 (96%) | 1.65 | ~1.0 |
| single pass (6 shards) | 580 | ~304 M | 1.057 | 340/365 (93%) | 1.79 | ~1.0 |

- κ ≈ 1.7 is **stable across run lengths** (not a short-range fit artifact, as I'd
  guessed). The distribution's **mode sits right at the κ=1 unit-accumulator line**,
  with a right-skewed tail of "faster" words pulling the median up.
- So the baseline ClimbMix GPT-2, on the single-pass distinct-token axis, is a
  near-unit accumulator at the mode — modestly above it overall. That's the
  starting point an acceleration hill-climb would push on.

## Pipeline (all idempotent; run from the autoresearch repo dir)

1. `adapt_a40.py` — autoresearch is **H100-designed** (FlashAttention-3 + torch
   2.9.1/cu128). This swaps FA3 → PyTorch SDPA (globally disabling the materializing
   math backend, which OOMs the A40) so it runs on Ampere; pair with our existing
   conda torch 2.4.1 and `DEVICE_BATCH_SIZE=64`. Architecture/math unchanged;
   per-layer sliding `window_size` ignored (full causal).
2. `cdi_words.txt` — 609 CDI words (Chang & Bergen set, = our sigmoid Token column).
3. `extract_cdi_contexts_ar.py` — single-token CDI words (365/609 survive the 8192
   BPE vocab; all ≥80 occurrences in the val shard) → fixed-length contexts JSONL.
4. `apply_probe_patch.py` + `cdi_probe.py` — per-word mean-NLL logger hooked into
   the nanochat `while` loop at log-spaced steps. Inert unless `CDI_PROBE=1`. Calls
   the model's `_orig_mod` (eager) to avoid recompiling the compiled training graph.
5. `single_pass_patch.py` — switch the stop condition from the 5-min wall-clock
   budget to one full pass (`epoch>=2`, with a step-cap backstop) and rescale the LR
   schedule to step-fraction. (Optional; the better axis for a κ measurement —
   controls *experience* rather than wall-clock, hardware-independent.)
6. `fit_kappa_ar.R` — 4-PL sigmoid per word on log10(step); κ = 0.434/scale.

## Reproduce (ccn2 A40 box)

```sh
# env: /data2/mcfrank/ladder/condaenv (torch 2.4.1 + tiktoken/rustbpe/pyarrow)
cd /data2/mcfrank/autoresearch
python prepare.py --num-shards 6           # ClimbMix shards + 8192 BPE tokenizer
python adapt_a40.py                         # FA3 -> SDPA
python extract_cdi_contexts_ar.py --cdi-words cdi_words.txt --out cdi_contexts.jsonl
python apply_probe_patch.py
python single_pass_patch.py
CDI_PROBE=1 CDI_PROBE_STEPS=700 CDI_SP_STEPS=700 CDI_PROBE_NPOINTS=30 \
  CDI_PROBE_CSV=word_surprisal_sp.csv CUDA_VISIBLE_DEVICES=0 python train.py
# then: Rscript fit_kappa_ar.R word_surprisal_sp.csv singlepass
```

## Caveats (why κ here isn't directly comparable to the paper's numbers)

- **Single-pass = undertrained**: at 580 steps the lower asymptote isn't fully
  reached, which can bias the per-word slope upward vs a *converged* measurement
  (our CHILDES ladder trains each budget to convergence; C&B trains 1M steps).
- **ClimbMix is web text, not CDS**; word frequencies/difficulty differ from child
  input. (Coverage is fine — CDI words are common — but it's not a child model.)
- 8192 BPE keeps only 365/609 CDI words as single tokens (the rest are split).
- This is a *step-axis* (within-run) measurement, distinct from the converged
  development-axis κ in the main paper.

Artifacts: `trajectory_{5min,singlepass}.csv` (raw per-word/step surprisals),
`word_surprisal_*_kappa_*.csv` (per-word fits), `*_kappa_singlepass.png` (density).
