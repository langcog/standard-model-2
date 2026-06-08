# Data-variance pilot on Marlowe — runbook

**Drafted 2026-06-08.** Goal: answer the gating design question for the LLM
variability experiment — *does the identity of the CHILDES training data move
the per-word Chang & Bergen sigmoid slope, holding everything else fixed?* —
by training GPT-2 small on **two disjoint 10M-word CHILDES chunks** and
comparing their per-word slope distributions.

This decides whether the main study's σ_data axis can use CHILDES subsamples
or must move to TinyDialogues / smaller subsamples.

## The design (what's held fixed, what varies)

| Variable | Value | Why |
|---|---|---|
| Architecture | GPT-2 small (124M, Feng §B config) | match the 3-seed Sherlock baseline |
| Seed | 42 (both chunks) | data is the *only* independent variable |
| Tokenizer | `GPT2_CHILDES` (trained on full corpus) | hold vocab fixed; isolate data identity |
| Eval set | one `cdi_contexts.jsonl` from the held-out val set | identical probe for both chunks; no leakage asymmetry |
| Epochs | 20 | match Feng spec → comparable plateau |
| **Training data** | **chunk A vs chunk B (disjoint)** | **the manipulation** |

Two clean outcomes:
- **Slopes ≈ identical across chunks** → σ_data is small even between
  systematically different halves of CHILDES → CHILDES subsamples are a
  weak source of variance; main study should use TinyDialogues (more disjoint
  draws) or accept small σ_data.
- **Slopes differ noticeably** → data identity matters at 10M scale →
  CHILDES is a viable σ_data substrate; scale up the n on CHILDES.

(Reference: the 3 same-data different-seed Sherlock runs gave medians 0.72 /
0.74 / 0.74 — so any cross-chunk spread materially above ~0.02 is a data effect,
not seed noise.)

## Data provenance — fully self-contained, no Sherlock/Steven dependency

Everything is in the **public** `styfeng/TinyDialogues` repo (Feng et al. 2024):
- `data/CHILDES_data.zip` → unzips to `CHILDES_{train,val}_ordered.txt`
  (~24.5M words ≈ ~47M GPT2_CHILDES BPE tokens)
- `tokenizers/GPT2_CHILDES/` — the trained tokenizer
- `tokenizers/GPT2-small_config/config.json` — the model config

> **Gotcha:** `CHILDES_data.zip` is a **Git-LFS** object. `git clone` alone
> fetches a pointer; you must `git lfs install` then `git lfs pull`. The staging
> script handles this.

The feng_eval pipeline (`train_gpt2_childes.py`, `surprisal_callback.py`,
`extract_cdi_contexts.py`, `fit_per_word_sigmoid.py`) and the 611-word
`cdi_words.txt` are on **master** of `langcog/standard-model-2`.

## New artifacts created this session

- `model/scripts/feng_eval/cdi_words.txt` — the 611 C&B CDI words (from the
  `Token` column of `data/chang_bergen_2022/bert_sigmoids.txt`)
- `model/scripts/feng_eval/make_disjoint_chunks.py` — splits the corpus into two
  disjoint chunks by BPE-token budget (modes: `contiguous` | `random`)
- `model/scripts/feng_eval/marlowe/stage_marlowe.sh` — one-shot staging on a
  login node (clone → LFS pull → unzip → venv → contexts → split)
- `model/scripts/feng_eval/marlowe/train_gpt2_childes.slurm` — preempt-partition
  array job, one task per chunk

## Run sequence on Marlowe

```bash
# on a LOGIN node
export PILOT_ROOT=/scratch/$USER/llm_var_pilot      # <-- your scratch path
bash standard-model-2/model/scripts/feng_eval/marlowe/stage_marlowe.sh
# inspect the split before committing GPU hours:
cat $PILOT_ROOT/chunks/split_manifest.json | python3 -m json.tool | head -20

# launch BOTH chunks (one GPU each)
cd $PILOT_ROOT
sbatch --array=0,1 \
  standard-model-2/model/scripts/feng_eval/marlowe/train_gpt2_childes.slurm
```

After both finish (~2h each on H100), fit sigmoids and compare:

```bash
source $PILOT_ROOT/venv/bin/activate
FENG=$PILOT_ROOT/standard-model-2/model/scripts/feng_eval
for C in A B; do
  python $FENG/fit_per_word_sigmoid.py \
    --surprisal_csv $PILOT_ROOT/runs/surprisal_gpt2_childes_chunk${C}_seed42.csv \
    --out_tsv $PILOT_ROOT/runs/chunk${C}_sigmoids.txt
done
# headline: compare median + IQR of 0.434/ParamScale between chunkA and chunkB
```

## Compute budget

- Per run: 10M words ≈ 19M BPE ≈ ~46k steps (batch 8 × 20 epochs) ≈ **~2h on H100**
- 2 runs, 1 GPU each, can run concurrently → **~2h wall, ~4 GPU-hours total**
- Trivial against the pilot allocation; fits the 12h preempt cap with huge margin

## Preemption handling (preempt partition)

`preempt` has a 12h wall and can evict a job within 15 min if higher-priority
work arrives. Default job mode is `--no_save` + **resubmit on eviction** (a ~2h
job usually survives). The surprisal CSV is append/resumable, so a restart-from-0
self-heals.

If preemption bites repeatedly, switch to checkpoint+resume (`ROBUST=1` in the
SLURM script) — this needs a one-time driver patch:

```python
# train_gpt2_childes.py: add the arg
ap.add_argument("--resume_if_present", action="store_true")
ap.add_argument("--save_steps", type=int, default=5000)
ap.add_argument("--save_strategy", default="epoch")
# ...thread save_strategy/save_steps into TrainingArguments, then:
import os, glob
ckpt = bool(glob.glob(os.path.join(args.output_dir, "checkpoint-*"))) and args.resume_if_present
trainer.train(resume_from_checkpoint=ckpt)
```

## Open questions for Mike (block the actual launch)

1. **Marlowe account string.** The job needs `-A marlowe-<ProjectID>`. With
   preempt/basic access this is usually `marlowe-<id>` (no `-pm01/-pl01`
   suffix). What's the ID? (Currently a `PROJECTID` placeholder in the SLURM.)
2. **Scratch path.** I defaulted `PILOT_ROOT=/scratch/$USER/llm_var_pilot`.
   Confirm Marlowe's scratch convention (the docs name DDN Lustre/Intelliflash
   but not the mount path).
3. **Python provisioning.** Staging assumes a venv on a `python/3.12` module.
   If Marlowe steers you to conda or Apptainer instead, I'll swap the env block.
4. **Chunk design.** Default is `contiguous` disjoint (first-10M-words vs
   last-10M-words, order preserved — "two developmental slices"). Alternative
   is `random` disjoint (pure subsample variance, removes the early/late
   content confound). Contiguous is the stronger "does data matter at all"
   test; random is the cleaner analog of "different children's draws." Easy to
   run both later — which for the first shot?

## Why this is the right first move

- One cheap, decisive experiment (~4 GPU-hours) that gates the entire main-study
  data axis.
- Exercises the full Sherlock→Marlowe port (env, modules, LFS, SLURM) on a
  low-stakes 2-job run before any 90-run campaign — exactly the "master Marlowe"
  value in the application.
- Reuses the proven pipeline unchanged; the only new code is the deterministic
  chunk splitter.
