---
name: gpt2-childes-training
description: Train GPT-2 (or similar small transformer) from scratch on CHILDES-scale corpora on academic SLURM clusters, logging per-CDI-word surprisal trajectories at log-spaced training steps. Use when extending Chang & Bergen (2022) / Feng et al. (2024, 2026) word-acquisition pipelines to new training data or new model variants. Covers GPU selection, Python env setup pitfalls, HF Trainer callback patterns, per-word surprisal eval design, and 4-PL sigmoid fitting. Generic across SLURM-based HPC clusters; Sherlock-specific bits are flagged.
---

# Training small transformer LMs on CHILDES-scale data with per-word surprisal trajectories

## When this skill applies

You're extending the Chang & Bergen (2022) word-acquisition pipeline — fitting per-CDI-word 4-PL sigmoids to mean-surprisal trajectories during LM training — to a new input distribution, new model, or new vocabulary. You want to:

- Train a small transformer (GPT-2 small / mini, ~125M params) from scratch on a CHILDES-scale corpus (~10–50M tokens), ~20 epochs, on a SLURM cluster.
- Log per-word surprisal at ~80 log-spaced training steps during the run.
- Fit 4-PL sigmoids to each trajectory and read off per-word slopes on `log10(steps)` (then convert to natural-log experience via `slope_natural = 0.434 / ParamScale`).

This skill captures the gotchas the SM2 feng_eval experiment ran into. Pipeline scripts and a one-shot runner live in `model/scripts/feng_eval/`; treat them as the reference implementation.

## The lessons

### 1. GPU selection: compute capability matters more than memory

GPT-2 small fits in 11 GB if you're careful, so cards from RTX 2080Ti up will *fit*. But cu121-built torch + bf16 / fp16 throughput varies wildly:

- **Volta (V100, sm_70)** is ~3× slower per training step than Ada (L40S). On a 114 k-step 20-epoch run that's the difference between fitting and not fitting in a 24 h wall.
- **Pascal (P100, TITAN Xp, sm_60/sm_61)** isn't reliably supported by recent torch wheels — cu121 ships PTX for sm_50 and sm_70+ by default. May OOM at trivial batch sizes due to lack of native fp16.
- **Ampere (A100, RTX 3090, sm_80/8.6)**, **Ada (L40S, sm_89)**, **Hopper (H100, sm_90)** all work well with bf16 and are ~3× faster than Volta per step.

**Use SLURM `--constraint` by GPU generation, not memory.** On Sherlock:

```bash
#SBATCH --constraint=GPU_GEN:AMP|GPU_GEN:LOV|GPU_GEN:HPR
```

Other clusters expose the constraint differently — look for `nvidia_smi --query-gpu=compute_cap` features in `scontrol show node`. Memory-only filters (`GPU_MEM:24GB|...`) admit V100 32GB, which will likely miss your wall.

A V100 is good enough for a smoke test but not for a 20-epoch real run. Specify the constraint on both the smoke and real SLURM scripts so jobs land on hardware that's *actually fast*, not just *big enough*.

### 2. Python env: a maze of pip pitfalls on glibc 2.17 hosts

CentOS 7 / RHEL 7 SLURM clusters often have glibc 2.17. This is too old for the latest miniconda — install will silently fail near the end. Options that *do* work:

- **Miniconda 24.x** (specifically `Miniconda3-py311_24.5.0-0`) — released before the glibc 2.28 bump.
- **A venv on top of a cluster Python module** (e.g., `python/3.12.1`). Cheaper than conda but requires the module to be loaded whenever the venv is activated (the venv's symlinked Python disappears if the module isn't loaded).

The Sherlock recipe that worked:

```bash
ml load python/3.12.1
python3 -m venv /scratch/$USER/feng_eval/venv
ml load cuda/12.6.1 cudnn/9.4.0
source /scratch/$USER/feng_eval/venv/bin/activate
export PIP_CACHE_DIR=/scratch/$USER/feng_eval/.pipcache
export TMPDIR=/scratch/$USER/feng_eval/.tmp
mkdir -p $PIP_CACHE_DIR $TMPDIR
# Downgrade pip to a stable version
pip install "pip<25"
# Install torch matched to cuda runtime
pip install "torch==2.4.*" --index-url https://download.pytorch.org/whl/cu121
# Stack — ALWAYS use --only-binary=:all: on academic clusters
pip install --only-binary=:all: \
    "numpy>=2.0,<2.3" "scipy>=1.13" "pandas>=2.2,<3" "pyarrow>=14" \
    "transformers>=4.40,<4.50" "datasets>=2.20,<4" \
    "accelerate>=0.30,<2" "tokenizers>=0.15,<0.22"
```

Critical lessons:

- **`--only-binary=:all:`** forces pip to use wheels and fail loudly if any package would need a source build. Without it, scipy/pyarrow/pandas will silently try to build from source — 15+ minute meson/cmake builds that often fail on minimal cluster toolchains.
- **`pip < 25`** — pip 25+ resolver has behaviors that make wrong source-build decisions on some constraint combinations.
- **Pin to known-working version ranges.** `numpy<2` looks innocent but forces an old scipy that has no cp312 wheels. `numpy>=2.0,<2.3 + pandas>=2.2,<3 + pyarrow>=14 + datasets>=2.20` all have wheels available and resolve together.
- **`--quiet` hides errors.** Don't use it during initial install; you won't see "package being built from source" warnings.
- **`$HOME` is small on most clusters** (Sherlock: 15 GB). Install everything in `$SCRATCH` — venv, pip cache, pip TMPDIR, HF cache. A torch install + transformers stack is ~5 GB.

### 3. SLURM gotchas

- **Source the LMOD init file in your sbatch script.** `/etc/profile.d/modules.sh` doesn't exist on Sherlock compute nodes; the correct path is `/share/software/user/open/lmod/lmod/init/bash`. Find yours via `which module` → trace back from `$LMOD_PKG`.
- **`#SBATCH --output=...` paths are relative to the sbatch invocation directory.** Either `cd` into the right place before `sbatch`, or use absolute paths in `--output=%x_seed%a_%j.out`.
- **Array tasks share the parent job ID's queue listing**, but each task gets its own `%j` for log naming. `sacct -j <parent>` shows them all.
- **Don't set `--gpus=N` if `N=1`** — it works but some clusters charge differently than `--gres=gpu:1`. Pick the one your group's documentation recommends.

### 4. HuggingFace Trainer integration

For transformers 4.40+:

- `evaluation_strategy` was deprecated in favor of `eval_strategy`. Use the new name.
- The `tokenizer=` kwarg in `Trainer.__init__` was deprecated; the new name is `processing_class=`. Both still work but emit warnings.
- `Trainer._get_train_sampler` is the override point for "no in-epoch shuffling" (e.g., curriculum experiments) — return a `SequentialSampler` over `self.train_dataset`.
- **bf16 vs fp16**: enable `bf16=True` on Ampere+; fall back to `fp16=True` on Volta. Don't request both — bf16 has fewer numerical issues but isn't supported on sm_70.

The `--keep_linebreaks` argument is needed if your corpus uses `\n` as discourse boundary markers (CHILDES does).

### 5. Surprisal-eval callback design (the high-leverage piece)

The naive design — eval ALL 611 CDI words × ALL ~200 occurrences × ALL training steps — is too slow. Design choices that matter:

**Truncate context to the last ~128 tokens.** GPT-2's per-token NLL barely depends on context past the immediate sentence; 1024-token contexts cost 8× more compute for no measurable change in surprisal estimates. Pre-extract contexts as JSONL with `ctx[-128:]` and `pos = len(ctx) - 1`.

**Random-subsample to ~50 contexts per word.** Chang & Bergen used ~512 — overkill. 50 gives SE on mean NLL well under 0.1 nat, which is far below the between-word signal. Deterministic seed (`random.Random(2026)`) for reproducibility.

**Length-bucket batches and right-pad.** Sort the flat (word, ctx, pos) list by `len(ctx)` before batching. With 30k contexts at batch 128, ~234 forward passes per eval; on L40S that's <10 s.

**Right-pad (not left-pad)** so GPT-2's absolute position embeddings sit at 0..L-1 — the positions the model was trained on for short sequences.

**Vectorize the NLL gather.** Don't loop over batch with `.item()` calls — that's a 200 ms hit per batch from GPU→CPU syncs. Instead:

```python
row_idx = torch.arange(B, device=device)
picked_logits = logits[row_idx, tgt_pos, :]                # [B, V]
picked_log_probs = torch.log_softmax(picked_logits.float(), dim=-1)
picked_nll = -picked_log_probs.gather(1, tgt_id.unsqueeze(1)).squeeze(1)  # [B]
nlls_cpu = picked_nll.detach().cpu().numpy()               # one sync per batch
```

Going from per-element `.item()` to batched gather sped our eval up ~5-10× on V100.

**Trigger on log-spaced step targets, not every N steps.** Generate the schedule once at train-begin from `state.max_steps`. After integer rounding the early dense steps collapse — expect 73 unique targets from a nominal 80, not 80. Don't hard-code the eval count.

**Fire one eval at on_train_end** to capture the very last step in case it wasn't in your scheduled set.

**Set `model.eval()` before evaluating and `model.train()` after.** Don't return early on `pos < 1` items — guard them via a valid-mask, fill their NLL with 0 (or skip via a `gather` after-mask), and they fall out cleanly without corrupting the next batch.

### 6. CSV resumability

If your callback writes a CSV in append mode, you have to handle the "file exists but no header" case (e.g., you cancelled and resubmitted, and the old file got truncated or moved). Backfill the header on init if the first line doesn't match the expected schema:

```python
needs_header = True
if os.path.exists(self.out_csv) and os.path.getsize(self.out_csv) > 0:
    with open(self.out_csv) as f:
        first = f.readline().rstrip("\n")
    needs_header = (first != "step,epoch,word,n_occurrences,mean_nll,sum_nll")
if needs_header:
    with open(self.out_csv, "w", newline="") as f:
        csv.writer(f).writerow([...schema...])
```

Without this, `csv.DictReader` will treat the first data row as the header and silently fit nothing.

### 7. 4-PL sigmoid fit interpretation

Use scipy `curve_fit` with bounds:

```python
def four_pl(x, upper, lower, xmid, scale):
    return lower + (upper - lower) / (1.0 + np.exp((x - xmid) / scale))
```

with `x = log10(steps)`, `bounds=([lo-5, lo-5, x.min()-3, 1e-3], [hi+5, hi+5, x.max()+3, 50])`.

**Filter** to keep only fits with `0.01 < ParamScale < 10` and `(ParamUpper - ParamLower) > 1 nat`. Without filtering you'll get spurious slope outliers from words with no learning signal.

**Partial-fit bias is real.** Mid-training fits **systematically overestimate ParamScale (underestimate slope)** because the lower asymptote is still falling; the fitter doesn't yet see where the trajectory plateaus. Don't read slope distributions until training has plateaued. The C&B comparand is *fully-trained*; your fit must be too. In the SM2 feng_eval run, the seed 42 median slope evolved:

| evals done | step | median slope |
|---|---|---|
| 23 | 62 | 3.22 |
| 33 | 364 | 2.36 |
| 53 | 5,177 | 1.35 |
| 73 (final) | 114,520 | **0.74** |

A 4× drop from intermediate to final. If you publish partial-training results, expect the final number to be much lower.

### 8. Monitoring the run

Don't poll every minute via `tail`/`wc`. Set up a single milestone-based watcher:

- Fire on every 25 evals (or some natural breakpoint) — not on every eval.
- Fire on terminal state transitions (`squeue` output transitions away from RUNNING).
- Fire on error patterns in stderr: `grep -qE "Traceback|OutOfMemoryError|RuntimeError|FAILED|Killed"`.

For a 73-eval log-spaced training run that takes ~8 hours, milestone-firing every 25 evals = 3 events from start to finish per seed, plus the completion signal. That's the right rate.

**Don't trust `nvidia-smi` snapshots showing 0% GPU util** — eval has CPU-bound phases (file reading, tokenization) where the GPU briefly sits idle. Confirm liveness by checking process state (`ps -p $pid -o state` should show R or D, not Z).

### 9. Picking total_steps and batch size

The HF data collator's `group_texts` chunks the tokenized corpus into `block_size`-token blocks (we use 1024). For 24.5M CHILDES tokens, that's ~24k blocks. At per-device batch 8, ~3000 steps/epoch × 20 epochs = ~60k steps in principle.

**In practice you'll see more** — the HF data pipeline pads/concatenates differently than your back-of-envelope calc. Our 9342 CHILDES conversations turned into 45,807 1024-token blocks → 114,520 steps total at batch 8 × 20 epochs. Use the value the Trainer logs at train_begin to size your log-spaced eval schedule.

### 10. Multi-seed array jobs

Use `--array=42,0,123` to launch 3 seeds in parallel as separate jobs. Each gets its own GPU. Treat each seed's slope distribution as a sample from "per-word slope distribution given (architecture, training distribution)" — i.e., the LM-side analog of between-instance variability.

3 seeds is the minimum for reliable inference; 5 is better if you can afford the compute. In the SM2 feng_eval run the three seed medians were 0.72, 0.74, 0.74 — strong reproducibility, confirming we have enough seeds.

## Reference implementation

`model/scripts/feng_eval/` in this repo:

- `extract_cdi_contexts.py` — one-shot context pre-extraction from a held-out val set
- `surprisal_callback.py` — HF TrainerCallback with the log-spaced schedule, batched gather, length-bucketing, header backfill
- `train_gpt2_childes.py` — self-contained training driver matching Feng et al. 2026 CHILDES condition
- `fit_per_word_sigmoid.py` — 4-PL fits via scipy curve_fit, outputs in Chang & Bergen's TSV schema
- `feng_chang_bergen_comparison.R` — kids + CHILDES + BookCorpus density plot
- `finalize.sh` — end-to-end one-command runner (rsync → fit → plot → summary)

SLURM scripts in `sherlock/`:

- `feng_smoke.slurm` — 30-min smoke test on a small CHILDES subset
- `feng_train_gpt2.slurm` — 24-h full training, array=42,0,123

## Quick reproduction checklist

Adapt to a new cluster by:

1. Replace the LMOD init path in the SLURM scripts (`source <path>`).
2. Replace the GPU constraint with whatever your cluster exposes for Ampere+ generations.
3. Replace `/scratch/$USER` with your cluster's scratch path.
4. If your cluster's Python module isn't 3.12, adjust the venv creation step and pip-wheel constraints accordingly.
5. Replace the SLURM `--constraint` with the equivalent for your scheduler.
6. The Python code, HF callback, sigmoid fit, and plot logic are cluster-agnostic.

A clean smoke test takes ~30 min and validates the entire pipeline before you commit to multi-hour real runs. Always smoke test on a new cluster.
