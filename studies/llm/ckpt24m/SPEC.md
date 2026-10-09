# 24M-rung retrain with log-spaced checkpoints (3 seeds) — spec

**Status:** spec, not yet run (2026-10-09).

## Why
The HF ladder (`mcxfrank/childes-gpt2-ladder`) holds only each run's best-val endpoint.
Training kept one checkpoint at a time (`--save_total_limit 1`), so no intermediate weights
exist anywhere. That covers ccn2 `/data2`, Oak `flatlr_grammar-archive` (final step only)
and Sherlock scratch; Marlowe was not checked because it needs Duo. This run recreates the
weight trajectory of the largest rung so analyses that need *weights*, not just per-word
surprisal, become possible.

## What (fixed: identical to L6 / L6-archive)
| | |
|---|---|
| seeds | **42, 0, 123** (the original L1 seeds; each has an L6 24M endpoint on HF) |
| data | `chunks/ladder_seed{S}/CHILDES_ladder_seed{S}_24M.txt` (byte-identical L6 chunks) |
| val / tokenizer / config | `CHILDES_val_ordered.txt`, `GPT2_CHILDES`, `GPT2-small_config` |
| optimizer | AdamW 1e-4, linear decay over **40 epochs = 224,880 steps**, no warmup, batch 8, seq 1024, no shuffle, bf16 |
| stopping | early stopping, patience 4 on per-epoch `eval_loss` (`--max_eval_blocks 500`), load best |
| script | `train_gpt2_childes.py` (deployed copy on ccn2 verified identical to the repo copy) |

The command line is `worker.sh`'s, verbatim, plus one new flag. **Do not change epochs/LR**:
the linear decay is defined over the 40-epoch horizon, so ES stops while the LR is still about
75% of its peak. Changing `num_train_epochs` would change the trajectory being recorded.

Reference from L6-archive logs: 5,622 steps/epoch; best val at epoch 6 (seeds 0, 42) / 7
(seed 123); ES stop at epoch 10 / 10 / 11 (step 56,220–61,842); 4.3–4.7 h per run on an A40.

## The one code change (additive, default-off)
`--traj_save_steps FILE` → a `TrainerCallback` that, at each listed global step (on
`on_step_end`), and at **every epoch end**, writes
`{output_dir}/trajectory/step-{NNNNNN}/` via `model.save_pretrained(..., safe_serialization=True)`
cast to **bf16** (same dtype convention as the HF ladder), plus the tokenizer and a
`state.json` (step, epoch, current LR, last eval_loss).
- It writes outside `checkpoint-*`, so the Trainer's `save_total_limit` rotation and ES
  best-model logic are untouched.
- Weights only, no optimizer state (fp32 + Adam ≈ 1.5 GB/ckpt; not needed for analysis).
- Step 0 = initialization, saved in `on_train_begin`.

**Step grid** (`make_steps.py`): step 0, plus 64 log-spaced integer steps from 1 to
**67,464 (12 epochs)**, de-duplicated, plus the epoch ends 5,622·k. The horizon is 12 epochs,
not 40, because every seed stops by epoch 11; spacing out to 224,880 would waste about 40% of
the points on steps never reached. Expect ~75 checkpoints/seed × ~250 MB ≈ **19 GB/seed,
~56 GB total** (ccn2 `/data2` has 2.1 TB free).

Surprisal is computed **offline from the saved checkpoints** (reusing the
`surprisal_callback` scoring code), not in-loop. That keeps training identical to L6 and
puts weights and readout on exactly the same steps.

## Verification gates
1. **Smoke** (0.5M rung, seed 42, ~25 min): checkpoints appear at the expected steps, reload
   cleanly, and the run's final best checkpoint matches an untouched run of the same config.
   Weights won't be bit-identical (GPU nondeterminism), so compare per-word surprisal: r ≥ 0.98.
2. **Fidelity, per seed:** best-val checkpoint surprisal vs `fits/llm/ladder_bestval_finer.csv`
   (24M row for that seed): r ≥ 0.98, mean within 0.05 nats (the L6-archive smoke got
   r = 0.988 / 0.04).
3. **Consistency:** the trajectory checkpoint at the best epoch end matches the ES-selected
   model the script reports (same step; surprisal agrees).

## Outputs
- **HF:** same repo, `trajectory/seed{S}/rung24M/step-{NNNNNN}/` + a `trajectory/README.md`
  (grid, config, best/stop steps per seed, and a note that post-best checkpoints are past the
  val minimum, i.e. overfitting). One `upload_large_folder` per seed after its gates pass.
- **Oak:** `rsync` each seed's `trajectory/` to `oak-dtn:/oak/stanford/groups/mcfrank/childes-ladder-trajectory/`
  (no `--delete`; verify dry-run is empty before freeing `/data2`).
- **Repo:** `fits/llm/traj24m_surprisal.csv` (seed, step, epoch, word, surprisal) and a
  journal entry in `journal/experiments_llm.md`.

## Compute
3 runs × ~4.5 h ≈ **14 A40-hours**, 3 GPUs in parallel ≈ 5 h wall, plus ~1 h offline
surprisal and upload. Each run needs ~31 GB of GPU memory. As of 2026-10-09 all 8 A40s show
17–24 GB in use by other jobs, so check `nvidia-smi` for ≥ 3 GPUs under 8 GB before launch.
Launch detached (`setsid nohup`) with an atomic-claim worker modeled on `worker.sh`.
