#!/usr/bin/env bash
#SBATCH --job-name=llmvar_pilot
#SBATCH --output=%x_chunk%a_%j.out
#SBATCH --error=%x_chunk%a_%j.err
#SBATCH --partition=preempt
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -A marlowe-PROJECTID          # <-- EDIT: your Marlowe account (-A marlowe-xxxx)
#
# Data-variance pilot: train GPT-2 small on two DISJOINT 10M-word CHILDES chunks,
# SAME seed, logging per-CDI-word surprisal at ~80 log-spaced steps. Submit as a
# 2-task array, one task per chunk:
#
#   sbatch --array=0,1 train_gpt2_childes.slurm
#
# Array task 0 -> chunk A, task 1 -> chunk B. Seed fixed (data is the only IV).
#
# PREEMPT partition notes (12h wall, can be preempted within 15 min):
#   - A 10M-word run is ~2h on an H100, so it usually completes inside one
#     allocation well under the 12h cap. Default mode is --no_save (matches the
#     original Sherlock run): if a task is preempted, just resubmit that one
#     array index. The surprisal CSV is append-mode, so a resubmit that restarts
#     from step 0 self-heals via the callback's header/dedup handling.
#   - For full robustness against repeated preemption, switch to checkpoint+
#     resume: set ROBUST=1 below. NOTE this requires the small --resume_if_present
#     patch to train_gpt2_childes.py documented in outputs/marlowe_pilot_runbook.md.
set -euo pipefail

PILOT_ROOT="${PILOT_ROOT:?set PILOT_ROOT to your Marlowe scratch path}"
FENG="$PILOT_ROOT/standard-model-2/studies/llm"
ROBUST="${ROBUST:-0}"

# ---- map array task -> chunk ----
TASK="${SLURM_ARRAY_TASK_ID:-0}"
case "$TASK" in
  0) CHUNK=A ;;
  1) CHUNK=B ;;
  *) echo "unknown array task $TASK (use --array=0,1)"; exit 1 ;;
esac
TRAIN_FILE="$PILOT_ROOT/chunks/CHILDES_chunk${CHUNK}.txt"
SEED=42                                # fixed across chunks: data is the only IV

RUN_DIR="$PILOT_ROOT/runs/gpt2_childes_chunk${CHUNK}_seed${SEED}"
SURPRISAL_CSV="$PILOT_ROOT/runs/surprisal_gpt2_childes_chunk${CHUNK}_seed${SEED}.csv"
mkdir -p "$RUN_DIR" "$PILOT_ROOT/runs"

# ---- modules + env (EDIT module names to match Marlowe) ----
module load slurm 2>/dev/null || true
module load nvhpc 2>/dev/null || true
module load cudnn/cuda12/9.3.0.75 2>/dev/null || true
module load python/3.12 2>/dev/null || module load python 2>/dev/null || true
source "$PILOT_ROOT/venv/bin/activate"
export HF_HOME="$PILOT_ROOT/hf_cache"; mkdir -p "$HF_HOME"

echo "=== chunk=$CHUNK seed=$SEED ==="
echo "node: $(hostname)"; echo "GPU: $(nvidia-smi -L 2>/dev/null || echo none)"
echo "train_file: $TRAIN_FILE"; echo "start: $(date)"

# ---- save/resume flags ----
# Default: --no_save (proven Sherlock path). Set ROBUST=1 only AFTER applying the
# --resume_if_present driver patch (see runbook) to survive repeated preemption.
SAVE_ARGS="--no_save"
if [ "$ROBUST" = "1" ]; then
  SAVE_ARGS="--save_strategy steps --save_steps 5000 --save_total_limit 2 --resume_if_present"
fi

python "$FENG/train_gpt2_childes.py" \
    --train_file "$TRAIN_FILE" \
    --val_file "$PILOT_ROOT/TinyDialogues/data/CHILDES/CHILDES_val_ordered.txt" \
    --tokenizer_dir "$PILOT_ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --config_file "$PILOT_ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json" \
    --output_dir "$RUN_DIR" \
    --cdi_contexts_jsonl "$PILOT_ROOT/cdi_contexts.jsonl" \
    --surprisal_csv "$SURPRISAL_CSV" \
    --num_train_epochs 20 \
    --per_device_batch_size 8 \
    --learning_rate 1e-4 \
    --seed "$SEED" \
    --n_log_points 80 \
    --eval_callback_batch_size 128 \
    --eval_max_ctx 128 \
    --eval_max_per_word 50 \
    $SAVE_ARGS

echo "end: $(date)"
