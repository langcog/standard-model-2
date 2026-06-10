#!/usr/bin/env bash
#SBATCH --job-name=feng_gpt2
#SBATCH --output=%x_seed%a_%j.out
#SBATCH --error=%x_seed%a_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --constraint=GPU_GEN:AMP|GPU_GEN:LOV|GPU_GEN:HPR
# Ampere (RTX 3090, A100) / Lovelace (L40S) / Hopper (H100). Excludes Volta
# (V100) and older — V100 was ~3x slower per step on this workload and could
# not finish in the 24h wall.
#
# Train one seed of GPT-2 small on CHILDES, logging per-CDI-word surprisal at
# ~80 log-spaced training steps. Seed is set per array task.
#
# Submit:
#   sbatch --array=42,0,123 sherlock/feng_train_gpt2.slurm
# or with a custom output prefix:
#   sbatch --array=42 sherlock/feng_train_gpt2.slurm
#
# Required env (set in script): FENG_ROOT, FENG_OUT
set -euo pipefail

FENG_ROOT="${FENG_ROOT:-/scratch/users/mcfrank/feng_eval}"
SCRIPTS_DIR="${SCRIPTS_DIR:-$FENG_ROOT/scripts}"
FENG_OUT="${FENG_OUT:-$FENG_ROOT/runs}"

SEED="${SLURM_ARRAY_TASK_ID:-42}"
RUN_DIR="$FENG_OUT/gpt2_childes_seed${SEED}"
SURPRISAL_CSV="$FENG_OUT/surprisal_gpt2_childes_seed${SEED}.csv"
mkdir -p "$RUN_DIR" "$FENG_OUT"

cd "$FENG_ROOT/TinyDialogues"

# Modules + venv
source /share/software/user/open/lmod/lmod/init/bash
ml reset >/dev/null 2>&1 || true
ml load python/3.12.1
ml load cuda/12.6.1 cudnn/9.4.0
source "$FENG_ROOT/venv/bin/activate"

# HF cache to scratch
export HF_HOME="$FENG_ROOT/hf_cache"
export TRANSFORMERS_CACHE="$HF_HOME/transformers"
export HF_DATASETS_CACHE="$HF_HOME/datasets"
mkdir -p "$HF_HOME"

echo "=== seed=$SEED ==="
echo "node: $(hostname)"
echo "GPU: $(nvidia-smi -L 2>/dev/null || echo none)"
echo "start: $(date)"

python "$SCRIPTS_DIR/train_gpt2_childes.py" \
    --train_file "$FENG_ROOT/TinyDialogues/data/CHILDES/CHILDES_train_ordered.txt" \
    --val_file "$FENG_ROOT/TinyDialogues/data/CHILDES/CHILDES_val_ordered.txt" \
    --tokenizer_dir "$FENG_ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --config_file "$FENG_ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json" \
    --output_dir "$RUN_DIR" \
    --cdi_contexts_jsonl "$FENG_ROOT/cdi_contexts.jsonl" \
    --surprisal_csv "$SURPRISAL_CSV" \
    --num_train_epochs 20 \
    --per_device_batch_size 8 \
    --learning_rate 1e-4 \
    --seed "$SEED" \
    --n_log_points 80 \
    --eval_callback_batch_size 128 \
    --eval_max_ctx 128 \
    --eval_max_per_word 50 \
    --no_save

echo "end: $(date)"
