#!/usr/bin/env bash
#SBATCH --job-name=llmvar_grid
#SBATCH --output=%x_%a_%j.out
#SBATCH --error=%x_%a_%j.err
#SBATCH --partition=preempt
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -A marlowe-m000102
#
# Developmental-ladder GRID: 5 seeds x 11 nested CHILDES budgets, each trained
# to convergence via early-stopping (patience 4) with best-val model loaded, so
# the final surprisal eval = competence at convergence for that rung.
#
# array idx -> (seed, rung):  seed = SEEDS[idx / 11],  rung = RUNGS[idx % 11]
# Submit:  sbatch --array=0-54 train_ladder_grid.slurm
set -euo pipefail
ROOT=/scratch/m000102/mcfrank/llm_var_pilot
PY="$ROOT/condaenv/bin/python"
FENG="$ROOT/standard-model-2/studies/llm"

module load conda/24.3.0-0 cuda12.9/toolkit/12.9.1 cudnn/cuda12/9.3.0.75 2>/dev/null || true

SEEDS=(42 0 123 7 99)
RUNGS=(0.5M 1M 1.5M 2M 3M 4M 6M 8M 12M 16M 24M)
NR=${#RUNGS[@]}
idx="${SLURM_ARRAY_TASK_ID:?use --array=0-54}"
SEED="${SEEDS[$(( idx / NR ))]}"
TAG="${RUNGS[$(( idx % NR ))]}"

TRAIN_FILE="$ROOT/chunks/ladder_seed${SEED}/CHILDES_ladder_seed${SEED}_${TAG}.txt"
RUN_DIR="$ROOT/runs_grid/seed${SEED}_${TAG}"
SURPRISAL_CSV="$ROOT/runs_grid/surprisal_seed${SEED}_${TAG}.csv"
mkdir -p "$RUN_DIR" "$ROOT/runs_grid"
export HF_HOME="$ROOT/hf_cache"

echo "=== grid idx=$idx seed=$SEED budget=$TAG node=$(hostname) ==="
nvidia-smi -L || true
echo "train_file=$TRAIN_FILE"; echo "start=$(date +%s) ($(date))"

# Early-stopping converges each rung; ceiling of 40 epochs (small rungs need many,
# big rungs stop well before). n_log_points modest (best-val read from final eval).
$PY "$FENG/train_gpt2_childes.py" \
    --train_file "$TRAIN_FILE" \
    --val_file "$ROOT/TinyDialogues/data/CHILDES_val_ordered.txt" \
    --tokenizer_dir "$ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --config_file "$ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json" \
    --output_dir "$RUN_DIR" \
    --cdi_contexts_jsonl "$ROOT/cdi_contexts.jsonl" \
    --surprisal_csv "$SURPRISAL_CSV" \
    --num_train_epochs 40 --per_device_batch_size 8 --learning_rate 1e-4 \
    --seed "$SEED" --n_log_points 30 \
    --early_stopping_patience 4 --save_total_limit 1 \
    --eval_callback_batch_size 128 --eval_max_ctx 128 --eval_max_per_word 50
echo "end=$(date +%s) ($(date))"
# tidy: drop the saved checkpoint to save scratch (we only needed it for load-best)
rm -rf "$RUN_DIR"/checkpoint-* 2>/dev/null || true
