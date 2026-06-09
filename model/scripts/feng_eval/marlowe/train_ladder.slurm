#!/usr/bin/env bash
#SBATCH --job-name=llmvar_ladder
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
# Developmental ladder, ONE seed: GPT-2 small trained on nested CHILDES budgets
# [1,2,4,8,16,24]M words. Per-budget epoch caps keep every run <=~90 min so it
# fits preempt's 4h cap under hero-run contention. eval_loss is logged per epoch
# (read competence at the val-loss minimum post-hoc).
#
# Submit:  sbatch --array=0-5 train_ladder.slurm
set -euo pipefail
ROOT=/scratch/m000102/mcfrank/llm_var_pilot
PY="$ROOT/condaenv/bin/python"
FENG="$ROOT/standard-model-2/model/scripts/feng_eval"
SEED="${LADDER_SEED:-42}"

module load conda/24.3.0-0 cuda12.9/toolkit/12.9.1 cudnn/cuda12/9.3.0.75 2>/dev/null || true

# array index -> (budget tag, epochs). Epochs capped so big budgets stay <~90 min.
TAGS=(1M 2M 4M 8M 16M 24M)
EPOCHS=(20 20 20 20 15 10)
i="${SLURM_ARRAY_TASK_ID:?use --array=0-5}"
TAG="${TAGS[$i]}"; EP="${EPOCHS[$i]}"

TRAIN_FILE="$ROOT/chunks/ladder_seed${SEED}/CHILDES_ladder_seed${SEED}_${TAG}.txt"
RUN_DIR="$ROOT/runs_ladder/seed${SEED}_${TAG}"
SURPRISAL_CSV="$ROOT/runs_ladder/surprisal_seed${SEED}_${TAG}.csv"
mkdir -p "$RUN_DIR" "$ROOT/runs_ladder"
export HF_HOME="$ROOT/hf_cache"

echo "=== ladder seed=$SEED budget=$TAG epochs=$EP node=$(hostname) ==="
nvidia-smi -L || true
echo "train_file=$TRAIN_FILE"; echo "start=$(date +%s) ($(date))"

$PY "$FENG/train_gpt2_childes.py" \
    --train_file "$TRAIN_FILE" \
    --val_file "$ROOT/TinyDialogues/data/CHILDES_val_ordered.txt" \
    --tokenizer_dir "$ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --config_file "$ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json" \
    --output_dir "$RUN_DIR" \
    --cdi_contexts_jsonl "$ROOT/cdi_contexts.jsonl" \
    --surprisal_csv "$SURPRISAL_CSV" \
    --num_train_epochs "$EP" --per_device_batch_size 8 --learning_rate 1e-4 \
    --seed "$SEED" --n_log_points 80 \
    --eval_callback_batch_size 128 --eval_max_ctx 128 --eval_max_per_word 50 \
    --no_save
echo "end=$(date +%s) ($(date))"
