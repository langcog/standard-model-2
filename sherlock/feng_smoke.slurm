#!/usr/bin/env bash
#SBATCH --job-name=feng_smoke
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --constraint=GPU_MEM:24GB|GPU_MEM:32GB|GPU_MEM:40GB|GPU_MEM:48GB|GPU_MEM:80GB
#
# Smoke test: run a tiny ~2-3 minute training on a small CHILDES subset to
# verify the full pipeline (HF Trainer + WordSurprisalCallback + per-word CSV
# logging) works end-to-end on a GPU.

set -euo pipefail

FENG_ROOT="${FENG_ROOT:-/scratch/users/mcfrank/feng_eval}"
SCRIPTS_DIR="$FENG_ROOT/scripts"

# Modules + venv
source /share/software/user/open/lmod/lmod/init/bash
ml reset >/dev/null 2>&1 || true
ml load python/3.12.1
ml load cuda/12.6.1 cudnn/9.4.0
source "$FENG_ROOT/venv/bin/activate"

export HF_HOME="$FENG_ROOT/hf_cache"
export TRANSFORMERS_CACHE="$HF_HOME/transformers"
export HF_DATASETS_CACHE="$HF_HOME/datasets"
mkdir -p "$HF_HOME"

# Build mini train/val by taking the first lines.
DATA_DIR="$FENG_ROOT/TinyDialogues/data/CHILDES"
MINI_DIR="$FENG_ROOT/smoke_mini"
mkdir -p "$MINI_DIR"
head -n 400 "$DATA_DIR/CHILDES_train_ordered.txt" > "$MINI_DIR/train.txt"
head -n 100 "$DATA_DIR/CHILDES_val_ordered.txt" > "$MINI_DIR/val.txt"

# Ensure CDI contexts JSONL exists (cheap; recompute if missing)
if [ ! -s "$FENG_ROOT/cdi_contexts.jsonl" ]; then
  echo "=== building cdi_contexts.jsonl ==="
  python "$SCRIPTS_DIR/extract_cdi_contexts.py" \
    --tokenizer "$FENG_ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --val_file "$DATA_DIR/CHILDES_val_ordered.txt" \
    --cdi_words "$FENG_ROOT/cdi_words.txt" \
    --out_jsonl "$FENG_ROOT/cdi_contexts.jsonl" \
    --out_coverage "$FENG_ROOT/cdi_coverage.csv" \
    --max_per_word 200 --max_ctx 1024
fi

OUT_DIR="$FENG_ROOT/smoke_out"
SURPRISAL_CSV="$FENG_ROOT/smoke_surprisal.csv"
rm -rf "$OUT_DIR" "$SURPRISAL_CSV"

echo "=== running smoke training ==="
echo "GPU: $(nvidia-smi -L 2>/dev/null || echo none)"

python "$SCRIPTS_DIR/train_gpt2_childes.py" \
    --train_file "$MINI_DIR/train.txt" \
    --val_file "$MINI_DIR/val.txt" \
    --tokenizer_dir "$FENG_ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --config_file "$FENG_ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json" \
    --output_dir "$OUT_DIR" \
    --cdi_contexts_jsonl "$FENG_ROOT/cdi_contexts.jsonl" \
    --surprisal_csv "$SURPRISAL_CSV" \
    --num_train_epochs 2 \
    --per_device_batch_size 4 \
    --learning_rate 1e-4 \
    --seed 42 \
    --n_log_points 8 \
    --eval_callback_batch_size 32 \
    --no_save

echo "=== smoke test output ==="
ls -la "$SURPRISAL_CSV"
echo
echo "first 5 + last 5 surprisal rows:"
head -5 "$SURPRISAL_CSV"
echo "..."
tail -5 "$SURPRISAL_CSV"
echo
echo "unique steps logged: $(tail -n +2 "$SURPRISAL_CSV" | cut -d, -f1 | sort -u | wc -l)"
echo "unique words logged: $(tail -n +2 "$SURPRISAL_CSV" | cut -d, -f3 | sort -u | wc -l)"
