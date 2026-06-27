#!/bin/bash
# Warm-start continuation on the nonce corpus from one saved checkpoint. Shared LR across all arms.
# Usage: run_continue.sh GPU TAG INIT_CHECKPOINT   (env: LR, NPTS)
set -u
ROOT=/data2/mcfrank/ladder
FENG=$ROOT/standard-model-2/model/scripts/feng_eval
GPU=$1; TAG=$2; INIT=$3
RD=$ROOT/flatlr_grammar/nonce/runs/$TAG
mkdir -p "$RD"; cd $ROOT || exit 1
env HF_HUB_OFFLINE=1 CUDA_VISIBLE_DEVICES=$GPU PYTHONDONTWRITEBYTECODE=1 \
    PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  $ROOT/condaenv/bin/python -B $FENG/continue_train.py \
    --init_from "$ROOT/flatlr_grammar/runs/$INIT/final" \
    --train_file $ROOT/flatlr_grammar/nonce/continuation_train.txt \
    --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
    --output_dir "$RD" --nonce_probe $ROOT/flatlr_grammar/nonce/nonce_probe.jsonl \
    --nonce_csv "$RD/nonce.csv" --cdi_contexts_jsonl $ROOT/cdi_contexts.jsonl \
    --surprisal_csv "$RD/cdi.csv" --num_train_epochs 1 \
    --learning_rate "${LR:-1e-5}" --n_points "${NPTS:-30}" --seed 0 > "$RD/run.log" 2>&1
rc=$?; [ $rc -eq 0 ] && touch "$RD/DONE" || touch "$RD/FAILED_rc$rc"
