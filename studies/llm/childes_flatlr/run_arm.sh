#!/bin/bash
# Launch ONE CHILDES flat-vs-decay LR arm. Robust on-box launcher (detached via setsid by caller).
# Usage: run_arm.sh GPU TAG TRAIN_FILE LR SCHED EPOCHS
# Env overrides: SEED NLOG NGRAM GMAXP NOBLIMP(=--no_blimp if set)
set -u
ROOT=/data2/mcfrank/ladder
PY=$ROOT/condaenv/bin/python
FENG=$ROOT/standard-model-2/model/scripts/feng_eval
ZORRO=$ROOT/flatlr_grammar/data/Zorro
GPU=$1; TAG=$2; TF=$3; LR=$4; SCHED=$5; EP=$6
RD=$ROOT/flatlr_grammar/runs/$TAG
mkdir -p "$RD"
cd $ROOT || exit 1
env HF_HOME=$ROOT/hf_cache CUDA_VISIBLE_DEVICES=$GPU PYTHONDONTWRITEBYTECODE=1 \
    PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  $PY -B $FENG/train_gpt2_childes_flatlr.py \
    --train_file "$TF" --val_file $ROOT/TinyDialogues/data/CHILDES_val_ordered.txt \
    --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
    --config_file $ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json \
    --output_dir "$RD" --cdi_contexts_jsonl $ROOT/cdi_contexts.jsonl \
    --surprisal_csv "$RD/cdi.csv" --grammar_csv "$RD/grammar.csv" \
    --num_train_epochs "$EP" --per_device_batch_size "${BSZ:-16}" --learning_rate "$LR" \
    --lr_scheduler_type "$SCHED" --seed "${SEED:-42}" \
    --n_log_points "${NLOG:-40}" --n_grammar_points "${NGRAM:-12}" \
    --grammar_max_pairs "${GMAXP:-150}" --zorro_dir $ZORRO \
    --max_eval_blocks 500 ${NOBLIMP:+--no_blimp} > "$RD/run.log" 2>&1
rc=$?
[ $rc -eq 0 ] && touch "$RD/DONE" || touch "$RD/FAILED_rc$rc"
