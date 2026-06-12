#!/bin/bash
# Finer ladder: existing 5 seeds get 7 INSERT rungs (combine w/ their 11 in runs_grid);
# 5 NEW seeds get all 18 rungs. Goal: clean per-word kappa distribution (18 budgets)
# + sigma_kappa from 10 seeds. Output -> runs_finer/.
set -u
ROOT=/data2/mcfrank/ladder; PY=$ROOT/condaenv/bin/python
FENG=$ROOT/standard-model-2/model/scripts/feng_eval
export HF_HOME=$ROOT/hf_cache TMPDIR=$ROOT/.tmp
mkdir -p $ROOT/runs_finer
EXIST=(42 0 123 7 99); INSERTS=(0.75M 1.25M 1.75M 2.5M 3.5M 5M 7M)
NEW=(1 2 3 4 5); FULL=(0.5M 0.75M 1M 1.25M 1.5M 1.75M 2M 2.5M 3M 3.5M 4M 5M 6M 7M 8M 12M 16M 24M)
read -r -a GPUS <<< "${GPUS:-0 1 2 3 4 5 6 7}"; NG=${#GPUS[@]}

run_one() {  # gpu seed tag
  local gpu=$1 seed=$2 tag=$3
  local tf=$ROOT/chunks/ladder_seed${seed}/CHILDES_ladder_seed${seed}_${tag}.txt
  local rd=$ROOT/runs_finer/seed${seed}_${tag}
  local sc=$ROOT/runs_finer/surprisal_seed${seed}_${tag}.csv
  [ -f "$rd/DONE" ] && { echo "[$(date +%T)] skip $seed/$tag" >> $ROOT/runs_finer/dispatch.log; return 0; }
  mkdir -p "$rd"
  echo "[$(date +%T)] GPU$gpu START seed$seed $tag" >> $ROOT/runs_finer/dispatch.log
  CUDA_VISIBLE_DEVICES=$gpu $PY $FENG/train_gpt2_childes.py \
    --train_file "$tf" --val_file $ROOT/TinyDialogues/data/CHILDES_val_ordered.txt \
    --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
    --config_file $ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json \
    --output_dir "$rd" --cdi_contexts_jsonl $ROOT/cdi_contexts.jsonl --surprisal_csv "$sc" \
    --num_train_epochs 40 --per_device_batch_size 8 --learning_rate 1e-4 --seed $seed \
    --n_log_points 3 --early_stopping_patience 4 --save_total_limit 1 --max_eval_blocks 500 \
    --eval_callback_batch_size 128 --eval_max_ctx 128 --eval_max_per_word 50 > "$rd/run.out" 2>&1
  local rc=$?; rm -rf "$rd"/checkpoint-* 2>/dev/null; [ $rc -eq 0 ] && touch "$rd/DONE"
  echo "[$(date +%T)] GPU$gpu END   seed$seed $tag (exit $rc)" >> $ROOT/runs_finer/dispatch.log
}
WORK=()
for r in "${INSERTS[@]}"; do for s in "${EXIST[@]}"; do WORK+=("$s:$r"); done; done
for r in "${FULL[@]}"; do for s in "${NEW[@]}"; do WORK+=("$s:$r"); done; done
echo "total work items: ${#WORK[@]}" >> $ROOT/runs_finer/dispatch.log
for g in $(seq 0 $((NG-1))); do
( gpu=${GPUS[$g]}; i=$g
  while [ $i -lt ${#WORK[@]} ]; do IFS=':' read -r s t <<< "${WORK[$i]}"; run_one "$gpu" "$s" "$t"; i=$((i+NG)); done
) &
done
wait
echo "[$(date +%T)] FINER_ALL_DONE" >> $ROOT/runs_finer/dispatch.log
