#!/bin/bash
# Disjoint-CHILDES-halves ladder (Option 1): tests whether "individuals converge
# with input" survives when individuals train on GENUINELY DISJOINT data.
# 2 pools (A,B = disjoint random halves of CHILDES) x 3 seeds x 8 rungs (0.5M..12M).
# Pool membership is fixed; seed varies the within-pool init+shuffle. A n B = empty
# at every budget -> no overlap confound.
#
# Usage: GPUS="0 2 3 4 5 6 7" bash run_disjoint_ladder_ccn2.sh
set -u
ROOT=/data2/mcfrank/ladder
PY=$ROOT/condaenv/bin/python
FENG=$ROOT/standard-model-2/model/scripts/feng_eval
export HF_HOME=$ROOT/hf_cache TMPDIR=$ROOT/.tmp
mkdir -p "$ROOT/runs_disjoint"

POOLS=(A B)
SEEDS=(42 0 123)
RUNGS=(0.5M 1M 2M 3M 4M 6M 8M 12M)
read -r -a GPUS <<< "${GPUS:-0 2 3 4 5 6 7}"
NG=${#GPUS[@]}

run_one() {  # $1=gpu $2=pool $3=seed $4=tag
  local gpu=$1 pool=$2 seed=$3 tag=$4
  local tf=$ROOT/chunks/pool${pool}_seed${seed}/CHILDES_ladder_seed${seed}_${tag}.txt
  local rd=$ROOT/runs_disjoint/pool${pool}_seed${seed}_${tag}
  local sc=$ROOT/runs_disjoint/surprisal_pool${pool}_seed${seed}_${tag}.csv
  local log=$ROOT/runs_disjoint/log_pool${pool}_seed${seed}_${tag}.out
  [ -f "$rd/DONE" ] && { echo "[$(date +%T)] skip $pool/$seed/$tag (DONE)" >> "$ROOT/runs_disjoint/dispatch.log"; return 0; }
  mkdir -p "$rd"
  echo "[$(date +%T)] GPU$gpu START pool$pool seed$seed $tag" >> "$ROOT/runs_disjoint/dispatch.log"
  CUDA_VISIBLE_DEVICES=$gpu $PY "$FENG/train_gpt2_childes.py" \
    --train_file "$tf" \
    --val_file "$ROOT/TinyDialogues/data/CHILDES_val_ordered.txt" \
    --tokenizer_dir "$ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --config_file "$ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json" \
    --output_dir "$rd" --cdi_contexts_jsonl "$ROOT/cdi_contexts.jsonl" \
    --surprisal_csv "$sc" \
    --num_train_epochs 40 --per_device_batch_size 8 --learning_rate 1e-4 \
    --seed "$seed" --n_log_points 3 --early_stopping_patience 4 --save_total_limit 1 \
    --max_eval_blocks 500 \
    --eval_callback_batch_size 128 --eval_max_ctx 128 --eval_max_per_word 50 \
    > "$log" 2>&1
  local rc=$?
  rm -rf "$rd"/checkpoint-* 2>/dev/null
  [ $rc -eq 0 ] && touch "$rd/DONE"
  echo "[$(date +%T)] GPU$gpu END   pool$pool seed$seed $tag (exit $rc)" >> "$ROOT/runs_disjoint/dispatch.log"
}

# work list: rung-outer so round-robin balances big/small across GPUs
WORK=()
for r in "${RUNGS[@]}"; do for p in "${POOLS[@]}"; do for s in "${SEEDS[@]}"; do WORK+=("$p:$s:$r"); done; done; done

for g in $(seq 0 $((NG-1))); do
( gpu=${GPUS[$g]}; i=$g
  while [ $i -lt ${#WORK[@]} ]; do
    IFS=':' read -r pool seed tag <<< "${WORK[$i]}"
    run_one "$gpu" "$pool" "$seed" "$tag"
    i=$((i+NG))
  done
) &
done
wait
echo "[$(date +%T)] DISJOINT_ALL_DONE" >> "$ROOT/runs_disjoint/dispatch.log"
echo DISJOINT_ALL_DONE
