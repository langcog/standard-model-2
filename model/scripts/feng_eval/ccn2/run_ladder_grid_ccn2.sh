#!/bin/bash
# Developmental-ladder GRID on the ccn2 A40 box (NO SLURM): round-robin 55 runs
# (5 seeds x 11 nested rungs) across 4 free GPUs, each rung early-stopped to
# convergence (best-val loaded). Completion markers ($rundir/DONE) make it
# resumable and let the smoke run count as grid run #1.
#
# Usage:
#   GPUS="2 3 6 7" bash run_ladder_grid_ccn2.sh          # full grid
#   ONLY="42:0.5M" GPUS="2" bash run_ladder_grid_ccn2.sh # smoke one run
set -u
ROOT=/data2/mcfrank/ladder
PY=$ROOT/condaenv/bin/python
FENG=$ROOT/standard-model-2/model/scripts/feng_eval
export HF_HOME=$ROOT/hf_cache TMPDIR=$ROOT/.tmp
mkdir -p "$ROOT/runs_grid"

SEEDS=(42 0 123 7 99)
RUNGS=(0.5M 1M 1.5M 2M 3M 4M 6M 8M 12M 16M 24M)
read -r -a GPUS <<< "${GPUS:-2 3 6 7}"
NG=${#GPUS[@]}

run_one() {  # $1=gpu $2=seed $3=tag
  local gpu=$1 seed=$2 tag=$3
  local tf=$ROOT/chunks/ladder_seed${seed}/CHILDES_ladder_seed${seed}_${tag}.txt
  local rd=$ROOT/runs_grid/seed${seed}_${tag}
  local sc=$ROOT/runs_grid/surprisal_seed${seed}_${tag}.csv
  local log=$ROOT/runs_grid/log_seed${seed}_${tag}.out
  [ -f "$rd/DONE" ] && { echo "[$(date +%T)] skip seed$seed $tag (DONE)" >> "$ROOT/runs_grid/dispatch.log"; return 0; }
  mkdir -p "$rd"
  echo "[$(date +%T)] GPU$gpu START seed$seed $tag" >> "$ROOT/runs_grid/dispatch.log"
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
  if [ $rc -eq 0 ]; then touch "$rd/DONE"; fi
  echo "[$(date +%T)] GPU$gpu END   seed$seed $tag (exit $rc)" >> "$ROOT/runs_grid/dispatch.log"
}

# build work list: rung-outer so round-robin balances big/small across GPUs
WORK=()
for r in "${RUNGS[@]}"; do for s in "${SEEDS[@]}"; do WORK+=("$s:$r"); done; done

# single-run mode (smoke)
if [ -n "${ONLY:-}" ]; then
  IFS=':' read -r s t <<< "$ONLY"; run_one "${GPUS[0]}" "$s" "$t"; exit $?
fi

# launch one sequential worker per GPU; worker g takes items g, g+NG, g+2NG, ...
for g in $(seq 0 $((NG-1))); do
( gpu=${GPUS[$g]}; i=$g
  while [ $i -lt ${#WORK[@]} ]; do
    IFS=':' read -r seed tag <<< "${WORK[$i]}"
    run_one "$gpu" "$seed" "$tag"
    i=$((i+NG))
  done
) &
done
wait
echo "[$(date +%T)] GRID_ALL_DONE" >> "$ROOT/runs_grid/dispatch.log"
echo GRID_ALL_DONE
