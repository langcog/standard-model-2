#!/bin/bash
# Fan CHILDES flat-vs-decay arms across GPUs, one per GPU with a per-GPU SEQUENTIAL queue
# (so #arms > #GPUs is fine; no double-booking). Detach the whole thing via setsid.
# Usage:  GPUS_LIST="0 1 2 3 4 5" dispatch_arms.sh MANIFEST
# Manifest line: tag|train_file|lr|sched|epochs|seed   (# comments ok)
set -u
ROOT=/data2/mcfrank/ladder
MANIFEST=$1
mapfile -t ARMS < <(grep -vE '^#|^[[:space:]]*$' "$MANIFEST")
read -r -a GPUS <<< "${GPUS_LIST:-$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | awk '$2<8000{print $1}' | tr '\n' ' ')}"
NG=${#GPUS[@]}
echo "[dispatch] $(date) arms=${#ARMS[@]} gpus=(${GPUS[*]})"
for g in $(seq 0 $((NG-1))); do
( gpu=${GPUS[$g]}; idx=$g
  while [ $idx -lt ${#ARMS[@]} ]; do
    IFS='|' read -r tag tf lr sched ep seed <<< "${ARMS[$idx]}"
    RD=$ROOT/flatlr_grammar/runs/$tag
    if [ ! -f "$RD/DONE" ]; then
      echo "[dispatch] GPU$gpu START $tag (lr=$lr sched=$sched ep=$ep seed=$seed) $(date +%T)"
      SEED=$seed NLOG=40 NGRAM=12 GMAXP=150 BSZ=16 \
        bash $ROOT/flatlr_grammar/run_arm.sh "$gpu" "$tag" "$tf" "$lr" "$sched" "$ep"
      echo "[dispatch] GPU$gpu END   $tag $(date +%T)"
    fi
    idx=$((idx+NG))
  done
) &
done
wait
echo "[dispatch] ALL DONE $(date)"
