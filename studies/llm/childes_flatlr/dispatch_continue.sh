#!/bin/bash
# Fan continuation runs across GPUs (per-GPU sequential queue). Detach via setsid.
# Usage: GPUS_LIST="0 1 2 3 4 5" dispatch_continue.sh MANIFEST   (manifest line: TAG|INIT_CHECKPOINT)
set -u
ROOT=/data2/mcfrank/ladder
MANIFEST=$1
mapfile -t ARMS < <(grep -vE '^#|^[[:space:]]*$' "$MANIFEST")
read -r -a GPUS <<< "${GPUS_LIST:-$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | awk '$2<8000{print $1}' | tr '\n' ' ')}"
NG=${#GPUS[@]}
echo "[cont-dispatch] $(date) arms=${#ARMS[@]} gpus=(${GPUS[*]})"
for g in $(seq 0 $((NG-1))); do
( gpu=${GPUS[$g]}; idx=$g
  while [ $idx -lt ${#ARMS[@]} ]; do
    IFS='|' read -r tag init <<< "${ARMS[$idx]}"
    RD=$ROOT/flatlr_grammar/nonce/runs/$tag
    if [ ! -f "$RD/DONE" ]; then
      echo "[cont-dispatch] GPU$gpu START $tag <- $init $(date +%T)"
      LR=1e-5 NPTS=30 bash $ROOT/flatlr_grammar/run_continue.sh "$gpu" "$tag" "$init"
      echo "[cont-dispatch] GPU$gpu END   $tag $(date +%T)"
    fi
    idx=$((idx+NG))
  done
) &
done
wait
echo "[cont-dispatch] ALL DONE $(date)"
