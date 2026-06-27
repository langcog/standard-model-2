#!/bin/bash
# On-box sequencer: wait until all 8 1x arms are SETTLED (DONE or FAILED), then launch the
# 4x sweep on whatever GPUs are free, then mark ALL_DONE. Fully detached (setsid) so it
# survives the laptop closing / ssh disconnect. Owns the 4x launch (local controller must NOT).
set -u
ROOT=/data2/mcfrank/ladder/flatlr_grammar
cd $ROOT
ARMS_1X="std_s42_1x flat1e4_s42_1x flat3e5_s42_1x flat1e5_s42_1x std_s1_1x flat1e4_s1_1x flat3e5_s1_1x flat1e5_s1_1x"
echo "[seq] $(date) waiting for 1x to settle" >> sequence_4x.log
while :; do
  s=0
  for t in $ARMS_1X; do { [ -f runs/$t/DONE ] || ls runs/$t/FAILED_rc* >/dev/null 2>&1; } && s=$((s+1)); done
  [ $s -ge 8 ] && break
  sleep 180
done
echo "[seq] $(date) 1x settled ($s/8); DONE=$(ls runs/*_1x/DONE 2>/dev/null|wc -l). launching 4x" >> sequence_4x.log
mkdir -p runs/std_s42_4x   # claim, so any stray local controller skips its 4x launch
# no GPUS_LIST -> dispatch_arms.sh auto-picks free GPUs (<8000MB). 4x = 4 arms, 60 epochs.
bash dispatch_arms.sh manifest_4x.txt >> dispatch_4x.log 2>&1
echo "[seq] $(date) 4x dispatcher returned" >> sequence_4x.log
touch ALL_DONE
echo "[seq] $(date) ALL_DONE" >> sequence_4x.log
