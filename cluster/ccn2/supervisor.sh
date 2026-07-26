#!/bin/bash
# Robust supervisor for the composition-control run. Auto-queues after the CHILDES
# archive sweep, smoke-gates ONE run before fanning out, keeps 8 workers alive,
# reaps stale claims, writes a STATUS file, fully resumable (just relaunch it).
# Launch: cd /data2/mcfrank/ladder && setsid nohup bash supervisor.sh </dev/null >/dev/null 2>&1 &
set -u
ROOT=/data2/mcfrank/ladder; OUT=$ROOT/register_dev; mkdir -p "$OUT/claims"
GPUS=(0 1 2 3 4 5 6 7); NTOTAL=576
slog(){ echo "[$(date +%m-%dT%T)] $*" >> $OUT/supervisor.log; }
ndone(){ ls -d $OUT/*/UPLOADED 2>/dev/null|wc -l; }
ngave(){ ls -d $OUT/*/GAVEUP 2>/dev/null|wc -l; }
nclaim(){ ls -d $OUT/claims/* 2>/dev/null|wc -l; }
worker_alive(){ local pf=$OUT/gpu$1.pid; [ -f "$pf" ] && kill -0 "$(cat $pf 2>/dev/null)" 2>/dev/null; }
launch(){ cd $ROOT && setsid nohup bash register_worker.sh $1 </dev/null >>$OUT/worker_$1.log 2>&1 & slog "launched worker GPU$1"; }
status(){ echo "$(date) | UPLOADED $(ndone)/$NTOTAL | gaveup $(ngave) | claimed $(nclaim) | workers: $(for g in ${GPUS[@]}; do worker_alive $g && printf %s $g || printf .; done)" > $OUT/STATUS; }

slog "supervisor start pid $$"
# 1. wait for CHILDES archive sweep
slog "waiting for CHILDES sweep (180/180 or CHILDES_DONE)..."
while [ "$(ls -d $ROOT/retrain_dev/*/UPLOADED 2>/dev/null|wc -l)" -lt 180 ] && [ ! -f $ROOT/retrain_dev/CHILDES_DONE ]; do status; sleep 300; done
slog "CHILDES done -> composition control"
# 2. smoke-gate: one worker, require first UPLOADED before fan-out
if [ "$(ndone)" -eq 0 ]; then
  slog "SMOKE: GPU0 worker, awaiting first UPLOADED"; launch 0; t=0
  while [ "$(ndone)" -eq 0 ]; do
    sleep 60; t=$((t+1)); status
    [ "$(ls -d $OUT/*/FAILED_* 2>/dev/null|wc -l)" -ge 3 ] && { slog "SMOKE FAILED - HALT, needs attention"; touch $OUT/NEEDS_ATTENTION; exit 1; }
    [ "$t" -ge 40 ] && { slog "SMOKE timeout 40m - HALT"; touch $OUT/NEEDS_ATTENTION; exit 1; }
  done
  slog "SMOKE OK -> fan out"
fi
# 3+4. fan out + supervise
while true; do
  status; done=$(ndone); gave=$(ngave)
  [ $((done+gave)) -ge $NTOTAL ] && { slog "COMPLETE: $done uploaded, $gave gaveup"; touch $OUT/REGISTER_ALL_DONE; break; }
  # reap stale claims (idle >300m, not uploaded)
  now=$(date +%s)
  for c in $OUT/claims/*; do
    [ -d "$c" ] || continue; item=$(basename "$c")
    [ -f "$OUT/$item/UPLOADED" ] && { rmdir "$c" 2>/dev/null; continue; }
    ro=$OUT/$item/run.out; [ -f "$ro" ] || continue
    idle=$(( (now - $(stat -c %Y "$ro")) / 60 ))
    [ "$idle" -gt 300 ] && { slog "reap stale $item (idle ${idle}m)"; rm -rf "$OUT/$item"/checkpoint-* 2>/dev/null; rmdir "$c" 2>/dev/null; }
  done
  # restart dead workers only while unclaimed work remains
  unclaimed=$(( NTOTAL - done - gave - $(nclaim) ))
  if [ "$unclaimed" -gt 0 ]; then for g in ${GPUS[@]}; do worker_alive $g || { launch $g; sleep 5; }; done; fi
  sleep 300
done
slog "supervisor exit"; status
