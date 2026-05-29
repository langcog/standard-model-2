#!/usr/bin/env bash
# Unattended weekend run of the two gamma (input-on-slope) variants of the
# pooled IO model on this laptop (M2 Max: 8 P-cores). Variants run
# SEQUENTIALLY, each as 4 chains x 2 threads (= 8 P-cores), so RAM stays low
# and each fit gets the full machine. caffeinate keeps the Mac awake.
#
# Launch detached (survives terminal close):
#   nohup ./run_gamma_weekend.sh > logs/gamma_weekend.log 2>&1 &
#
# Keep the laptop PLUGGED IN. Check progress: tail -f logs/io_pooled_gamma_*.log
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs

export STAN_CHAINS=4
export STAN_THREADS_PER_CHAIN=2
export STAN_WARMUP=500
export STAN_ITER=1000

for v in add mult; do
  echo "===== gamma-$v START $(date) ====="
  caffeinate -dimsu Rscript model/scripts/fit_io_pooled_gamma.R "$v" \
    > "logs/io_pooled_gamma_${v}.log" 2>&1
  echo "===== gamma-$v END $(date) (exit $?) ====="
done
echo "ALL DONE $(date)"
