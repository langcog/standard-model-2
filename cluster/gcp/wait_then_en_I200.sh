#!/bin/bash
# Wait for the NO warmstart fit pipeline to finish, then launch a
# fresh M_best fit on the I=200 J=198 English bundle. Quick to run
# (N=145K, ~10-15 min/chain expected), so it fits comfortably after
# the NO refit and before VM shutdown.
#
# Usage: nohup bash gcp/wait_then_en_I200.sh > wait_en_I200.log 2>&1 &

set -uo pipefail
cd "$HOME/standard_model_2"
export STANDARD_MODEL_ROOT="$PWD"

LOG="$PWD/wait_en_I200.log"
echo "==== wait_then_en_I200 started $(date -u) ====" | tee -a "$LOG"

# Poll for the warmstart pipeline to finish.
while ! grep -q 'wait_then_no_warmstart done' "$PWD/wait_no_warmstart.log" 2>/dev/null; do
  sleep 60
done
echo "NO warmstart pipeline done at $(date -u). Launching EN I=200." | tee -a "$LOG"

rm -f fits/long_no_freq_slopes_english_I200.rds
rm -rf fits/csvs_long_no_freq_slopes_english_I200

export STAN_ITER=2000
export STAN_WARMUP=1000
export STAN_ADAPT_DELTA=0.95
export STAN_THREADS_PER_CHAIN=16
export STAN_MAX_TREEDEPTH=10

echo "==== $(date -u): running fit_longitudinal.R long_no_freq_slopes english_I200" | tee -a "$LOG"
time Rscript model/scripts/fit_longitudinal.R long_no_freq_slopes english_I200 2>&1 | tee -a "$LOG"

TAG=long_no_freq_slopes_english_I200
if [ -f "fits/${TAG}.rds" ]; then
  echo "==== $(date -u): .rds present, running extracts" | tee -a "$LOG"
  Rscript sherlock/extract_summary_table_only.R "$TAG" 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_scalar_draws.R     "$TAG" 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_delta_j_slim.R     "$TAG" 2>&1 | tee -a "$LOG"
else
  echo "==== $(date -u): no .rds, recovering from CSVs" | tee -a "$LOG"
  Rscript sherlock/recover_from_csvs.R "$TAG" 2>&1 | tee -a "$LOG"
fi

echo "==== wait_then_en_I200 done $(date -u) ====" | tee -a "$LOG"
