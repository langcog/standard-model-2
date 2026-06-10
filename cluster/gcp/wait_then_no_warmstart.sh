#!/bin/bash
# Wait for the IO/Peekbank queue to finish, then refit NO M_best with
# warm-start from the previous run's summary medians. Iter/warmup
# unchanged (2000/1000) -- warm-start should resolve the slow mixing
# on sigma_alpha (Rhat=1.23, ess=12 in the cold-start fit).
#
# Usage: nohup bash gcp/wait_then_no_warmstart.sh > wait_no_warmstart.log 2>&1 &

set -uo pipefail
cd "$HOME/standard_model_2"
export STANDARD_MODEL_ROOT="$PWD"

LOG="$PWD/wait_no_warmstart.log"
echo "==== wait_then_no_warmstart started $(date -u) ====" | tee -a "$LOG"

# Poll for the queue to finish. Detect the "all done" marker in its log
# (more reliable than process-existence, because the watcher itself may
# have died between sessions).
while ! grep -q 'queue_io_proc all done' "$PWD/queue_io_proc.log" 2>/dev/null; do
  sleep 60
done
echo "IO/Peekbank queue done at $(date -u). Launching NO refit." | tee -a "$LOG"

# Force a fresh fit: remove any cached .rds + persistent CSV dir from
# the previous failed run.
rm -f fits/long_no_freq_slopes_norwegian.rds
rm -rf fits/csvs_long_no_freq_slopes_norwegian

# Warm-start from the recovered scalar summary medians (built by
# build_init_from_summary in helpers.R). Same iter/warmup as before.
export STAN_INIT_FROM="long_no_freq_slopes_norwegian"
export STAN_ITER=2000
export STAN_WARMUP=1000
export STAN_ADAPT_DELTA=0.95
export STAN_THREADS_PER_CHAIN=16
export STAN_MAX_TREEDEPTH=10

echo "==== $(date -u): running fit_longitudinal.R long_no_freq_slopes norwegian (warm-start)" | tee -a "$LOG"
time Rscript model/scripts/fit_longitudinal.R long_no_freq_slopes norwegian 2>&1 | tee -a "$LOG"

# If save_object succeeded, run the normal extracts.
# If it failed (no .rds), fall back to recovery from persistent CSVs.
if [ -f fits/long_no_freq_slopes_norwegian.rds ]; then
  echo "==== $(date -u): .rds present, running extracts" | tee -a "$LOG"
  Rscript sherlock/extract_summary_table_only.R long_no_freq_slopes_norwegian 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_scalar_draws.R     long_no_freq_slopes_norwegian 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_delta_j_slim.R     long_no_freq_slopes_norwegian 2>&1 | tee -a "$LOG"
else
  echo "==== $(date -u): no .rds (save_object failed), recovering from CSVs" | tee -a "$LOG"
  Rscript sherlock/recover_from_csvs.R long_no_freq_slopes_norwegian 2>&1 | tee -a "$LOG"
fi

echo "==== wait_then_no_warmstart done $(date -u) ====" | tee -a "$LOG"
