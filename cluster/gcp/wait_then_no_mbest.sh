#!/bin/bash
# Wait for the EN M_best fit (long_no_freq_slopes english) to finish
# its full pipeline (sampling + extracts), then launch the equivalent
# fit on the Norwegian bundle. Lets the compute node stay busy without
# manual intervention.
#
# Usage: nohup bash gcp/wait_then_no_mbest.sh > wait_no_mbest.log 2>&1 &

set -uo pipefail

cd "$HOME/standard_model_2"

LOG="$PWD/wait_no_mbest.log"

echo "==== wait_then_no_mbest: started $(date -u) ====" | tee -a "$LOG"

# Poll for the EN run_fit.sh script to exit. Looks for the whole
# pipeline (sampling + extracts), not just the sampling R process,
# since the extracts run sequentially after the sample.
while pgrep -af 'run_fit.sh long_no_freq_slopes english' \
        | grep -v 'wait_then_no_mbest\|grep' > /dev/null; do
  sleep 60
done

echo "EN fit pipeline finished at $(date -u). Launching NO." | tee -a "$LOG"

# Match the env used for EN.
exec env STAN_ITER=2000 STAN_WARMUP=1000 STAN_ADAPT_DELTA=0.95 \
         STAN_THREADS_PER_CHAIN=16 STAN_MAX_TREEDEPTH=10 \
     bash gcp/run_fit.sh long_no_freq_slopes norwegian
