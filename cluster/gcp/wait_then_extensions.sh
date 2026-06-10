#!/bin/bash
# Wait for the EN family to finish (specifically: long_no_freq_slopes_si_signed.rds
# appears on disk), then launch run_extensions.sh.
#
# Run in background:
#   nohup bash gcp/wait_then_extensions.sh > gcp_wait_extensions.log 2>&1 < /dev/null & disown

set -euo pipefail

cd "$HOME/standard_model_2"

SENTINEL="fits/long_no_freq_slopes_si_signed.rds"
echo "[$(date)] waiting for $SENTINEL ..."

# Poll every 5 min.
while [ ! -f "$SENTINEL" ]; do
  sleep 300
done

echo "[$(date)] $SENTINEL detected. Sleeping 60s to let extracts settle."
sleep 60

echo "[$(date)] launching gcp/run_extensions.sh ..."
bash gcp/run_extensions.sh

echo "[$(date)] extensions done."
