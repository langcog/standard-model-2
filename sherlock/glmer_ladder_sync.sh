#!/usr/bin/env bash
## Pull ladder fit summaries (and optionally fit RDS) back from Sherlock.
## Wraps sync_from_remote.sh with the ladder/ pattern.
##
## Usage:
##   bash sherlock/glmer_ladder_sync.sh             # just summaries (small)
##   bash sherlock/glmer_ladder_sync.sh --with-fits # also pull fit RDS (large)

set -euo pipefail

cd "$(dirname "$0")/.."

WITH_FITS=0
if [ "${1:-}" = "--with-fits" ]; then
  WITH_FITS=1
fi

mkdir -p fits/glmer_ladder

# always pull the small summary CSVs
bash sherlock/sync_from_remote.sh 'ladder/summary_*.csv'

if [ "$WITH_FITS" = "1" ]; then
  bash sherlock/sync_from_remote.sh 'ladder/fit_*.rds'
fi

echo
echo "After sync, run:"
echo "  Rscript model/scripts/aggregate_ladder.R"
