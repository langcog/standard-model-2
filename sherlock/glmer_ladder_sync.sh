#!/usr/bin/env bash
## Pull glmer_ladder fit outputs back from Sherlock.
## Wraps sync_from_remote.sh with the glmer_ladder/ pattern.
##
## Usage:
##   bash sherlock/glmer_ladder_sync.sh               # summaries only (small)
##   bash sherlock/glmer_ladder_sync.sh --with-fits   # also pull fit RDS (big; needed for ranef analysis)
##
## Fit RDS files are essential if you want to extract per-child random
## effects later (xi_i and zeta_i BLUPs for demographic regressions).
## They are also large — ~tens of MB to a couple of GB total across 49
## cells, depending on language size — so we keep them as an opt-in flag.

set -euo pipefail

cd "$(dirname "$0")/.."

WITH_FITS=0
if [ "${1:-}" = "--with-fits" ]; then
  WITH_FITS=1
fi

mkdir -p fits/glmer_ladder

# always pull the small per-fit summary CSVs + per-kid BLUP CSVs
# (the ranef_*.csv files are what you need for demographic regressions
# on alpha_i / zeta_i across all 7 languages).
#
# Note: sync_from_remote.sh expands wildcards server-side and rsync's
# default flattens them into the local destination root. So pulled
# files land in fits/ regardless of their remote subdirectory. We
# move them into fits/glmer_ladder/ after the pull.
bash sherlock/sync_from_remote.sh 'glmer_ladder/summary_*.csv'
bash sherlock/sync_from_remote.sh 'glmer_ladder/ranef_*.csv'

if [ "$WITH_FITS" = "1" ]; then
  echo
  echo "Pulling fit RDS files (this may take a while) ..."
  bash sherlock/sync_from_remote.sh 'glmer_ladder/fit_*.rds'
fi

# Re-home the pulled files into fits/glmer_ladder/.
mv fits/summary_*.csv  fits/glmer_ladder/ 2>/dev/null || true
mv fits/ranef_*.csv    fits/glmer_ladder/ 2>/dev/null || true
mv fits/fit_*.rds      fits/glmer_ladder/ 2>/dev/null || true

echo
echo "After sync, run:"
echo "  Rscript model/scripts/glmer_ladder/03_aggregate.R"
