#!/bin/bash
# Driver to run the 5-variant family build for one language sequentially.
# Each fit warm-starts from the prior fit (smart init), saving ~10-15%
# warmup time per fit.
#
# Usage:
#   gcp/run_family.sh [english|norwegian]
# Default: english.

set -euo pipefail

DATASET="${1:-english}"

cd "$HOME/standard_model_2"

# Order: small -> large. Each fit inits from the previous one's summary
# (param subset; missing params are drawn from the prior in cmdstanr).
VARIANTS=(
  "long_demo_pure"
  "long_demo_alpha"
  "long_demo_kappa_pop"
  "long_no_freq_slopes"
  "long_no_freq_slopes_si_signed"
)

PREV_TAG=""
for V in "${VARIANTS[@]}"; do
  TAG="$V"
  [ "$DATASET" != "english" ] && TAG="${V}_${DATASET}"
  echo ""
  echo "============================================================"
  echo "Family build: $TAG  (init from: ${PREV_TAG:-<none>})"
  echo "============================================================"

  bash gcp/run_fit.sh "$V" "$DATASET" "$PREV_TAG"

  PREV_TAG="$TAG"
done

echo ""
echo "============================================================"
echo "Family build complete: $DATASET"
echo "============================================================"
