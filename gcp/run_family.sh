#!/bin/bash
# Driver to run the 5-variant family build for one language sequentially.
#
# Warm-start policy (revised after EN demo_alpha treedepth disaster, 2026-05-18):
#   - demo_alpha and demo_kappa_pop: NO warm-start. Initing from demo_pure
#     parks all 4 chains near sigma_alpha = 0 (the boundary), which leads
#     to 100% max-treedepth saturation and Rhat ~ 2.7. Cold-start is safer
#     and the warmup cost is small on these cheaper variants.
#   - slopes and signed_si: keep warm-start. Source variants (kappa_pop,
#     slopes) have freed the relevant parameters, so the init is in a
#     valid region of the posterior, and the warmup savings on these
#     ~14h fits matter.
#
# (Previously also bumped adapt_delta=0.97 and max_treedepth=12 for
# demo_alpha/demo_kappa_pop as belt-and-suspenders. Reverted because
# they caused ~8x per-iter slowdown without being necessary -- cold-start
# alone fixes the boundary problem. Knobs back at defaults.)
#
# Usage:
#   gcp/run_family.sh [english|norwegian]
# Default: english.

set -euo pipefail

DATASET="${1:-english}"

cd "$HOME/standard_model_2"

# Order: small -> large.
VARIANTS=(
  "long_demo_pure"
  "long_demo_alpha"
  "long_demo_kappa_pop"
  "long_no_freq_slopes"
  "long_no_freq_slopes_si_signed"
)

# Variants that get cold-started (no init_from prior fit).
COLD_START_VARIANTS=(
  "long_demo_alpha"
  "long_demo_kappa_pop"
)

contains() {  # args: needle "${haystack[@]}"
  local needle="$1"; shift
  for x in "$@"; do [ "$x" = "$needle" ] && return 0; done
  return 1
}

PREV_TAG=""
for V in "${VARIANTS[@]}"; do
  TAG="$V"
  [ "$DATASET" != "english" ] && TAG="${V}_${DATASET}"

  # Skip if fit already on disk; advance the chain anyway so subsequent
  # warm-starts find the source.
  if [ -f "fits/${TAG}.rds" ]; then
    echo ""
    echo "============================================================"
    echo "Family build: $TAG -- already on disk, SKIPPING"
    echo "============================================================"
    PREV_TAG="$TAG"
    continue
  fi

  # Decide init_from
  if contains "$V" "${COLD_START_VARIANTS[@]}"; then
    INIT_FROM=""
  else
    INIT_FROM="$PREV_TAG"
  fi

  echo ""
  echo "============================================================"
  echo "Family build: $TAG"
  echo "  init from:     ${INIT_FROM:-<cold start>}"
  echo "============================================================"

  bash gcp/run_fit.sh "$V" "$DATASET" "$INIT_FROM"

  PREV_TAG="$TAG"
done

echo ""
echo "============================================================"
echo "Family build complete: $DATASET"
echo "============================================================"
