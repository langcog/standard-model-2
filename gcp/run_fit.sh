#!/bin/bash
# Single-fit driver for the GCP VM.
# Runs one variant on one dataset, extracts summaries + delta_j + LOO,
# saves to fits/ and the summaries dir. Intended to be called in a
# loop for the 5-variant family build.
#
# Usage:
#   gcp/run_fit.sh <variant> <dataset> [init_from_tag]
#
# Env (passed to fit_longitudinal.R):
#   STAN_ITER          (default 4000)
#   STAN_WARMUP        (default 2500)
#   STAN_ADAPT_DELTA   (default 0.95)
#   STAN_THREADS_PER_CHAIN  (default 8)

set -euo pipefail

VARIANT="${1:?Usage: $0 <variant> <dataset> [init_from_tag]}"
DATASET="${2:?Usage: $0 <variant> <dataset> [init_from_tag]}"
INIT_FROM="${3:-}"

cd "$HOME/standard_model_2"

# Tell config.R where the project root is (defaults to mcfrank's Mac path).
export STANDARD_MODEL_ROOT="$PWD"

export STAN_ITER="${STAN_ITER:-4000}"
export STAN_WARMUP="${STAN_WARMUP:-2500}"
export STAN_ADAPT_DELTA="${STAN_ADAPT_DELTA:-0.95}"
export STAN_MAX_TREEDEPTH="${STAN_MAX_TREEDEPTH:-10}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-8}"
[ -n "$INIT_FROM" ] && export STAN_INIT_FROM="$INIT_FROM"

# Direct extract scripts to write to <repo>/fits/summaries (where the
# rest of the local pipeline reads them). The Sherlock-side scripts
# default to $SCRATCH/standard_model_2/summaries; set SCRATCH so that
# resolves to a location inside our repo.
export STANDARD_MODEL_FITS_DIR="$PWD/fits"
export SCRATCH="$PWD/_local"
mkdir -p fits/summaries
mkdir -p "$SCRATCH/standard_model_2"
# Symlink so extracts written to $SCRATCH/standard_model_2/summaries
# land in fits/summaries/. Previous version raced mkdir -p in front of
# the symlink check, leaving a real directory that swallowed output.
SUMM_LINK="$SCRATCH/standard_model_2/summaries"
if [ -e "$SUMM_LINK" ] && [ ! -L "$SUMM_LINK" ]; then
  echo "Migrating pre-existing $SUMM_LINK directory into fits/summaries/ ..."
  mv -n "$SUMM_LINK"/* fits/summaries/ 2>/dev/null || true
  rmdir "$SUMM_LINK" 2>/dev/null || rm -rf "$SUMM_LINK"
fi
ln -sfn "$PWD/fits/summaries" "$SUMM_LINK"

TAG="${VARIANT}"
[ "$DATASET" != "english" ] && TAG="${VARIANT}_${DATASET}"

LOG="$PWD/gcp_${TAG}.log"
echo "==== Starting fit: tag=$TAG ====" | tee -a "$LOG"
echo "iter=$STAN_ITER warmup=$STAN_WARMUP adapt_delta=$STAN_ADAPT_DELTA threads=$STAN_THREADS_PER_CHAIN" | tee -a "$LOG"
[ -n "$INIT_FROM" ] && echo "init_from=$INIT_FROM" | tee -a "$LOG"
echo "start: $(date)" | tee -a "$LOG"

time Rscript model/scripts/fit_longitudinal.R "$VARIANT" "$DATASET" 2>&1 | tee -a "$LOG"

echo "end fit: $(date)" | tee -a "$LOG"
echo "==== Extracting: $TAG ====" | tee -a "$LOG"

Rscript sherlock/extract_summary_table_only.R "$TAG" 2>&1 | tee -a "$LOG"
Rscript sherlock/extract_scalar_draws.R "$TAG" 2>&1 | tee -a "$LOG"
Rscript sherlock/extract_delta_j_slim.R "$TAG" 2>&1 | tee -a "$LOG"
Rscript sherlock/extract_loo_thinned.R "$TAG" 2>&1 | tee -a "$LOG"

echo "==== Done: $TAG ====" | tee -a "$LOG"
echo "end: $(date)" | tee -a "$LOG"
