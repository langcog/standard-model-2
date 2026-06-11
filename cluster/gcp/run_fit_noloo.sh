#!/bin/bash
# LOO-free single-fit driver for the GCP VMs. Identical to run_fit.sh but:
#   * always skips save_object (slim summaries recovered from CSVs), and
#   * NEVER runs extract_loo_thinned.R (as_cmdstan_fit OOMs on big CSVs;
#     LOO isn't needed for the sigma_r sensitivity / variance-decomp runs).
# Use for validation pins and any run where you only need scalar posteriors.
#
# Usage:  gcp/run_fit_noloo.sh <variant> <dataset>
set -euo pipefail

VARIANT="${1:?Usage: $0 <variant> <dataset>}"
DATASET="${2:?Usage: $0 <variant> <dataset>}"

cd "$HOME/standard_model_2"
export STANDARD_MODEL_ROOT="$PWD"
export STAN_ITER="${STAN_ITER:-2000}"
export STAN_WARMUP="${STAN_WARMUP:-1000}"
export STAN_ADAPT_DELTA="${STAN_ADAPT_DELTA:-0.95}"
export STAN_MAX_TREEDEPTH="${STAN_MAX_TREEDEPTH:-10}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-16}"
export STANDARD_MODEL_FITS_DIR="$PWD/fits"
export SCRATCH="$PWD/_local"
mkdir -p fits/summaries "$SCRATCH/standard_model_2"
SUMM_LINK="$SCRATCH/standard_model_2/summaries"
if [ -e "$SUMM_LINK" ] && [ ! -L "$SUMM_LINK" ]; then
  mv -n "$SUMM_LINK"/* fits/summaries/ 2>/dev/null || true
  rmdir "$SUMM_LINK" 2>/dev/null || rm -rf "$SUMM_LINK"
fi
ln -sfn "$PWD/fits/summaries" "$SUMM_LINK"

TAG="${VARIANT}"
[ "$DATASET" != "english" ] && TAG="${VARIANT}_${DATASET}"
# reflect the sigma_r override in the tag the extract step looks for
if [ -n "${STAN_SIGMA_R_OVERRIDE:-}" ]; then
  TAG="${TAG}_sigmaR_$(printf '%.2f' "$STAN_SIGMA_R_OVERRIDE" | sed 's/\./p/')"
fi

LOG="$PWD/gcp_${TAG}.log"
echo "==== Starting LOO-free fit: tag=$TAG ====" | tee -a "$LOG"
echo "iter=$STAN_ITER warmup=$STAN_WARMUP threads=$STAN_THREADS_PER_CHAIN sigma_r=${STAN_SIGMA_R_OVERRIDE:-bundle}" | tee -a "$LOG"
echo "start: $(date)" | tee -a "$LOG"

STAN_SKIP_SAVE_OBJECT=1 time Rscript model/scripts/fit_longitudinal.R "$VARIANT" "$DATASET" 2>&1 | tee -a "$LOG"

echo "==== Recovering slim summaries (no LOO): $TAG ====" | tee -a "$LOG"
Rscript sherlock/recover_from_csvs.R "$TAG" 2>&1 | tee -a "$LOG"

SUMM_FILE="$STANDARD_MODEL_FITS_DIR/summaries/${TAG}.summary.rds"
CSV_DIR="$STANDARD_MODEL_FITS_DIR/csvs_$TAG"
if [ -f "$SUMM_FILE" ] && [ -d "$CSV_DIR" ]; then
  echo "Reclaiming $(du -sh "$CSV_DIR" 2>/dev/null | cut -f1) from $CSV_DIR (summary present)" | tee -a "$LOG"
  rm -rf "$CSV_DIR"
elif [ -d "$CSV_DIR" ]; then
  echo "NOT reclaiming $CSV_DIR: $SUMM_FILE missing -- CSVs preserved" | tee -a "$LOG"
fi
echo "==== Done: $TAG ($(date)) ====" | tee -a "$LOG"
