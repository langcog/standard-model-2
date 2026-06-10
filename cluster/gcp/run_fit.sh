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

if [ "${STAN_SKIP_SAVE_OBJECT:-0}" = "1" ]; then
  # save_object was skipped (no <tag>.rds): extract the slim summaries
  # straight from the persistent CSVs. recover_from_csvs streams the scalar
  # columns (incl. gamma_in) + delta_j with `cut` -- low RAM at any size.
  Rscript sherlock/recover_from_csvs.R "$TAG" 2>&1 | tee -a "$LOG"
  # extract_loo_thinned still builds a cmdstanr fit object (as_cmdstan_fit),
  # which OOMs above ~150 GB of CSV (Norwegian). Skip LOO for those -- the
  # variance decomposition doesn't need it.
  CSV_GB=$(du -s "$STANDARD_MODEL_FITS_DIR/csvs_$TAG" 2>/dev/null | awk '{print int($1/1048576)}')
  if [ "${CSV_GB:-999}" -lt 150 ]; then
    Rscript sherlock/extract_loo_thinned.R "$TAG" 2>&1 | tee -a "$LOG"
  else
    echo "Skipping LOO: csvs_$TAG is ${CSV_GB} GB (> 150 GB threshold, would OOM)" | tee -a "$LOG"
  fi
else
  Rscript sherlock/extract_summary_table_only.R "$TAG" 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_scalar_draws.R "$TAG" 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_delta_j_slim.R "$TAG" 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_loo_thinned.R "$TAG" 2>&1 | tee -a "$LOG"
fi

# Reclaim the raw sampler CSVs (tens to ~210 GB per fit) -- but ONLY if the
# slim summary was actually written. A failed/empty extraction must keep
# the CSVs as the sole recovery copy; deleting them regardless is how EN D'
# was lost. Without reclaiming at all, a D -> D' chain accumulates two CSV
# sets and fills the boot disk mid-fit.
SUMM_FILE="$STANDARD_MODEL_FITS_DIR/summaries/${TAG}.summary.rds"
CSV_DIR="$STANDARD_MODEL_FITS_DIR/csvs_$TAG"
if [ -f "$SUMM_FILE" ] && [ -d "$CSV_DIR" ]; then
  echo "Reclaiming $(du -sh "$CSV_DIR" 2>/dev/null | cut -f1) from $CSV_DIR (summary present)" | tee -a "$LOG"
  rm -rf "$CSV_DIR"
elif [ -d "$CSV_DIR" ]; then
  echo "NOT reclaiming $CSV_DIR: $SUMM_FILE missing -- CSVs preserved for recovery" | tee -a "$LOG"
fi

echo "==== Done: $TAG ====" | tee -a "$LOG"
echo "end: $(date)" | tee -a "$LOG"
