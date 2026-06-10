#!/bin/bash
# Extension fits: BabyView, SEEDLingS, Peekbank/Stanford-linked.
# Runs two waves:
#   Wave A (no si_signed): refits the existing variants at s=6 from
#     config.R DEFAULT_PRIORS. Gets us refreshed sigma_alpha, pi_alpha,
#     kappa_pop, sigma_zeta for the §19 cross-sample table.
#   Wave B (si_signed):    refits with per-child trajectory phase added
#     via the new log_irt_long_*_si_signed.stan files. Gets us sigma_s
#     for the same samples, completing the apples-to-apples §19 table.
#
# Each fit auto-skips if its .rds is already on disk (see skip-if-done
# logic). So a partial completion can be resumed by re-running.
#
# Usage:
#   gcp/run_extensions.sh

# Note: deliberately NOT using `set -e`. We want individual fit failures
# (e.g., missing bundle, runtime Stan error) to be logged and skipped
# rather than killing the entire queue. The skip-if-done check at the
# top of run_one means a re-launch will pick up where we left off.
set -uo pipefail

cd "$HOME/standard_model_2"

# Standard env (matches run_fit.sh; defaults are GCP-VM appropriate)
export STAN_ITER="${STAN_ITER:-2000}"
export STAN_WARMUP="${STAN_WARMUP:-1000}"
export STAN_ADAPT_DELTA="${STAN_ADAPT_DELTA:-0.95}"
export STAN_MAX_TREEDEPTH="${STAN_MAX_TREEDEPTH:-10}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-16}"
export STAN_BACKEND="${STAN_BACKEND:-cmdstanr}"
export STANDARD_MODEL_ROOT="$PWD"
export STANDARD_MODEL_FITS_DIR="$PWD/fits"
export SCRATCH="$PWD/_local"

mkdir -p fits/summaries
mkdir -p "$SCRATCH/standard_model_2"
SUMM_LINK="$SCRATCH/standard_model_2/summaries"
if [ -e "$SUMM_LINK" ] && [ ! -L "$SUMM_LINK" ]; then
  mv -n "$SUMM_LINK"/* fits/summaries/ 2>/dev/null || true
  rmdir "$SUMM_LINK" 2>/dev/null || rm -rf "$SUMM_LINK"
fi
ln -sfn "$PWD/fits/summaries" "$SUMM_LINK"

# Fit specs: (variant, dataset, fitter_script, tag).
# Note: fit_io.R names babyview fits with bare variant; non-babyview fits
# get a "_<dataset>" suffix appended. fit_stanford_linked.R uses the
# variant as the tag directly.
#
# Wave A: existing variants at s=6 (no si_signed)
WAVE_A=(
  "io_no_freq_slopes|babyview|model/scripts/fit_io.R|io_no_freq_slopes"
  "io_comp_no_freq_slopes|seedlings|model/scripts/fit_io.R|io_comp_no_freq_slopes_seedlings"
  "long_proc_no_freq_slopes|stanford_linked|model/scripts/fit_stanford_linked.R|long_proc_no_freq_slopes"
)

# Wave B: si_signed extensions
WAVE_B=(
  "io_no_freq_slopes_si_signed|babyview|model/scripts/fit_io.R|io_no_freq_slopes_si_signed"
  "io_comp_no_freq_slopes_si_signed|seedlings|model/scripts/fit_io.R|io_comp_no_freq_slopes_si_signed_seedlings"
  "long_proc_no_freq_slopes_si_signed|stanford_linked|model/scripts/fit_stanford_linked.R|long_proc_no_freq_slopes_si_signed"
)

run_one() {
  local variant="$1"
  local dataset="$2"
  local script="$3"
  local tag="$4"

  echo ""
  echo "============================================================"
  echo "Extension fit: $tag"
  echo "  variant: $variant"
  echo "  dataset: $dataset"
  echo "  script:  $script"
  echo "============================================================"

  if [ -f "fits/${tag}.rds" ]; then
    echo "  Already on disk; SKIPPING."
    return 0
  fi

  local log="$PWD/gcp_${tag}.log"
  echo "start: $(date)" | tee -a "$log"

  # fit_stanford_linked.R takes only one positional arg (variant).
  if [[ "$script" == *fit_stanford_linked.R ]]; then
    time Rscript "$script" "$variant" 2>&1 | tee -a "$log"
  else
    time Rscript "$script" "$variant" "$dataset" 2>&1 | tee -a "$log"
  fi

  echo "fit done: $(date)" | tee -a "$log"

  # Extracts (same set as run_fit.sh; small N here so all extracts cheap)
  echo "==== Extracting: $tag ====" | tee -a "$log"
  Rscript sherlock/extract_summary_table_only.R "$tag" 2>&1 | tee -a "$log"
  Rscript sherlock/extract_scalar_draws.R "$tag"        2>&1 | tee -a "$log"
  # delta_j extract: many extension samples are too small for stable per-item
  # estimates, but the script handles that gracefully (skip empty draws).
  Rscript sherlock/extract_delta_j_slim.R "$tag"        2>&1 | tee -a "$log" || true
  # LOO: small-N fits don't generate big log_lik matrices, but use the
  # chunked extractor anyway (it auto-falls back to single-pass at small N).
  Rscript sherlock/extract_loo_thinned.R "$tag"         2>&1 | tee -a "$log" || true

  echo "end: $(date)" | tee -a "$log"
}

echo ""
echo "############################################################"
echo "# Extensions wave A: no si_signed (s=6 refits)"
echo "############################################################"
for spec in "${WAVE_A[@]}"; do
  IFS='|' read -r variant dataset script tag <<< "$spec"
  run_one "$variant" "$dataset" "$script" "$tag" || \
    echo "[$(date)] FIT FAILED: $tag (continuing to next)" >&2
done

echo ""
echo "############################################################"
echo "# Extensions wave B: + si_signed"
echo "############################################################"
for spec in "${WAVE_B[@]}"; do
  IFS='|' read -r variant dataset script tag <<< "$spec"
  run_one "$variant" "$dataset" "$script" "$tag" || \
    echo "[$(date)] FIT FAILED: $tag (continuing to next)" >&2
done

echo ""
echo "============================================================"
echo "Extensions complete."
echo "============================================================"
