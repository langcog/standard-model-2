#!/bin/bash
# Queue: wait for the NO M_best fit (long_no_freq_slopes norwegian) to
# finish, then refit all IO + Peekbank/proc variants with the new
# DEFAULT_PRIORS (s pinned at 0, no s_i). Each refit also runs the
# scalar / summary / delta_j extracts. LOO extract is skipped -- it's
# the slowest step and Mike isn't using LOO for these big models.
#
# Usage: nohup bash gcp/queue_io_proc.sh > queue_io_proc.log 2>&1 &

set -uo pipefail

cd "$HOME/standard_model_2"
export STANDARD_MODEL_ROOT="$PWD"
export STAN_ITER="${STAN_ITER:-2000}"
export STAN_WARMUP="${STAN_WARMUP:-1000}"
export STAN_ADAPT_DELTA="${STAN_ADAPT_DELTA:-0.95}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-16}"
export STAN_MAX_TREEDEPTH="${STAN_MAX_TREEDEPTH:-10}"
export STANDARD_MODEL_FITS_DIR="$PWD/fits"
export SCRATCH="$PWD/_local"
mkdir -p fits/summaries "$SCRATCH/standard_model_2"
SUMM_LINK="$SCRATCH/standard_model_2/summaries"
[ -e "$SUMM_LINK" ] && [ ! -L "$SUMM_LINK" ] && rm -rf "$SUMM_LINK"
ln -sfn "$PWD/fits/summaries" "$SUMM_LINK"

LOG="$PWD/queue_io_proc.log"
echo "==== queue_io_proc started $(date -u) ====" | tee -a "$LOG"

# Poll for the NO run_fit.sh to exit (covers sampling + extracts).
echo "Waiting for NO M_best pipeline to finish..." | tee -a "$LOG"
while pgrep -af 'run_fit.sh long_no_freq_slopes norwegian' \
        | grep -v 'queue_io_proc\|grep' > /dev/null; do
  sleep 60
done
echo "NO pipeline finished at $(date -u). Starting IO/proc refits." | tee -a "$LOG"

# ---- One-fit runner --------------------------------------------------
# Args: script  variant  dataset_or_blank  tag
run_one() {
  local script=$1
  local variant=$2
  local dataset=$3   # empty string means "no dataset arg" (stanford)
  local tag=$4

  echo "==== $(date -u): Refit $tag (script=$script variant=$variant dataset=$dataset) ====" | tee -a "$LOG"

  # Force a fresh fit by deleting the cached .rds first. fit_variant_cmdstanr
  # short-circuits with "already fit, loading" if the file exists.
  rm -f "fits/${tag}.rds"

  if [ -n "$dataset" ]; then
    time Rscript "model/scripts/${script}" "$variant" "$dataset" 2>&1 | tee -a "$LOG"
  else
    time Rscript "model/scripts/${script}" "$variant" 2>&1 | tee -a "$LOG"
  fi

  echo "==== $(date -u): Extracts for $tag ====" | tee -a "$LOG"
  Rscript sherlock/extract_summary_table_only.R "$tag" 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_scalar_draws.R     "$tag" 2>&1 | tee -a "$LOG"
  Rscript sherlock/extract_delta_j_slim.R     "$tag" 2>&1 | tee -a "$LOG"
  # Skipping extract_loo_thinned.R (slowest extract; not used in current
  # headline analyses per Mike, 2026-05-22).
  echo "==== $(date -u): Done $tag ====" | tee -a "$LOG"
}

# ---- IO / BabyView ---------------------------------------------------
run_one fit_io.R io_no_freq_slopes babyview io_no_freq_slopes

# ---- IO / SEEDLingS --------------------------------------------------
run_one fit_io.R io_no_freq_slopes               seedlings io_no_freq_slopes_seedlings
run_one fit_io.R io_comp_no_freq_slopes          seedlings io_comp_no_freq_slopes_seedlings
run_one fit_io.R io_std_no_freq_slopes           seedlings io_std_no_freq_slopes_seedlings
run_one fit_io.R io_comp_std_no_freq_slopes      seedlings io_comp_std_no_freq_slopes_seedlings

# ---- Peekbank / Stanford linked --------------------------------------
run_one fit_stanford_linked.R long_proc_no_freq_slopes "" long_proc_no_freq_slopes

echo "==== queue_io_proc all done $(date -u) ====" | tee -a "$LOG"
