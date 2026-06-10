#!/usr/bin/env bash
#SBATCH --job-name=io_pooled_widedelta
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=normal

# Diagnostic refit of pooled IO with delta prior widened from N(0,0.5)
# to N(0,10) -- effectively uninformative over the plausible kappa range.
# Tests whether the default tight delta prior was the deflator.

set -euo pipefail
cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"
for f in "$HOME/standard_model_2/fits"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"
  dst="$RESULTS_DIR/fits/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done

export RSTAN_AUTO_WRITE_CACHE_DIR="$RESULTS_DIR/stan_cache"
mkdir -p "$RSTAN_AUTO_WRITE_CACHE_DIR"
ml R

echo "===== SLURM job $SLURM_JOB_ID ====="
echo "Host: $(hostname)"
echo "Job: pooled IO with widened delta prior"
echo "Start: $(date)"
echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STAN_CHAINS="${STAN_CHAINS:-4}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-4}"

Rscript model/scripts/fit_io_pooled_widedelta.R

echo
echo "End: $(date)"
