#!/usr/bin/env bash
#SBATCH --job-name=io_pooled_unanchored
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=normal

# Diagnostic refit of the pooled IO model with delta_j anchor RELAXED
# to SD = 5.0 for all items (instead of the default 0.10 for items
# matched to the EN longitudinal anchor). Tests whether the anchor is
# what's deflating per-study kappa.
#
# Usage:
#   sbatch sherlock/io_pooled_unanchored_fit.slurm

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
echo "Job: pooled IO (delta_j anchor relaxed)"
echo "Start: $(date)"
echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STAN_CHAINS="${STAN_CHAINS:-4}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-4}"

Rscript model/scripts/fit_io_pooled_unanchored.R

echo
echo "End: $(date)"
