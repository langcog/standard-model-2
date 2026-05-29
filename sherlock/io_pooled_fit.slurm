#!/usr/bin/env bash
#SBATCH --job-name=io_pooled
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=normal

# Fit the pooled hierarchical IO model over all four datasets
# (BabyView + SEEDLingS + AM2018 + FMW2013).
#
# The CDI likelihood (~4e5 Bernoulli obs) is threaded via reduce_sum:
# 4 chains x 4 threads/chain = 16 cpus (matches --cpus-per-task).
#
# Usage:
#   sbatch sherlock/io_pooled_fit.slurm
#
# Stan iter overridable:
#   sbatch --export=ALL,STAN_ITER=1500,STAN_WARMUP=750 sherlock/io_pooled_fit.slurm

set -euo pipefail

cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"

# Symlink repo-side input bundles into scratch fits dir (race-tolerant)
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
echo "Job: pooled hierarchical IO"
echo "Start: $(date)"
echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STAN_CHAINS="${STAN_CHAINS:-4}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-4}"

Rscript model/scripts/fit_io_pooled.R

echo
echo "End: $(date)"
