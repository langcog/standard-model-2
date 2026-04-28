#!/usr/bin/env bash
#SBATCH --job-name=io_fit
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --partition=normal

# Usage:
#   sbatch sherlock/io_fit.slurm <variant> [dataset]
# Examples:
#   sbatch sherlock/io_fit.slurm io_2pl_slopes              # BabyView default
#   sbatch sherlock/io_fit.slurm io_baseline babyview
#
# Stan iter overridable:
#   sbatch --export=ALL,STAN_ITER=1500,STAN_WARMUP=750 sherlock/io_fit.slurm ...

set -euo pipefail
VARIANT="${1:-io_2pl_slopes}"
DATASET="${2:-babyview}"

cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"

# Symlink repo-side input bundles into scratch fits dir (race-tolerant)
for f in "$HOME/standard_model_2/model/fits"/*.rds; do
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
echo "Variant: $VARIANT  Dataset: $DATASET"
echo "Start: $(date)"
echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

Rscript model/scripts/fit_io.R "$VARIANT" "$DATASET"

echo
echo "End: $(date)"
