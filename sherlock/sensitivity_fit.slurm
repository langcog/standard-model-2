#!/usr/bin/env bash
#SBATCH --job-name=sens
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --partition=normal

# One sigma_r sensitivity fit at a single value.
# Usage:
#   sbatch sherlock/sensitivity_fit.slurm <variant> <sigma_r>
# Example:
#   sbatch sherlock/sensitivity_fit.slurm 2pl 0.8

set -euo pipefail
VARIANT="${1:-2pl}"
SIGMA_R="${2:-0.53}"

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
echo "Variant: $VARIANT  sigma_r=$SIGMA_R"; echo "Start: $(date)"; echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

Rscript model/scripts/sensitivity_sigma_r.R "$VARIANT" "$SIGMA_R"

echo; echo "End: $(date)"
