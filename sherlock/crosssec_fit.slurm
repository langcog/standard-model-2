#!/usr/bin/env bash
#SBATCH --job-name=cs_fit
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal

# Usage:
#   sbatch sherlock/crosssec_fit.slurm <variant>
# Variants: baseline | fix_delta | fix_s | both_fixed | 2pl | 2pl_fix_delta
#           | slopes | 2pl_slopes

set -euo pipefail
VARIANT="${1:-baseline}"

cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"

# Symlink repo-side input bundles into scratch fits dir
for f in "$HOME/standard_model_2/model/fits"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"
  [ -e "$RESULTS_DIR/fits/$base" ] || ln -s "$f" "$RESULTS_DIR/fits/$base"
done

export RSTAN_AUTO_WRITE_CACHE_DIR="$RESULTS_DIR/stan_cache"
mkdir -p "$RSTAN_AUTO_WRITE_CACHE_DIR"

ml R

echo "===== SLURM job $SLURM_JOB_ID ====="
echo "Host: $(hostname)"; echo "Variant: $VARIANT"; echo "Start: $(date)"; echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

Rscript model/scripts/fit_variant.R "$VARIANT"

echo; echo "End: $(date)"
