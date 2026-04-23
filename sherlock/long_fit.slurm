#!/usr/bin/env bash
#SBATCH --job-name=long_fit
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=normal

# Usage:
#   sbatch sherlock/long_fit.slurm <variant>
# Examples:
#   sbatch sherlock/long_fit.slurm long_2pl_slopes           # English
#   sbatch sherlock/long_fit.slurm long_2pl_slopes_nor       # Norwegian
#   sbatch sherlock/long_fit.slurm long_baseline             # Rasch, no slopes

set -euo pipefail

VARIANT="${1:-long_2pl_slopes}"

# Move into the project root (assumes $HOME/standard_model_2)
cd "$HOME/standard_model_2"

# Results directory on $SCRATCH (persistent enough; 90-day purge)
RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"

# Point rstan's auto_write cache at $SCRATCH so we don't rebuild each time
export RSTAN_AUTO_WRITE_CACHE_DIR="$RESULTS_DIR/stan_cache"
mkdir -p "$RSTAN_AUTO_WRITE_CACHE_DIR"

# Load R
ml R

echo "===== SLURM job $SLURM_JOB_ID ====="
echo "Host: $(hostname)"
echo "Variant: $VARIANT"
echo "Start: $(date)"
echo

# Override PATHS$fits_dir for this job so outputs go to $SCRATCH
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

Rscript model/scripts/fit_longitudinal.R "$VARIANT"

echo
echo "End: $(date)"
