#!/usr/bin/env bash
#SBATCH --job-name=io_pooled_gamma_widedelta
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=normal

# Wide-delta version of the gamma variant fits (add / mult).
# Usage:
#   sbatch sherlock/io_pooled_gamma_widedelta_fit.slurm add
#   sbatch sherlock/io_pooled_gamma_widedelta_fit.slurm mult

set -euo pipefail
VARIANT="${1:-add}"
if [[ "$VARIANT" != "add" && "$VARIANT" != "mult" ]]; then
  echo "ERROR: variant must be 'add' or 'mult', got '$VARIANT'" >&2
  exit 2
fi

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
echo "Job: pooled IO gamma-$VARIANT wide-delta"
echo "Start: $(date)"
echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STAN_CHAINS="${STAN_CHAINS:-4}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-4}"

Rscript model/scripts/fit_io_pooled_gamma_widedelta.R "$VARIANT"

echo
echo "End: $(date)"
