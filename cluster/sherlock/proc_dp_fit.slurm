#!/usr/bin/env bash
#SBATCH --job-name=proc_dp
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --partition=normal

# D'0-D'3 processing regression ladder (log_irt_long_proc_dp.stan).
#
# Usage:
#   sbatch sherlock/proc_dp_fit.slurm <rung> [dataset_tag]
#     rung: 0,1,2,3   dataset_tag: all (default) | adams
# Run the whole ladder:
#   for r in 0 1 2 3; do sbatch sherlock/proc_dp_fit.slurm $r all; done
#
# Stan overridable: --export=ALL,STAN_ITER=1000,STAN_WARMUP=750,STAN_CHAINS=4

set -euo pipefail
RUNG="${1:-0}"
DSTAG="${2:-all}"

cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"

# Symlink repo-side input bundles into scratch fits dir (race-tolerant)
for f in "$HOME/standard_model_2/fits"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"; dst="$RESULTS_DIR/fits/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done

export RSTAN_AUTO_WRITE_CACHE_DIR="$RESULTS_DIR/stan_cache"
mkdir -p "$RSTAN_AUTO_WRITE_CACHE_DIR"
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

ml R

echo "===== SLURM job $SLURM_JOB_ID ====="
echo "Host: $(hostname)"
echo "Rung: D'$RUNG  Dataset: $DSTAG"
echo "Start: $(date)"
echo

Rscript model/scripts/fit_proc_dp.R "$RUNG" "$DSTAG"

echo
echo "End: $(date)"
