#!/usr/bin/env bash
#SBATCH --job-name=proc_fit
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --partition=normal

# Joint vocab + LWL processing model on the Stanford-linked dataset.
#
# Usage:
#   sbatch sherlock/proc_fit.slurm <variant>
# Examples:
#   sbatch sherlock/proc_fit.slurm long_proc_slopes               # default lean
#   sbatch sherlock/proc_fit.slurm long_proc_2pl_slopes           # adds 2PL
#
# Stan iter overridable: --export=ALL,STAN_ITER=1500,STAN_WARMUP=750

set -euo pipefail
VARIANT="${1:-long_proc_slopes}"

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
echo "Variant: $VARIANT"
echo "Start: $(date)"
echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

Rscript model/scripts/fit_stanford_linked.R "$VARIANT"

echo
echo "End: $(date)"
