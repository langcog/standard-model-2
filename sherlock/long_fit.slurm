#!/usr/bin/env bash
#SBATCH --job-name=long_fit
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=normal
#
# 16 cpus = 4 chains x 4 threads/chain (reduce_sum within-chain
# parallelism). Override via STAN_THREADS_PER_CHAIN env var if needed.

# Usage:
#   sbatch sherlock/long_fit.slurm <variant> [dataset]
# Examples:
#   sbatch sherlock/long_fit.slurm long_slopes                   # English (default; lean + slopes)
#   sbatch sherlock/long_fit.slurm long_slopes norwegian
#   sbatch sherlock/long_fit.slurm long_2pl_slopes english       # adds 2PL on top
#   sbatch sherlock/long_fit.slurm long_baseline english         # lean Rasch, no slopes
# Dataset keys are defined in model/R/datasets.R.
#
# Stan iteration count is overridable via env vars:
#   sbatch --export=STAN_ITER=1000,STAN_WARMUP=500 sherlock/long_fit.slurm ...

set -euo pipefail

VARIANT="${1:-long_slopes}"
DATASET="${2:-english}"

# Move into the project root (assumes $HOME/standard_model_2)
cd "$HOME/standard_model_2"

# Results directory on $SCRATCH (persistent enough; 90-day purge)
RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"

# Make repo input bundles visible from $SCRATCH/.../fits (where Stan reads).
# Symlink any .rds in the repo's model/fits that isn't already in scratch.
for f in "$HOME/standard_model_2/model/fits"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"
  dst="$RESULTS_DIR/fits/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done

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

# cmdstanr threading: each chain gets cpus / chains threads.
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-$((SLURM_CPUS_PER_TASK / 4))}"
echo "Threads per chain: $STAN_THREADS_PER_CHAIN"

Rscript model/scripts/fit_longitudinal.R "$VARIANT" "$DATASET"

echo
echo "End: $(date)"
