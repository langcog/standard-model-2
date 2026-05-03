#!/usr/bin/env bash
#SBATCH --job-name=smoke
#SBATCH --output=smoke_%j.out
#SBATCH --error=smoke_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=normal

# Smoke test for the full Sherlock pipeline: recovery sim on synthetic
# data. Writes fits/recovery.rds to $SCRATCH. Rsync it back and
# inspect to confirm the full cycle works.

set -euo pipefail
cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"
export RSTAN_AUTO_WRITE_CACHE_DIR="$RESULTS_DIR/stan_cache"
mkdir -p "$RSTAN_AUTO_WRITE_CACHE_DIR"

ml R

echo "===== SLURM job $SLURM_JOB_ID (smoke test) ====="
echo "Host: $(hostname)"
echo "Start: $(date)"
echo

export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STANDARD_MODEL_FIGS_DIR="$RESULTS_DIR/figs"

Rscript model/scripts/01_recovery.R

echo
echo "End: $(date)"
echo
echo "Output files:"
ls -la "$RESULTS_DIR/fits/"
