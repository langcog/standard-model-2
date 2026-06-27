#!/usr/bin/env bash
#SBATCH --job-name=verify_simplify
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=normal
# Phase-0 verification: fit the CURRENT lean model vs the RASCH-simplified model on the
# 50x50 subsample bundle; the simplification is OK iff key scalars match.
set -euo pipefail
cd "$HOME/standard_model_2"
RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits/summaries" "$HOME/standard_model_2/fits/summaries"
for f in "$HOME/standard_model_2/fits"/*.rds; do
  [ -e "$f" ] || continue
  ln -sfn "$f" "$RESULTS_DIR/fits/$(basename "$f")" 2>/dev/null || true
done
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STAN_MM_BUNDLE=joint_io_proc_sub50x50.rds
export STAN_CHAINS=4 STAN_THREADS_PER_CHAIN=4
export STAN_WARMUP=500 STAN_ITER=500 STAN_ADAPT_DELTA=0.95 STAN_DELTA_SD=10
ml R
echo "===== LEAN baseline  $(date) ====="
Rscript model/scripts/fit_joint_io_proc_lean.R 2
echo "===== RASCH simplified  $(date) ====="
STAN_MODEL=model/stan/log_irt_long_proc_rasch.stan Rscript model/scripts/fit_joint_io_proc_lean.R 2
# stage summaries back to $HOME for scp
cp "$RESULTS_DIR/fits/summaries/"joint_io_proc_lean*_sub50x50.summary.rds "$HOME/standard_model_2/fits/summaries/" 2>/dev/null || true
echo "End: $(date)"
