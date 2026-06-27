#!/usr/bin/env bash
#SBATCH --job-name=verify_bi
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=normal
set -euo pipefail
cd "$HOME/standard_model_2"
RESULTS_DIR="$SCRATCH/standard_model_2"; mkdir -p "$RESULTS_DIR/fits/summaries" fits/summaries
for f in "$HOME/standard_model_2/fits"/*.rds; do [ -e "$f" ] && ln -sfn "$f" "$RESULTS_DIR/fits/$(basename "$f")" 2>/dev/null || true; done
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STAN_MM_BUNDLE=joint_io_proc_sub50x50_bi.rds STAN_MODEL=model/stan/log_irt_long_proc_bilingual.stan
export STAN_CHAINS=4 STAN_THREADS_PER_CHAIN=4 STAN_WARMUP=500 STAN_ITER=500 STAN_DELTA_SD=10
ml R
Rscript model/scripts/fit_joint_io_proc_lean.R 2
cp "$RESULTS_DIR/fits/summaries/"joint_io_proc_lean_rasch_d2_sub50x50_bi.summary.rds "$HOME/standard_model_2/fits/summaries/" 2>/dev/null || true
