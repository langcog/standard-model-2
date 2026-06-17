#!/usr/bin/env bash
#SBATCH --job-name=joint_io_proc_lean
#SBATCH --output=%x_d%a_%j.out
#SBATCH --error=%x_d%a_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --partition=normal

# Joint io+proc MEASUREMENT-MODEL ladder (log_irt_long_proc_dp_joint_mm.stan,
# reduce_sum). Usage:
#   sbatch sherlock/joint_io_proc_lean_fit.slurm <rung>          (single rung 0-3)
#   sbatch --array=0-3 sherlock/joint_io_proc_lean_fit.slurm     (whole ladder)
# Needs fits/joint_io_proc_lean_subset_data.rds present (scp'd; gitignored).
set -euo pipefail
RUNG="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"
for f in "$HOME/standard_model_2/fits"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"; dst="$RESULTS_DIR/fits/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

# 4 chains x 8 threads = 32 cores (reduce_sum within-chain threading)
export STAN_CHAINS=4 STAN_THREADS_PER_CHAIN=8
export STAN_WARMUP=1000 STAN_ITER=1000 STAN_ADAPT_DELTA=0.95 STAN_MAX_TREEDEPTH=10

ml R
echo "===== SLURM $SLURM_JOB_ID  rung D'$RUNG  host $(hostname)  $(date) ====="
Rscript model/scripts/fit_joint_io_proc_lean.R "$RUNG"
echo "End: $(date)"
