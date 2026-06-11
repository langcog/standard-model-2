#!/usr/bin/env bash
#SBATCH --job-name=joint_io_proc
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=normal

# Joint input+processing model (log_irt_long_proc_dp_joint.stan) on all 5
# datasets. Usage: sbatch sherlock/joint_io_proc_fit.slurm <rung>   (rung 0-3)
# Whole ladder: for r in 0 1 2 3; do sbatch sherlock/joint_io_proc_fit.slurm $r; done
set -euo pipefail
RUNG="${1:-1}"
cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits" "$RESULTS_DIR/logs"
for f in "$HOME/standard_model_2/fits"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"; dst="$RESULTS_DIR/fits/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done
export RSTAN_AUTO_WRITE_CACHE_DIR="$RESULTS_DIR/stan_cache"; mkdir -p "$RSTAN_AUTO_WRITE_CACHE_DIR"
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

ml R
echo "===== SLURM $SLURM_JOB_ID  rung D'$RUNG  host $(hostname)  start $(date) ====="
Rscript model/scripts/fit_joint_io_proc.R "$RUNG"
echo "End: $(date)"
