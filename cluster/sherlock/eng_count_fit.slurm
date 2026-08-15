#!/usr/bin/env bash
#SBATCH --job-name=eng_count
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=normal

# English-only + sumscore count branch stability checks. Usage:
#   sbatch sherlock/eng_count_fit.slurm 2 proc      (WITH processing: bi-lean model, rung D'2)
#   sbatch sherlock/eng_count_fit.slurm 0 noproc    (NO processing: io+count model)
#   sbatch sherlock/eng_count_fit.slurm 2 proconly  (NO input: proc+count model, mirror of noproc)
# Needs fits/joint_io_proc_english_count_subset_data.rds (scp'd; gitignored).
set -euo pipefail
RUNG="${1:-2}"; VARIANT="${2:-proc}"
cd "$HOME/standard_model_2"

RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits/summaries" "$RESULTS_DIR/logs"
for f in "$HOME/standard_model_2/fits"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"; dst="$RESULTS_DIR/fits/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

if [ "$VARIANT" = "noproc" ]; then
  export STAN_MODEL=model/stan/log_irt_long_io_count.stan
elif [ "$VARIANT" = "proconly" ]; then
  export STAN_MODEL=model/stan/log_irt_long_proc_count.stan
else
  export STAN_MODEL=model/stan/log_irt_long_proc_bilingual.stan
fi
export STAN_MM_BUNDLE=joint_io_proc_english_count_subset_data.rds
export STAN_CHAINS=4 STAN_THREADS_PER_CHAIN=8
export STAN_WARMUP=1000 STAN_ITER=1000 STAN_ADAPT_DELTA=0.95 STAN_MAX_TREEDEPTH=10 STAN_DELTA_SD=10

ml R
echo "===== SLURM $SLURM_JOB_ID  ENG+COUNT $VARIANT  rung D'$RUNG  host $(hostname)  $(date) ====="
Rscript model/scripts/fit_joint_io_proc_lean.R "$RUNG"
echo "End: $(date)"
