#!/usr/bin/env bash
#SBATCH --job-name=joint_io_proc_bi
#SBATCH --output=%x_d%a_%j.out
#SBATCH --error=%x_d%a_%j.err
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=normal

# BILINGUAL ("bi-lean") io+proc fit: English item-level + Spanish item-level (SLENA)
# + the sumscore count branch (ELENA-WS, HABLA). Model log_irt_long_proc_bilingual.stan,
# bundle joint_io_proc_bilingual_subset_data.rds (scp'd; gitignored). Usage:
#   sbatch sherlock/joint_io_proc_bi_fit.slurm 2            (full model, rung D'2)
#   sbatch --array=0-2 sherlock/joint_io_proc_bi_fit.slurm  (ladder D'0-D'2)
set -euo pipefail
RUNG="${SLURM_ARRAY_TASK_ID:-${1:-2}}"
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

export STAN_MM_BUNDLE=joint_io_proc_bilingual_subset_data.rds
export STAN_MODEL=model/stan/log_irt_long_proc_bilingual.stan
export STAN_CHAINS=4 STAN_THREADS_PER_CHAIN=8
export STAN_WARMUP=1000 STAN_ITER=1000 STAN_ADAPT_DELTA=0.95 STAN_MAX_TREEDEPTH=10 STAN_DELTA_SD=10

ml R
echo "===== SLURM $SLURM_JOB_ID  BILINGUAL rung D'$RUNG  host $(hostname)  $(date) ====="
Rscript model/scripts/fit_joint_io_proc_lean.R "$RUNG"
echo "End: $(date)"
