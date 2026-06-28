#!/usr/bin/env bash
#SBATCH --job-name=eng_count_refit
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=normal

# Convergence refit of the EN+count +proc fit: 2000/2000 + adapt_delta 0.97 to mix the weakly
# identified gamma_in. Distinct tag (_2k) so it does NOT clobber the current enct summary.
#   sbatch sherlock/eng_count_refit.slurm
set -euo pipefail
cd "$HOME/standard_model_2"
RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits/summaries" "$RESULTS_DIR/logs"
for f in "$HOME/standard_model_2/fits"/*.rds; do
  [ -e "$f" ] || continue; base="$(basename "$f")"; dst="$RESULTS_DIR/fits/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"
export STAN_MODEL=model/stan/log_irt_long_proc_bilingual.stan
export STAN_MM_BUNDLE=joint_io_proc_english_count_subset_data.rds
export STAN_TAG_SFX=_2k
export STAN_CHAINS=4 STAN_THREADS_PER_CHAIN=8
export STAN_WARMUP=2000 STAN_ITER=2000 STAN_ADAPT_DELTA=0.97 STAN_MAX_TREEDEPTH=12 STAN_DELTA_SD=10
ml R
echo "===== SLURM $SLURM_JOB_ID  EN+count +proc REFIT 2000/2000 adapt0.97  $(date) ====="
Rscript model/scripts/fit_joint_io_proc_lean.R 2
echo "End: $(date)"
