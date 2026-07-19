#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
# Fit the mega-model (m_pool) on the combined cross-dataset bundle.
# Owner node:  sbatch -p mcfrank --job-name=blong_pool studies/bayes_long/pool.slurm _a3
set -euo pipefail
SFX="${1:-_a3}"
cd "$HOME/standard_model_2"
ml R
export STAN_CHAINS="${STAN_CHAINS:-4}" STAN_WARMUP="${STAN_WARMUP:-1000}" STAN_ITER="${STAN_ITER:-1000}" STAN_ADAPT_DELTA="${STAN_ADAPT_DELTA:-0.9}"
export STAN_THREADS="${STAN_THREADS:-6}"
echo "===== pool$SFX  job $SLURM_JOB_ID  host $(hostname)  start $(date) ====="
Rscript studies/bayes_long/01b_fit_pool.R "$SFX"
echo "===== end $(date) ====="
