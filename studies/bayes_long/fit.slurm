#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=normal
# normal QOS caps at 48h (partition physical max is 7d). Be generous with --time
# (a job only uses what it needs; over-requesting costs nothing, a short cap kills
# mid-sample). For >48h fits add: sbatch --qos=long --time=4-00:00:00 ...
#
# Fit one (dataset, model) of the bayes_long ladder on Sherlock.
# Usage: sbatch --job-name=blong_<slug>_<model> studies/bayes_long/fit.slurm <slug> <model>
# Full sweep:
#   for s in japanese thal smith marchman norwegian; do for m in m0 m1 m2 m3; do
#     sbatch --job-name=blong_${s}_${m} studies/bayes_long/fit.slurm $s $m; done; done
set -euo pipefail
SLUG="$1"; MODEL="$2"
cd "$HOME/standard_model_2"
ml R
# 4 parallel chains x 8 threads/chain = 32 cores (reduce_sum within-chain threading)
export STAN_CHAINS=4 STAN_WARMUP=1000 STAN_ITER=1000 STAN_THREADS=8 STAN_ADAPT_DELTA=0.9
echo "===== job $SLURM_JOB_ID  $SLUG/$MODEL  host $(hostname)  start $(date) ====="
Rscript studies/bayes_long/01_fit.R "$SLUG" "$MODEL"
echo "===== end $(date) ====="
