#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=normal
# Forward (prospective) cross-validation: fit M2/M3 on each child's EARLIER
# administrations, score their held-out LAST one. See 01c_fit_forward_cv.R for why
# item-level LOO cannot do this test.
#
# Only ONE refit per (dataset, model) is needed -- every child's last administration
# is held out simultaneously -- so the whole campaign is 10 jobs, not K-fold.
#
# Prereq (run once, locally or on the cluster):
#   Rscript studies/bayes_long/00c_prepare_forward_cv.R _a3
#
# Usage: sbatch --job-name=fcv_<slug>_<model> studies/bayes_long/fcv.slurm <slug> <model>
# Full campaign (owner node for the long pole, owners for the rest):
#   sbatch -p mcfrank --cpus-per-task=24 --export=ALL,STAN_THREADS=6 \
#          --job-name=fcv_norwegian_m3 studies/bayes_long/fcv.slurm norwegian m3
#   for s in japanese smith marchman thal; do for m in m2 m3; do
#     sbatch -p owners --requeue --job-name=fcv_${s}_${m} \
#            studies/bayes_long/fcv.slurm $s $m; done; done
#
# Partition notes are identical to fit.slurm: -p mcfrank = owned node sh04-17n29
# (24 cores -> --cpus-per-task=24 + STAN_THREADS=6); -p owners = preemptible, add
# --requeue; -p normal = shared fairshare pool.
set -euo pipefail
SLUG="$1"; MODEL="$2"
cd "$HOME/standard_model_2"
ml R
export STAN_CHAINS="${STAN_CHAINS:-4}" STAN_WARMUP="${STAN_WARMUP:-1000}" STAN_ITER="${STAN_ITER:-1000}" STAN_ADAPT_DELTA="${STAN_ADAPT_DELTA:-0.9}"
export STAN_THREADS="${STAN_THREADS:-8}"
export CV_DRAWS="${CV_DRAWS:-400}"
echo "===== job $SLURM_JOB_ID  fcv $SLUG/$MODEL  host $(hostname)  start $(date) ====="
Rscript studies/bayes_long/01c_fit_forward_cv.R "$SLUG" "$MODEL"
echo "===== end $(date) ====="
