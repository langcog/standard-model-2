#!/usr/bin/env bash
#SBATCH --job-name=ladder
#SBATCH --output=logs/ladder_%A_%a.out
#SBATCH --error=logs/ladder_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=normal
# 8 CPUs per task: glmer itself is single-threaded, but RcppEigen / BLAS
# inside the IRLS Cholesky updates auto-parallelizes via OpenMP. Observed
# ~4.5x speedup on NO D_log without SLURM requesting more cores; this
# bumps the request to match and adds a few more for the bigger Cholesky
# updates on D_log at I=1840 (English). Smaller languages won't use all
# 8 cores but the partition has capacity.
#
# Model-ladder array. One task per (language, model_id) cell. The
# manifest below (cluster/sherlock/) is read in by SLURM_ARRAY_TASK_ID (1-indexed line).
#
# Submit with:   sbatch --array=1-N cluster/sherlock/glmer_ladder.slurm
#                where N = number of rows in sherlock/glmer_ladder_manifest.csv
#                (excluding the header).
#
# The manifest must live at sherlock/glmer_ladder_manifest.csv with columns:
#   task_id, lang_slug, model_id
#
# Output goes to $SCRATCH/standard_model_2/fits/glmer_ladder/...

set -euo pipefail

cd "$HOME/standard_model_2"

# scratch dirs
RESULTS_DIR="$SCRATCH/standard_model_2"
mkdir -p "$RESULTS_DIR/fits/glmer_ladder" "$RESULTS_DIR/logs"

# point R at scratch for fit outputs (and inputs)
export STANDARD_MODEL_FITS_DIR="$RESULTS_DIR/fits"

# make local data .rds files visible in scratch via symlink
for f in "$HOME/standard_model_2/fits/glmer_ladder"/data_*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"
  dst="$RESULTS_DIR/fits/glmer_ladder/$base"
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  ln -sfn "$f" "$dst" 2>/dev/null || true
done

# read the manifest line for this array task
MANIFEST="$HOME/standard_model_2/cluster/sherlock/glmer_ladder_manifest.csv"
if [ ! -f "$MANIFEST" ]; then
  echo "ERROR: missing $MANIFEST"
  exit 1
fi

# array task id is 1-indexed; line 1 of file is header; data lines start at 2
LINE=$((SLURM_ARRAY_TASK_ID + 1))
ROW=$(sed -n "${LINE}p" "$MANIFEST")
LANG_SLUG=$(echo "$ROW" | cut -d, -f2)
MODEL_ID=$(echo  "$ROW" | cut -d, -f3)

echo "===== SLURM task $SLURM_ARRAY_TASK_ID (job $SLURM_JOB_ID) ====="
echo "Host:      $(hostname)"
echo "LANG:      $LANG_SLUG"
echo "MODEL:     $MODEL_ID"
echo "Start:     $(date)"
echo

ml R

# Tell BLAS / Eigen how many threads to use. Cap at our SLURM allocation
# so we don't oversubscribe (oversubscription on Sherlock can get a job
# throttled). Setting both OPENBLAS and OMP env vars covers both common
# backends.
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
echo "Threads: $SLURM_CPUS_PER_TASK"

Rscript studies/glmer_ladder/02_fit_one.R "$LANG_SLUG" "$MODEL_ID"

echo
echo "End: $(date)"
