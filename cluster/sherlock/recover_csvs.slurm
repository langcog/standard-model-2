#!/usr/bin/env bash
#SBATCH --job-name=recover
#SBATCH --output=logs/recover_%j.out
#SBATCH --error=logs/recover_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=normal
#
# Stream a finished fit's slim summary + draws out of its cmdstan CSVs with
# recover_from_csvs.R (grep | cut over ~200 GB per big fit: I/O-bound, RAM
# stays tiny). For fits whose post-sampling R step died (OOM) before the
# inline recovery in long_fit.slurm could run -- or any time the CSVs are
# the only surviving artifact.
#
# Usage:  sbatch cluster/sherlock/recover_csvs.slurm <tag>
#         e.g. sbatch cluster/sherlock/recover_csvs.slurm no_freq_slopes_sigmaR_0p35_2k
# Reads   $SCRATCH/standard_model_2/fits/csvs_<tag>/
# Writes  $SCRATCH/standard_model_2/fits/summaries/<tag>.{summary,draws}.rds (+ _psi.csv)

set -euo pipefail
TAG="${1:?usage: recover_csvs.slurm <tag>}"
cd "$HOME/standard_model_2"
export STANDARD_MODEL_FITS_DIR="$SCRATCH/standard_model_2/fits"
ml R

echo "===== recover $TAG  host $(hostname)  $(date) ====="
# idempotent: several recovery routes may be queued for the same tag
if [ -s "$STANDARD_MODEL_FITS_DIR/summaries/$TAG.summary.rds" ]; then
  echo "SKIP: summary already exists"; exit 0
fi
Rscript cluster/sherlock/recover_from_csvs.R "$TAG"
echo "End: $(date)"
