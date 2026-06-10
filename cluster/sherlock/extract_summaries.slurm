#!/usr/bin/env bash
#SBATCH --job-name=extract_summaries
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=normal

# Usage:
#   sbatch sherlock/extract_summaries.slurm <tag1> [<tag2> ...]
#
# Runs sherlock/extract_summaries.R on each fit tag, producing small
# summary RDS files in $SCRATCH/standard_model_2/summaries/.
#
# 24G RAM is enough for one 4 GB rstan fit at a time plus the loo
# computation. The script processes tags sequentially so memory only
# needs to hold one fit at a time.

set -euo pipefail
cd "$HOME/standard_model_2"

module load R

if [ "$#" -eq 0 ]; then
  echo "Usage: sbatch sherlock/extract_summaries.slurm <tag1> [<tag2> ...]"
  exit 1
fi

echo "Extracting summaries for tags: $*"
Rscript sherlock/extract_summaries.R "$@"
