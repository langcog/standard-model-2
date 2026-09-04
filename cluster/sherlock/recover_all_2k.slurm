#!/usr/bin/env bash
#SBATCH --job-name=recover6
#SBATCH --output=logs/recover6_%j.out
#SBATCH --error=logs/recover6_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --exclusive
#SBATCH --partition=mcfrank
#
# Recover all six imputed _2k refits from their cmdstan CSVs in parallel on
# the owned node. Why one whole-node job instead of six small ones: small
# (3-CPU) jobs on the mcfrank partition sat PENDING behind owners-partition
# jobs for 20+ min on 2026-09-03 -- Slurm only preempts lower-tier jobs when
# the higher-tier job cannot otherwise be scheduled, and it kept waiting for a
# 3-CPU gap. A 16-24 CPU request forces preemption (the refits themselves
# started within minutes that way). Each recovery is a grep|cut stream over
# ~200 GB (I/O-bound, ~1 GB RAM), so six in parallel is fine.
#
# Usage: sbatch cluster/sherlock/recover_all_2k.slurm

set -uo pipefail
cd "$HOME/standard_model_2"
export STANDARD_MODEL_FITS_DIR="$SCRATCH/standard_model_2/fits"
ml R
mkdir -p logs

TAGS=(no_freq_slopes_sigmaR_0p35_2k no_freq_slopes_sigmaR_0p44_2k no_freq_slopes_sigmaR_0p58_2k
      no_freq_slopes_norwegian_sigmaR_0p35_2k no_freq_slopes_norwegian_sigmaR_0p44_2k
      no_freq_slopes_norwegian_sigmaR_0p58_2k)

echo "===== recover6  host $(hostname)  $(date) ====="
for t in "${TAGS[@]}"; do
  ( Rscript cluster/sherlock/recover_from_csvs.R "$t" > "logs/recover_${t}.log" 2>&1 \
      && echo "$(date '+%T') DONE $t" || echo "$(date '+%T') FAILED $t (see logs/recover_${t}.log)" ) &
done
wait
echo "End: $(date)"
