#!/usr/bin/env bash
# Submit the full planned set of fits to Sherlock's queue.
# Run from $HOME/standard_model_2 on the Sherlock login node.
#
# Edit the VARIANT_LIST/SIGMA_R_LIST below to taste.

set -euo pipefail
cd "$HOME/standard_model_2"

echo "Submitting longitudinal variants..."
for v in long_baseline long_2pl long_slopes long_2pl_slopes \
         long_2pl_slopes_nor; do
  jid=$(sbatch --parsable sherlock/long_fit.slurm "$v")
  echo "  $v -> job $jid"
done

echo
echo "Submitting cross-sectional sensitivity sweep..."
for sr in 0.3 0.53 0.8 1.2; do
  jid=$(sbatch --parsable sherlock/sensitivity_fit.slurm 2pl "$sr")
  echo "  2pl sigma_r=$sr -> job $jid"
done

echo
echo "Queue status:"
squeue -u "$USER" -o "%.10i %.12j %.2t %.10M %.20S"
