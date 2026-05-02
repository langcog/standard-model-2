#!/usr/bin/env bash
# Submit LMM (linear-in-age) ablations for direct comparison against the
# log-linear baseline. The LMM uses a different Stan file
# (log_irt_long_lmm.stan); fit_longitudinal.R routes by variant name.
#
# Variants:
#   long_lmm         pure LMM, no per-child slopes
#   long_lmm_slopes  LMM + per-child slopes (matched to long_slopes)

set -euo pipefail
cd "$HOME/standard_model_2"

echo "Linking input bundles into \$SCRATCH..."
./sherlock/sync_inputs.sh

submit_long() {
  local variant="$1" dataset="$2"
  local jid
  jid=$(sbatch --parsable sherlock/long_fit.slurm "$variant" "$dataset")
  echo "  $variant/$dataset -> job $jid"
}

echo
echo "Submitting LMM vs log-linear comparison fits..."
for ds in english norwegian; do
  submit_long long_lmm_slopes "$ds"
done

echo
echo "Queue status:"
squeue -u "$USER" -o "%.10i %.18j %.2t %.10M %.20S"
