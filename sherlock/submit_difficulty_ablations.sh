#!/usr/bin/env bash
# Submit the difficulty-side ablation set: variants that probe the
# beta_j = psi_j - log p_j - log H term, complementing the existing
# ability-side ablations (slopes / pin delta / free s).
#
# Variants:
#   long_no_class_slopes    drop class hierarchy on psi_j (single global mu, tau)
#   long_class_beta_slopes  per-class slope beta_c on log p_j (free, prior N(1, 0.5))
#
# Both variants run on English and Norwegian for cross-linguistic
# replication, mirroring the ability-side ablations.

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
echo "Submitting difficulty-side ablations (English + Norwegian)..."
for ds in english norwegian; do
  submit_long long_no_class_slopes    "$ds"
  submit_long long_class_beta_slopes  "$ds"
done

echo
echo "Queue status:"
squeue -u "$USER" -o "%.10i %.12j %.2t %.10M %.20S"
