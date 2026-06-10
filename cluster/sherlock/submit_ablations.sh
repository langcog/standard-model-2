#!/usr/bin/env bash
# Launch the ablation set against the lean baseline (long_slopes).
# Run from $HOME/standard_model_2 on the Sherlock login node.
#
# Each ablation toggles ONE component off/on relative to long_slopes:
#   long_baseline         drops per-child slopes        (RQ4 component)
#   long_fix_delta_slopes pins delta = 0                (RQ3: no acceleration)
#   long_free_s_slopes    frees start time s            (RQ2)
#   io_slopes             lean BabyView (input-observed)
#
# long_2pl_slopes already on disk; long_slopes is the reference.
# Norwegian gets the same treatment for cross-linguistic replication.

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

submit_io() {
  local variant="$1" dataset="$2"
  local jid
  jid=$(sbatch --parsable sherlock/io_fit.slurm "$variant" "$dataset")
  echo "  $variant/$dataset -> job $jid"
}

echo
echo "Submitting longitudinal ablations (English + Norwegian)..."
for ds in english norwegian; do
  submit_long long_baseline         "$ds"
  submit_long long_fix_delta_slopes "$ds"
  submit_long long_free_s_slopes    "$ds"
done

echo
echo "Submitting lean BabyView (input-observed)..."
submit_io io_slopes babyview

echo
echo "Queue status:"
squeue -u "$USER" -o "%.10i %.12j %.2t %.10M %.20S"
