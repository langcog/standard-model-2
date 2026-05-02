#!/usr/bin/env bash
# Submit the new pieces of the M0..M5 nested family on English.
# M2/M3/M4 are already on disk under their established names
# (long_baseline / long_slopes / long_class_beta_slopes); we only need
# M0, M1, M5.
#
# Mapping:
#   M0  long_m0                     no time, no freq, no slopes, no 2PL
#   M1  long_m1                     unit time + freq, delta pinned, no slopes/2PL
#   M2  long_baseline               adds free delta
#   M3  long_slopes                 + per-child zeta
#   M4  long_class_beta_slopes      + class-specific beta_c
#   M5  long_m5                     + 2PL (lambda_j)
#
# Plus no_freq_slopes as a robustness comparison to M3 (drops log p_j).

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
echo "Submitting nested-family fits (English only)..."
submit_long long_m0              english
submit_long long_m1              english
submit_long long_m5              english
submit_long long_no_freq_slopes  english

echo
echo "Already on disk (skipping):"
echo "  M2 = long_baseline"
echo "  M3 = long_slopes"
echo "  M4 = long_class_beta_slopes"
echo
echo "Queue status:"
squeue -u "$USER" -o "%.10i %.22j %.2t %.10M %.20S"
