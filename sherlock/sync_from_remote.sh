#!/usr/bin/env bash
# Pull fit outputs + logs back from Sherlock to the local repo.
# Relies on `.env.local` for SHERLOCK_HOST and SHERLOCK_FITS_DIR.
# Auth is whatever your ~/.ssh/config says (SSH key + DUO 2FA).
#
# Usage:   ./sherlock/sync_from_remote.sh [pattern]
#          pattern = optional glob passed to rsync (default: '*')
#
# Examples:
#   ./sherlock/sync_from_remote.sh
#   ./sherlock/sync_from_remote.sh 'long_2pl_slopes_norwegian*'

set -euo pipefail

cd "$(dirname "$0")/.."

if [ -f .env.local ]; then
  # shellcheck disable=SC1091
  source .env.local
else
  echo "Warning: .env.local not found. Using defaults (SHERLOCK_HOST=sherlock)."
  SHERLOCK_HOST="${SHERLOCK_HOST:-sherlock}"
  SHERLOCK_FITS_DIR="${SHERLOCK_FITS_DIR:-\$SCRATCH/standard_model_2/fits}"
fi

PATTERN="${1:-*}"

# Use single quotes on the remote path so $SCRATCH expands on Sherlock.
echo "Pulling $PATTERN from $SHERLOCK_HOST:$SHERLOCK_FITS_DIR ..."
# -L dereferences symlinks: the Sherlock SLURM scripts symlink input
# bundles into $SCRATCH/.../fits, and we want their actual contents
# pulled back, not the (broken-on-laptop) symlink itself.
rsync -avzL --progress \
    "${SHERLOCK_HOST}:${SHERLOCK_FITS_DIR}/${PATTERN}" \
    fits/

echo
echo "Done. Files in fits/:"
ls -la fits/ | tail -20
