#!/usr/bin/env bash
# Convenience: run a command on Sherlock inside the project directory.
#
# Usage:  ./sherlock/remote_exec.sh "<shell command>"
# Example:
#   ./sherlock/remote_exec.sh "git pull && squeue -u \$USER"
#   ./sherlock/remote_exec.sh "sbatch sherlock/long_fit.slurm long_2pl_slopes norwegian"

set -euo pipefail
cd "$(dirname "$0")/.."

if [ -f .env.local ]; then
  # shellcheck disable=SC1091
  source .env.local
fi

SHERLOCK_HOST="${SHERLOCK_HOST:-sherlock}"
SHERLOCK_PROJECT_DIR="${SHERLOCK_PROJECT_DIR:-\$HOME/standard_model_2}"

CMD="${*:-}"
if [ -z "$CMD" ]; then
  echo "Usage: $0 \"<remote shell command>\""
  exit 1
fi

# Quote the whole remote command so $HOME / $SCRATCH expand on the server
ssh "$SHERLOCK_HOST" "cd ${SHERLOCK_PROJECT_DIR} && ${CMD}"
