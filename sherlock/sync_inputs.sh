#!/usr/bin/env bash
# Symlink all repo-side input bundles into $SCRATCH/standard_model_2/fits/
# so SLURM jobs find their data. Run once before sbatch — calling this
# from inside multiple concurrent SLURM jobs hits a race on `ln -f`.

set -euo pipefail

SRC="$HOME/standard_model_2/fits"
DST="$SCRATCH/standard_model_2/fits"
mkdir -p "$DST"

count=0
for f in "$SRC"/*.rds; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"
  dst="$DST/$base"
  # Skip if already a correct symlink
  if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$f" ]; then continue; fi
  rm -f "$dst" 2>/dev/null || true
  ln -s "$f" "$dst"
  count=$((count + 1))
done

echo "sync_inputs: linked $count files into $DST"
ls -la "$DST"
