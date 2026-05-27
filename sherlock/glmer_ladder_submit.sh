#!/usr/bin/env bash
## Build the glmer_ladder manifest and submit a SLURM array.
##
## Usage:
##   bash sherlock/glmer_ladder_submit.sh            # just build the manifest, print sbatch lines
##   bash sherlock/glmer_ladder_submit.sh smoke      # submit ONE small task (Japanese B_log)
##   bash sherlock/glmer_ladder_submit.sh nor_big    # submit Norwegian D_log alone (timing test)
##   bash sherlock/glmer_ladder_submit.sh all        # submit the full 49-cell array
##
## Manifest at sherlock/glmer_ladder_manifest.csv, one row per (lang, model).

set -euo pipefail

cd "$(dirname "$0")/.."

LANGS=(english_american norwegian finnish french_quebecois japanese spanish_mexican french_french)
MODELS=(A B_log B_lin C_log C_lin D_log D_lin)

MANIFEST="sherlock/glmer_ladder_manifest.csv"
echo "task_id,lang_slug,model_id" > "$MANIFEST"

i=1
SMOKE_TASK=""
NOR_BIG_TASK=""
for L in "${LANGS[@]}"; do
  for M in "${MODELS[@]}"; do
    echo "$i,$L,$M" >> "$MANIFEST"
    # remember which task indices are our two tests
    if [ "$L" = "japanese" ] && [ "$M" = "B_log" ]; then SMOKE_TASK="$i"; fi
    if [ "$L" = "norwegian" ] && [ "$M" = "D_log" ]; then NOR_BIG_TASK="$i"; fi
    i=$((i + 1))
  done
done
N=$((i - 1))

echo "Manifest with $N tasks written to $MANIFEST"
echo "  smoke (Japanese B_log)       = task $SMOKE_TASK"
echo "  nor_big (Norwegian D_log)    = task $NOR_BIG_TASK"
echo

ACTION="${1:-}"
case "$ACTION" in
  smoke)
    echo "Submitting infrastructure smoke test (Japanese B_log)..."
    sbatch --array="$SMOKE_TASK" sherlock/glmer_ladder.slurm
    ;;
  nor_big)
    echo "Submitting Norwegian D_log timing/size test..."
    sbatch --array="$NOR_BIG_TASK" sherlock/glmer_ladder.slurm
    ;;
  all)
    echo "Submitting full 49-cell array..."
    sbatch --array="1-$N" sherlock/glmer_ladder.slurm
    ;;
  "")
    echo "No action; manifest built only. Submit options:"
    echo "  bash sherlock/glmer_ladder_submit.sh smoke      # task $SMOKE_TASK"
    echo "  bash sherlock/glmer_ladder_submit.sh nor_big    # task $NOR_BIG_TASK"
    echo "  bash sherlock/glmer_ladder_submit.sh all        # full 1-$N"
    ;;
  *)
    echo "Unknown action: $ACTION"
    echo "Use: smoke | nor_big | all   (or no arg)"
    exit 1
    ;;
esac
