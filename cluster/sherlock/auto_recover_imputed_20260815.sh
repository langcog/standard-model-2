#!/usr/bin/env bash
# One-shot watcher for the 2026-08-15 imputed refits (jobs launched without
# STAN_SKIP_SAVE_OBJECT, so they OOM at the save step after clean sampling).
# Waits for each job to leave the queue, then streams its slim summary from
# the CSVs. Runs on the login node under nohup; exits when all are done.
# Log: $SCRATCH/standard_model_2/logs/auto_recover.log
cd "$HOME/standard_model_2"
export STANDARD_MODEL_FITS_DIR="$SCRATCH/standard_model_2/fits"
ml R >/dev/null 2>&1
declare -A TAGS=(
  [39218376]=no_freq_slopes_sigmaR_0p35
  [39218379]=no_freq_slopes_norwegian_sigmaR_0p35
  [39218380]=no_freq_slopes_sigmaR_0p44
  [39218382]=no_freq_slopes_sigmaR_0p58
  [39218385]=no_freq_slopes_norwegian_sigmaR_0p58
)
pending=("${!TAGS[@]}")
while [ ${#pending[@]} -gt 0 ]; do
  still=()
  for j in "${pending[@]}"; do
    if squeue -h -j "$j" 2>/dev/null | grep -q .; then
      still+=("$j")
    else
      tag="${TAGS[$j]}"
      state=$(sacct -j "$j" -n -X --format=State 2>/dev/null | head -1 | tr -d ' ')
      echo "$(date '+%F %T') job $j ($tag) left queue: state=$state -> recovering"
      if [ -s "$STANDARD_MODEL_FITS_DIR/summaries/$tag.summary.rds" ]; then
        echo "$(date '+%F %T')   summary already present, skipping"
      else
        Rscript cluster/sherlock/recover_from_csvs.R "$tag" \
          > "$SCRATCH/standard_model_2/logs/recover_${tag}.log" 2>&1 \
          && echo "$(date '+%F %T')   recovered $tag" \
          || echo "$(date '+%F %T')   RECOVERY FAILED for $tag (see recover_${tag}.log)"
      fi
    fi
  done
  pending=("${still[@]}")
  [ ${#pending[@]} -gt 0 ] && sleep 600
done
echo "$(date '+%F %T') all imputed refits recovered; done."
