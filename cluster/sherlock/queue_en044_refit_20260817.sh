#!/usr/bin/env bash
# When the 1000-iter EN sigma_r=0.44 job (39218380) leaves the queue and its
# summary is recovered (by auto_recover_imputed_20260815.sh), submit the
# 2000-iter convergence refit warm-started from it. One-shot; nohup on login.
cd "$HOME/standard_model_2"
SUMM="$SCRATCH/standard_model_2/fits/summaries/no_freq_slopes_sigmaR_0p44.summary.rds"
LOG="$SCRATCH/standard_model_2/logs/queue_en044_refit.log"
while squeue -h -j 39218380 2>/dev/null | grep -q .; do sleep 600; done
echo "$(date '+%F %T') 39218380 left queue; waiting for recovered summary" >> "$LOG"
for i in $(seq 1 60); do [ -s "$SUMM" ] && break; sleep 300; done
if [ -s "$SUMM" ]; then INIT="STAN_INIT_FROM=no_freq_slopes_sigmaR_0p44,"; echo "$(date '+%F %T') summary present -> warm start" >> "$LOG"
else INIT=""; echo "$(date '+%F %T') no summary after 5h -> cold start" >> "$LOG"; fi
sbatch -p mcfrank -t 96:00:00 --mem=64G \
  --export=ALL,STAN_SIGMA_R_OVERRIDE=0.44,${INIT}STAN_ITER=2000,STAN_WARMUP=1000,STAN_ADAPT_DELTA=0.97,STAN_TAG_SFX=_2k \
  cluster/sherlock/long_fit.slurm no_freq_slopes english >> "$LOG" 2>&1
echo "$(date '+%F %T') submitted EN 0.44 _2k refit" >> "$LOG"
