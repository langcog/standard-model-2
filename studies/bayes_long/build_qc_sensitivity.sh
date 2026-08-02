#!/usr/bin/env bash
# Build the exclusion-sensitivity bundles for the SI.
#
# WHY. A reviewer asked to see the analysis "with no exclusions" and "with several decline
# thresholds", on the concern that the filter could preferentially remove noisy or
# non-monotonic trajectories and thereby INFLATE estimated acceleration. That worry has a
# direction, so it is directly testable: if kappa is unchanged with the filter off, the
# criticism is answered; if kappa falls, it is not.
#
# The concern also bites harder than it first appeared. The exclusion rate quoted in the
# Methods is now computed at the main-text 3+ threshold, where Marchman loses 16.5% of
# eligible children -- not the 2.6% the 2+ table showed. Only 5.7% of Marchman
# administrations are actually flagged, but at a 3+ threshold a child who loses one drops
# out entirely, so the filter and the inclusion rule interact.
#
# Four conditions, all at MIN_ADMINS=3 to match the main text:
#   none    filter off entirely                      -> _a3_qcnone
#   loose   only extreme outliers                    -> _a3_qcloose
#   main    the reported setting (already built)     -> _a3
#   tight   aggressive                               -> _a3_qctight
#
# Usage:  bash studies/bayes_long/build_qc_sensitivity.sh
# Then:   fit m3 on each variant (see the sbatch loop printed at the end).
set -euo pipefail
cd "$(dirname "$0")/../.."

echo "=== none: no QC filter at all ==="
MIN_ADMINS=3 QC_OFF=1 QC_TAG=_qcnone Rscript studies/bayes_long/00_prepare_bundles.R

echo "=== loose: crater 40% below peak, jump 60 pts/mo ==="
MIN_ADMINS=3 QC_REL_TOL=0.40 QC_RATE_MAX=0.60 QC_TAG=_qcloose \
  Rscript studies/bayes_long/00_prepare_bundles.R

echo "=== tight: crater 15% below peak, jump 25 pts/mo ==="
MIN_ADMINS=3 QC_REL_TOL=0.15 QC_RATE_MAX=0.25 QC_TAG=_qctight \
  Rscript studies/bayes_long/00_prepare_bundles.R

echo
echo "built. sample sizes by condition:"
Rscript -e '
for (tag in c("", "_qcnone", "_qcloose", "_qctight")) {
  n <- vapply(c("thal","smith","marchman","norwegian","japanese"), function(s) {
    f <- sprintf("fits/bayes_long/bundle_%s_a3%s.rds", s, tag)
    if (file.exists(f)) readRDS(f)$stan_data$I else NA_integer_ }, integer(1))
  cat(sprintf("  %-9s %s  total %d\n", ifelse(tag=="","main",sub("_qc","",tag)),
              paste(sprintf("%5d", n), collapse=""), sum(n, na.rm=TRUE)))
}' 2>&1 | grep -viE "built under|Warning"

cat <<'EOF'

Now fit M3 on each variant (main is already fitted):
  for tag in _qcnone _qcloose _qctight; do
    for s in japanese smith marchman thal; do
      sbatch -p owners --requeue --job-name=qc${tag}_${s} \
        studies/bayes_long/fit.slurm ${s}_a3${tag} m3
    done
    sbatch -p mcfrank --cpus-per-task=24 --export=ALL,STAN_THREADS=6 \
      --job-name=qc${tag}_norwegian studies/bayes_long/fit.slurm norwegian_a3${tag} m3
  done
EOF
