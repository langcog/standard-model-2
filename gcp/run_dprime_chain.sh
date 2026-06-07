#!/bin/bash
# Run the D -> D' chain for one dataset on a GCP node, sequentially, on
# the deduped (one-admin-per-child-age) longitudinal bundle.
#
#   Node A:  gcp/run_dprime_chain.sh english
#   Node B:  gcp/run_dprime_chain.sh norwegian
#
# Fits, in order:
#   1. long_no_freq_slopes        (D / M_best)  -- free rho_xi_zeta
#   2. long_no_freq_slopes_dprime (D')          -- rho_xi_zeta pinned 0,
#                                                   + gamma_in * log_r_dev
#                                                   in the slope.
#
# Both fits land on the SAME corrected bundle so LOO(D) vs LOO(D') is a
# clean within-language comparison. D' is COLD-started (it's a new model;
# we want an init-independent baseline before any warm-starting in the
# sigma_r sweep).
#
# gamma_in = Cov(xi, kappa)/sigma_r^2, so it is sensitive to the external
# sigma_r pin. After this chain converges, map that dependence with the
# sigma_r sweep (see gcp/run_dprime_sigmar_sweep.sh).
#
# Run detached so the node can be left alone:
#   nohup bash gcp/run_dprime_chain.sh english > dprime_english.log 2>&1 &

set -euo pipefail

DATASET="${1:-english}"
cd "$HOME/standard_model_2"

# Match the env used for the EN M_best baseline (wait_then_no_mbest.sh).
# Override per-node if a node has more/fewer vCPUs.
export STAN_ITER="${STAN_ITER:-2000}"
export STAN_WARMUP="${STAN_WARMUP:-1000}"
export STAN_ADAPT_DELTA="${STAN_ADAPT_DELTA:-0.95}"
export STAN_THREADS_PER_CHAIN="${STAN_THREADS_PER_CHAIN:-16}"
export STAN_MAX_TREEDEPTH="${STAN_MAX_TREEDEPTH:-10}"

echo "==== [$DATASET] D (M_best, long_no_freq_slopes) starting $(date -u) ===="
bash gcp/run_fit.sh long_no_freq_slopes "$DATASET"

echo "==== [$DATASET] D' (long_no_freq_slopes_dprime) starting $(date -u) ===="
bash gcp/run_fit.sh long_no_freq_slopes_dprime "$DATASET"

echo "==== [$DATASET] D -> D' chain DONE $(date -u) ===="