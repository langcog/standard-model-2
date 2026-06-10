#!/bin/bash
# Startup script for GCP VM running Stan fits.
# Installs R + cmdstanr + cmdstan + relevant R packages.
# Runs once at boot; subsequent logins just use the installed env.
#
# This is set as the `startup-script` metadata when the VM is created.
# Output goes to /var/log/syslog (and /var/log/cloud-init-output.log
# on systems where that's enabled).
#
# DISK REQUIREMENTS: provision at least 250 GB boot disk for fits
# with N >= 1M observations. Stan writes log_lik for every (chain, iter,
# obs) tuple to /tmp CSV (~17 MB per iter at N=2M, ~17 GB per chain
# for 1000 sampling iters, ~68 GB per fit total). The 100 GB default
# is NOT enough -- a previously-good fit's read_cmdstan_csv silently
# failed because /tmp filled, losing 11h of compute. Use:
#   gcloud compute instances create ... --boot-disk-size=250GB

set -euo pipefail

LOG=/var/log/sm2_setup.log
exec >>"$LOG" 2>&1

echo "[$(date)] === SM2 GCP VM setup starting ==="

# --- System packages -----------------------------------------------
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
  build-essential gfortran git curl wget rsync \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev \
  libjpeg-dev libcairo2-dev libharfbuzz-dev libfribidi-dev \
  cmake tmux htop

# --- R --------------------------------------------------------------
# Use Posit's binary package manager for fast installs (binary packages,
# not source compile). Adds R itself + a binary repo for packages.
apt-get install -y --no-install-recommends r-base r-base-dev

# Configure R to use Posit Package Manager for binary packages
mkdir -p /etc/R
UBUNTU_CODENAME=$(lsb_release -cs)
cat > /etc/R/Rprofile.site <<EOF
options(repos = c(
  CRAN = "https://packagemanager.posit.co/cran/__linux__/${UBUNTU_CODENAME}/latest"
))
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(),
        R.version["platform"], R.version["arch"], R.version["os"])))
EOF

# --- R packages -----------------------------------------------------
Rscript -e '
pkgs <- c("cmdstanr", "posterior", "loo", "dplyr", "tidyr", "ggplot2",
          "patchwork", "readr", "quantregGrowth", "rstan",
          "lme4", "Matrix", "tibble")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) {
  # cmdstanr is on stan-dev R-universe, not CRAN
  if ("cmdstanr" %in% need) {
    install.packages("cmdstanr",
                     repos = c("https://stan-dev.r-universe.dev",
                               getOption("repos")))
    need <- setdiff(need, "cmdstanr")
  }
  if (length(need)) install.packages(need, Ncpus = parallel::detectCores())
}
cat("R packages installed.\n")
'

# --- cmdstan (compiled toolchain for Stan) --------------------------
# Install to /opt so all users can find it. ~10 min for the cmdstan
# build the first time.
mkdir -p /opt/cmdstan
chmod 777 /opt/cmdstan
Rscript -e '
options(cmdstanr_cmdstan_install_path = "/opt/cmdstan")
cmdstanr::install_cmdstan(dir = "/opt/cmdstan", cores = parallel::detectCores())
'

# Find the actual installed version (cmdstanr installs the latest by
# default; version moves over time) and set env vars accordingly.
CMDSTAN_PATH=$(ls -d /opt/cmdstan/cmdstan-*/ 2>/dev/null | head -1 | sed 's:/$::')
echo "Detected cmdstan at: $CMDSTAN_PATH"

# Chown so non-root users can compile Stan models (compile writes
# intermediate object files into the cmdstan tree).
chmod -R 777 /opt/cmdstan

echo "CMDSTAN=$CMDSTAN_PATH" >> /etc/environment
cat > /etc/profile.d/cmdstan.sh <<EOF
export CMDSTAN="$CMDSTAN_PATH"
EOF

# Find and chown cmdstan to the default user so compiles work without sudo.
CMDSTAN_PATH=$(ls -d /opt/cmdstan/cmdstan-*/ 2>/dev/null | head -1 | sed 's:/$::')
chmod -R 777 /opt/cmdstan
# Set thread support in cmdstan make/local. Do NOT add -march=native or
# -DNDEBUG -- they trigger "double free or corruption (out)" on Stan 2.38
# with AMD EPYC and -O3 (default). Stan's default flags are fine.
cat > "$CMDSTAN_PATH/make/local" <<EOF
STAN_THREADS=true
EOF

# IMPORTANT thread-count note for users: when launching fits, set
# STAN_THREADS_PER_CHAIN to the number of PHYSICAL cores divided by
# number of chains. GCP n2d-* machines have 1 physical core per 2 vCPUs
# (SMT). Using all vCPUs as threads causes FPU contention and slows
# Stan ~3x. e.g. n2d-standard-32 = 16 physical cores -> 4 chains x 4
# threads/chain = 16 (full saturation, no SMT pressure).

echo "[$(date)] === SM2 GCP VM setup DONE ==="
touch /var/run/sm2_setup_complete
