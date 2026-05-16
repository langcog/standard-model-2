#!/bin/bash
# Startup script for GCP VM running Stan fits.
# Installs R + cmdstanr + cmdstan + relevant R packages.
# Runs once at boot; subsequent logins just use the installed env.
#
# This is set as the `startup-script` metadata when the VM is created.
# Output goes to /var/log/syslog (and /var/log/cloud-init-output.log
# on systems where that's enabled).

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
chmod 755 /opt/cmdstan
Rscript -e '
options(cmdstanr_cmdstan_install_path = "/opt/cmdstan")
cmdstanr::install_cmdstan(dir = "/opt/cmdstan", cores = parallel::detectCores())
'

# Make cmdstanr know where to find the install
echo 'CMDSTAN=/opt/cmdstan/cmdstan-2.36.0' >> /etc/environment
echo 'export CMDSTAN=/opt/cmdstan/cmdstan-2.36.0' >> /etc/profile.d/cmdstan.sh

echo "[$(date)] === SM2 GCP VM setup DONE ==="
touch /var/run/sm2_setup_complete
