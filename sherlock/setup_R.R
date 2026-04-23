## One-time R-environment setup for Sherlock.
##
## Run inside an interactive dev session (sh_dev -c 4) after `ml R`.
## Installs rstan and every package our scripts import.
##
## Sherlock's $HOME is user-writable; install.packages goes there by
## default. No sudo needed.

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "rstan", "posterior", "dplyr", "tidyr", "ggplot2", "tibble",
  "patchwork", "MASS", "arrow", "remotes", "purrr"
)

installed <- rownames(installed.packages())
needed    <- setdiff(cran_pkgs, installed)
if (length(needed) > 0) {
  message("Installing CRAN packages: ", paste(needed, collapse = ", "))
  install.packages(needed)
}

gh_pkgs <- list(
  wordbankr = "langcog/wordbankr",
  childesr  = "langcog/childesr"
)
for (pkg in names(gh_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s from GitHub...", pkg))
    remotes::install_github(gh_pkgs[[pkg]], upgrade = "never")
  }
}

cat("\nPackages installed. rstan version:\n")
print(packageVersion("rstan"))

## Pre-compile the Stan model so the first fit doesn't pay compile cost.
if (file.exists("model/stan/log_irt.stan")) {
  suppressPackageStartupMessages(library(rstan))
  rstan_options(auto_write = TRUE)
  message("Pre-compiling Stan models...")
  stan_model("model/stan/log_irt.stan")
  if (file.exists("model/stan/log_irt_long.stan")) {
    stan_model("model/stan/log_irt_long.stan")
  }
}

cat("\nSetup complete.\n")
