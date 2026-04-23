## One-time R-environment setup for Sherlock.
##
## Run inside an interactive dev session (sh_dev -c 4) after `ml R`.
## Installs rstan and every package our scripts import.
##
## Sherlock's $HOME is user-writable; install.packages goes there by
## default. No sudo needed.

options(repos = c(CRAN = "https://cloud.r-project.org"))

# On Sherlock the default system library isn't writable; make sure
# the user library exists and is first on the path.
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) {
  rv <- paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  user_lib <- file.path("~/R",
                        paste0(R.version$platform, "-library"), rv)
  user_lib <- path.expand(user_lib)
}
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))
cat("User library:", user_lib, "\n")
cat("Library paths:\n"); print(.libPaths())

cran_pkgs <- c(
  "rstan", "posterior", "dplyr", "tidyr", "ggplot2", "tibble",
  "patchwork", "MASS", "arrow", "remotes", "purrr"
)

installed <- rownames(installed.packages())
needed    <- setdiff(cran_pkgs, installed)
if (length(needed) > 0) {
  message("Installing CRAN packages to ", user_lib, ": ",
          paste(needed, collapse = ", "))
  install.packages(needed, lib = user_lib)
}

gh_pkgs <- list(
  wordbankr = "langcog/wordbankr",
  childesr  = "langcog/childesr"
)
for (pkg in names(gh_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s from GitHub...", pkg))
    remotes::install_github(gh_pkgs[[pkg]], upgrade = "never",
                            lib = user_lib)
  }
}

cat("\nPackages installed. rstan version:\n")
print(packageVersion("rstan"))

## Pre-compile the Stan models so the first fit doesn't pay compile cost.
## This requires ~8 GB RAM per model; if the dev session doesn't have that
## memory (OOM -> "Killed signal terminated program cc1plus"), skip it and
## the first SLURM fit will compile (SLURM scripts reserve more RAM).
maybe_compile <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  message("Pre-compiling ", path, " ...")
  tryCatch({
    suppressPackageStartupMessages(library(rstan))
    rstan_options(auto_write = TRUE)
    stan_model(path)
    message("  OK")
  }, error = function(e) {
    message("  SKIPPED (", conditionMessage(e), ")")
    message("  This is usually an OOM in the dev session. It's fine:")
    message("  the first SLURM fit will compile the model and cache it.")
  })
}
maybe_compile("model/stan/log_irt.stan")
maybe_compile("model/stan/log_irt_long.stan")

cat("\nSetup complete.\n")
