## Slim extraction: scalar summaries + per-item psi medians, with
## minimal RAM footprint. Uses cmdstanr's per-variable accessors so we
## never materialize the full draws data frame (which is ~30 GB for
## these fits and OOMs anything with <128 GB allocated).
##
## Usage on Sherlock:
##   Rscript sherlock/extract_summaries_slim.R <tag1> [<tag2> ...]
##
## Output for each tag:
##   <tag>.summary.rds   posterior summary for scalar params
##   <tag>.draws.rds     scalar parameter draws (small)
##   <tag>_psi.csv       per-item psi median (small)
##   <tag>.loo.rds       loo from log_lik (only if log_lik present)
##
## All outputs go to $SCRATCH/standard_model_2/summaries/.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_summaries_slim.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(loo)
})

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = file.path(Sys.getenv("SCRATCH"),
                                          "standard_model_2/fits"))
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCALAR_PARS <- c(
  "sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
  "pi_alpha", "delta", "s", "sigma_s", "beta_age",
  "sigma_lambda",
  "beta_c[1]", "beta_c[2]", "beta_c[3]", "beta_c[4]",
  "mu_r", "sigma_r"
)

extract_one <- function(tag) {
  in_path <- file.path(FITS_DIR, paste0(tag, ".rds"))
  if (!file.exists(in_path)) { cat("MISSING:", in_path, "\n"); return(invisible(NULL)) }
  cat(sprintf("\n== %s ==\n", tag))
  cat("Reading fit (cmdstanr lazy-load) ...\n")
  fit <- readRDS(in_path)
  stopifnot(inherits(fit, "CmdStanMCMC"))

  # Which scalar pars are actually present? Use $metadata for the
  # variable list without pulling any draws.
  all_vars <- fit$metadata()$stan_variables
  pars_present <- intersect(SCALAR_PARS, c(all_vars,
                                            sapply(seq_len(4), function(i) sprintf("beta_c[%d]", i))))
  # Drop indexed pars whose base isn't in metadata
  bases <- sub("\\[.*", "", pars_present)
  pars_present <- pars_present[bases %in% all_vars]
  cat(sprintf("Scalar pars present (%d): %s\n",
              length(pars_present), paste(pars_present, collapse=", ")))

  # Pull scalar draws ONLY (much smaller than full draws)
  d_scalar <- fit$draws(variables = unique(bases), format = "draws_df")
  d_scalar <- posterior::subset_draws(d_scalar, variable = pars_present)

  summary_tbl <- posterior::summarise_draws(
    d_scalar, "mean", "median", "sd",
    q025 = ~ stats::quantile(.x, 0.025, names = FALSE),
    q975 = ~ stats::quantile(.x, 0.975, names = FALSE),
    "ess_bulk", "rhat"
  )
  saveRDS(summary_tbl, file.path(OUT_DIR, paste0(tag, ".summary.rds")))

  draws_small <- as.data.frame(d_scalar)
  saveRDS(draws_small, file.path(OUT_DIR, paste0(tag, ".draws.rds")))
  cat(sprintf("Wrote .summary.rds and .draws.rds (%d rows, %d pars).\n",
              nrow(draws_small), length(pars_present)))

  # ---- psi medians, one at a time so we never materialize all 200 ---- #
  if ("psi" %in% all_vars) {
    cat("Computing per-item psi medians ...\n")
    psi_draws <- fit$draws(variables = "psi", format = "draws_matrix")
    psi_med <- apply(psi_draws, 2, median)
    psi_mean <- apply(psi_draws, 2, mean)
    psi_csv <- data.frame(
      jj = seq_along(psi_med),
      psi_mean = psi_mean,
      psi_median = psi_med
    )
    write.csv(psi_csv, file.path(OUT_DIR, paste0(tag, "_psi.csv")),
              row.names = FALSE)
    rm(psi_draws, psi_med, psi_mean); gc(verbose = FALSE)
    cat(sprintf("Wrote _psi.csv (%d items)\n", nrow(psi_csv)))
  }

  # ---- LOO if log_lik present (small per-variable draws, then collapse) ----
  if ("log_lik" %in% all_vars) {
    cat("Computing LOO ...\n")
    loo_obj <- tryCatch(fit$loo(cores = 2),
                        error = function(e) { cat("LOO error:", conditionMessage(e), "\n"); NULL })
    if (!is.null(loo_obj)) {
      saveRDS(loo_obj, file.path(OUT_DIR, paste0(tag, ".loo.rds")))
      cat(sprintf("Wrote .loo.rds (elpd = %.1f +- %.1f).\n",
                  loo_obj$estimates["elpd_loo", "Estimate"],
                  loo_obj$estimates["elpd_loo", "SE"]))
    }
  } else {
    cat("No log_lik; skipping LOO.\n")
  }

  rm(fit, d_scalar, draws_small); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}

cat("\nAll done. Summaries in", OUT_DIR, "\n")
