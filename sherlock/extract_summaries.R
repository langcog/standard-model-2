## Sherlock-side: extract small summary objects from a (4 GB) fit, so we
## can sync just the summaries home instead of the whole fit.
##
## For each input fit <tag>.rds, writes:
##   <tag>.summary.rds  posterior summary table for headline scalars
##                      (mean, sd, 2.5%, 50%, 97.5%, Rhat, n_eff)
##   <tag>.draws.rds    a tibble of named scalar parameters' full draws
##                      (small; used for plotting forests & joint
##                       posteriors in R locally)
##   <tag>.loo.rds      loo() output computed from the model's log_lik;
##                      omitted if log_lik is not present.
##
## Usage on Sherlock:
##   Rscript sherlock/extract_summaries.R fit_tag1 fit_tag2 ...
##
## Each call processes its arguments in order; output goes to
## $SCRATCH/standard_model_2/summaries/.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_summaries.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(loo)
})

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = file.path(Sys.getenv("SCRATCH"),
                                          "standard_model_2/fits"))
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Scalar parameters worth keeping full draws for (small; used for forest
# plots, joint posterior contours, and PPC sampling locally).
SCALAR_PARS <- c(
  "sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
  "pi_alpha", "delta", "s", "beta_age",
  "sigma_lambda",
  "beta_c[1]", "beta_c[2]", "beta_c[3]", "beta_c[4]",
  "mu_r", "sigma_r", "beta_react", "sigma_within",
  "gamma_input", "sigma_zeta_resid", "sigma_zeta_marginal",
  "input_share_zeta",
  "gamma_rt", "mu_rtslope", "sigma_rtslope", "sigma_lwl"
)

extract_one <- function(tag) {
  in_path  <- file.path(FITS_DIR, paste0(tag, ".rds"))
  if (!file.exists(in_path)) {
    cat(sprintf("MISSING: %s\n", in_path)); return(invisible(NULL))
  }
  cat(sprintf("\n== %s ==\n", tag))
  cat(sprintf("Reading %s ...\n", in_path))
  fit <- readRDS(in_path)
  d   <- as_draws_df(fit)

  # ---- Summary: scalar posteriors ---- #
  pars_present <- intersect(SCALAR_PARS, names(d))
  cat(sprintf("Scalar pars present (%d): %s\n",
              length(pars_present),
              paste(pars_present, collapse = ", ")))

  summary_tbl <- summary(fit, pars = pars_present)$summary
  saveRDS(summary_tbl,
          file = file.path(OUT_DIR, paste0(tag, ".summary.rds")))

  # ---- Draws: scalar parameters only (small) ---- #
  draws_small <- d[, c("draw__", "chain__", "iteration__",
                       intersect(pars_present, names(d)))]
  saveRDS(draws_small,
          file = file.path(OUT_DIR, paste0(tag, ".draws.rds")))
  cat(sprintf("Wrote summary + draws (%d rows, %d pars).\n",
              nrow(draws_small), ncol(draws_small) - 3))

  # ---- LOO: ELPD + diagnostics + per-obs estimates ---- #
  if ("log_lik[1]" %in% names(d)) {
    cat("Computing LOO from log_lik ...\n")
    log_lik_arr <- extract_log_lik(fit, parameter_name = "log_lik",
                                    merge_chains = FALSE)
    r_eff <- relative_eff(exp(log_lik_arr))
    loo_obj <- loo(log_lik_arr, r_eff = r_eff, cores = 4)
    saveRDS(loo_obj,
            file = file.path(OUT_DIR, paste0(tag, ".loo.rds")))
    cat(sprintf("Wrote loo (elpd = %.1f +- %.1f).\n",
                loo_obj$estimates["elpd_loo", "Estimate"],
                loo_obj$estimates["elpd_loo", "SE"]))
  } else {
    cat("No log_lik in this fit; skipping LOO.\n")
  }

  # Free the giant fit object before processing the next one.
  rm(fit, d, draws_small)
  gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) {
             cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e)))
           })
}

cat("\nAll done. Summaries in", OUT_DIR, "\n")
