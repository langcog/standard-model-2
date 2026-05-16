## Scalars-only extraction. Skips psi and log_lik to fit in <16 GB.
## Use sherlock/extract_delta_j.R separately for psi medians (needs more
## memory because loading psi as draws_matrix instantiates ~12 GB).
##
## Usage: Rscript sherlock/extract_scalars_only.R <tag1> [<tag2> ...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_scalars_only.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = file.path(Sys.getenv("SCRATCH"),
                                          "standard_model_2/fits"))
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCALAR_PARS <- c(
  "sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
  "pi_alpha", "delta", "s", "sigma_s",
  "sigma_lambda",
  "beta_c[1]", "beta_c[2]", "beta_c[3]", "beta_c[4]",
  "mu_r", "sigma_r"
)

extract_one <- function(tag) {
  in_path <- file.path(FITS_DIR, paste0(tag, ".rds"))
  if (!file.exists(in_path)) { cat("MISSING:", in_path, "\n"); return(NULL) }
  cat(sprintf("\n== %s ==\n", tag))
  fit <- readRDS(in_path)
  all_vars <- fit$metadata()$stan_variables
  pars_present <- intersect(SCALAR_PARS, c(all_vars,
                  sapply(seq_len(4), function(i) sprintf("beta_c[%d]", i))))
  bases <- sub("\\[.*", "", pars_present)
  pars_present <- pars_present[bases %in% all_vars]
  cat(sprintf("Scalar pars present (%d): %s\n",
              length(pars_present), paste(pars_present, collapse=", ")))

  d_scalar <- fit$draws(variables = unique(bases), format = "draws_df")
  d_scalar <- posterior::subset_draws(d_scalar, variable = pars_present)

  summary_tbl <- posterior::summarise_draws(
    d_scalar, "mean", "median", "sd",
    q025 = ~ stats::quantile(.x, 0.025, names = FALSE),
    q975 = ~ stats::quantile(.x, 0.975, names = FALSE),
    "ess_bulk", "rhat"
  )
  saveRDS(summary_tbl, file.path(OUT_DIR, paste0(tag, ".summary.rds")))
  saveRDS(as.data.frame(d_scalar),
          file.path(OUT_DIR, paste0(tag, ".draws.rds")))
  cat(sprintf("Wrote .summary.rds and .draws.rds (%d draws, %d pars).\n",
              nrow(d_scalar), length(pars_present)))

  rm(fit, d_scalar); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}
cat("\nAll done.\n")
