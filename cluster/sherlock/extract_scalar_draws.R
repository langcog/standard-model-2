## Extract scalar parameter draws ONE VARIABLE AT A TIME to stay under
## the 16 GB memory limit. Combines them into a single draws_df at the
## end (small: ~2000 draws x ~10 params).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_scalar_draws.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = file.path(Sys.getenv("SCRATCH"),
                                          "standard_model_2/fits"))
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCALAR_PARS <- c(
  "sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
  "rho_xi_s", "rho_zeta_s",
  "pi_alpha", "delta", "s", "sigma_s", "sigma_lambda",
  "mu_r", "sigma_r"
)

extract_one <- function(tag) {
  in_path <- file.path(FITS_DIR, paste0(tag, ".rds"))
  if (!file.exists(in_path)) { cat("MISSING:", in_path, "\n"); return(NULL) }
  cat(sprintf("\n== %s ==\n", tag))
  fit <- readRDS(in_path)
  all_vars <- fit$metadata()$stan_variables
  pars_present <- intersect(SCALAR_PARS, all_vars)

  draws_list <- list()
  for (p in pars_present) {
    cat("  ", p, "...\n", sep="")
    dr <- fit$draws(variables = p, format = "draws_df")
    # Keep only the value column (and indexing cols on first iteration)
    if (length(draws_list) == 0) {
      draws_list$.chain     <- dr$.chain
      draws_list$.iteration <- dr$.iteration
      draws_list$.draw      <- dr$.draw
    }
    val_col <- setdiff(names(dr), c(".chain", ".iteration", ".draw"))
    for (vc in val_col) draws_list[[vc]] <- dr[[vc]]
    rm(dr); gc(verbose = FALSE)
  }
  out_df <- as.data.frame(draws_list)
  saveRDS(out_df, file.path(OUT_DIR, paste0(tag, ".draws.rds")))
  cat(sprintf("Wrote .draws.rds (%d draws, %d pars).\n",
              nrow(out_df), length(pars_present)))

  rm(fit, out_df, draws_list); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}
cat("\nAll done.\n")
