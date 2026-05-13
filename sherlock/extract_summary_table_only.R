## Minimal extraction: cmdstanr's $summary() per scalar variable,
## one at a time. Avoids materializing the full draws matrix that
## OOMs at 16 GB.
##
## NOTE: this version writes ONLY the summary table (no .draws.rds).
## We can pull .draws.rds separately when more memory is available.
##
## Usage: Rscript sherlock/extract_summary_table_only.R <tag1> [<tag2> ...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_summary_table_only.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({ library(cmdstanr); library(dplyr) })

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
  cat(sprintf("Pars present: %s\n", paste(pars_present, collapse=", ")))

  # cmdstanr $summary() per variable; quote = FALSE -> mean+sd+q5+q95+rhat+ess
  summary_rows <- lapply(pars_present, function(p) {
    cat("  ", p, "...\n", sep="")
    s <- fit$summary(variables = p,
                     "mean", "median", "sd",
                     q025 = ~ stats::quantile(.x, 0.025, names = FALSE),
                     q975 = ~ stats::quantile(.x, 0.975, names = FALSE),
                     "ess_bulk", "rhat")
    s
  })
  summary_tbl <- dplyr::bind_rows(summary_rows)
  saveRDS(summary_tbl, file.path(OUT_DIR, paste0(tag, ".summary.rds")))
  cat(sprintf("Wrote .summary.rds (%d rows).\n", nrow(summary_tbl)))

  rm(fit); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}
cat("\nAll done.\n")
