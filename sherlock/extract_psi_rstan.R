## Extract per-item delta_j medians from an rstan stanfit. Same output
## format as extract_psi_slim.R but works on stanfit (S4 class with
## different API than CmdStanMCMC).
##
## Usage: Rscript sherlock/extract_psi_rstan.R <tag1> [<tag2> ...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_psi_rstan.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({ library(rstan) })

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = file.path(Sys.getenv("SCRATCH"),
                                          "standard_model_2/fits"))
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

extract_one <- function(tag) {
  in_path <- file.path(FITS_DIR, paste0(tag, ".rds"))
  if (!file.exists(in_path)) { cat("MISSING:", in_path, "\n"); return(NULL) }
  cat(sprintf("\n== %s ==\n", tag))
  fit <- readRDS(in_path)
  stopifnot(inherits(fit, "stanfit"))

  # rstan summary() lazy-computes per-variable means/medians; doesn't
  # materialize the full draws matrix.
  s <- summary(fit, pars = "delta_j",
               probs = c(0.025, 0.5, 0.975))$summary
  if (is.null(s) || nrow(s) == 0) stop("No delta_j in this fit.")
  out <- data.frame(
    jj         = as.integer(sub("delta_j\\[(\\d+)\\]", "\\1", rownames(s))),
    delta_j_mean   = s[, "mean"],
    delta_j_median = s[, "50%"]
  )
  out <- out[order(out$jj), ]
  out_path <- file.path(OUT_DIR, paste0(tag, "_psi.csv"))
  write.csv(out, out_path, row.names = FALSE)
  cat(sprintf("Wrote %s (%d items)\n", out_path, nrow(out)))

  rm(fit, s); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}
cat("\nAll done.\n")
