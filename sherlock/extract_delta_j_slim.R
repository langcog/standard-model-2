## Extract per-item delta_j medians and means, one item at a time, to keep
## memory bounded under 16 GB. The naive approach of pulling all delta_j
## draws as a matrix instantiates ~12 GB which OOMs.
##
## Usage: Rscript sherlock/extract_delta_j_slim.R <tag1> [<tag2> ...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_delta_j_slim.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({ library(cmdstanr) })

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
  all_vars <- fit$metadata()$stan_variables
  if (!"delta_j" %in% all_vars) { cat("No delta_j in this fit; skipping.\n"); return(NULL) }

  # cmdstanr $summary handles delta_j[] as an array; one call returns 1 row
  # per index. This is the lightest path -- doesn't materialize draws.
  s <- fit$summary(variables = "delta_j", "mean", "median")
  s$jj <- as.integer(sub("delta_j\\[(\\d+)\\]", "\\1", s$variable))
  s <- s[order(s$jj), c("jj", "mean", "median")]
  names(s) <- c("jj", "delta_j_mean", "delta_j_median")
  out_path <- file.path(OUT_DIR, paste0(tag, "_psi.csv"))
  write.csv(s, out_path, row.names = FALSE)
  cat(sprintf("Wrote %s (%d items)\n", out_path, nrow(s)))

  rm(fit, s); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}
cat("\nAll done.\n")
