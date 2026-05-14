## Slim LOO extractor for very-large posteriors (6000+ draws).
## Thins log_lik to N_THIN draws before computing PSIS-LOO so peak
## memory stays bounded. PSIS-LOO is robust to draw thinning as long as
## effective sample size remains adequate (which 2000 well-mixed draws
## generally do).
##
## Usage: Rscript sherlock/extract_loo_thinned.R <tag1> [<tag2> ...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_loo_thinned.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({
  library(posterior)
  library(loo)
})

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = file.path(Sys.getenv("SCRATCH"),
                                          "standard_model_2/fits"))
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_THIN <- 2000L

extract_one <- function(tag) {
  in_path  <- file.path(FITS_DIR, paste0(tag, ".rds"))
  if (!file.exists(in_path)) { cat("MISSING:", in_path, "\n"); return(NULL) }
  cat(sprintf("\n== %s ==\n", tag))
  cat(sprintf("Reading %s ...\n", in_path))
  fit <- readRDS(in_path)
  stopifnot(inherits(fit, "CmdStanMCMC"))

  # Get log_lik as a [iters x chains x N] array, then thin and flatten.
  cat("Extracting log_lik draws (this is the memory-heavy step)...\n")
  ll_full <- fit$draws("log_lik", format = "draws_matrix")
  D <- nrow(ll_full)
  N <- ncol(ll_full)
  cat(sprintf("  log_lik shape: %d draws x %d obs (%.2f GB)\n",
              D, N, prod(dim(ll_full)) * 8 / 1e9))

  # Thin to N_THIN draws by even spacing
  if (D > N_THIN) {
    keep <- as.integer(seq(1, D, length.out = N_THIN))
    ll <- ll_full[keep, , drop = FALSE]
    cat(sprintf("Thinned to %d draws (%.2f GB)\n",
                nrow(ll), prod(dim(ll)) * 8 / 1e9))
    rm(ll_full); gc(verbose = FALSE)
  } else {
    ll <- ll_full
    rm(ll_full); gc(verbose = FALSE)
  }

  # Compute relative-eff approximately (1 chain assumption since we thinned;
  # PSIS-LOO is robust to misspecified r_eff).
  cat("Computing PSIS-LOO...\n")
  r_eff <- rep(1.0, ncol(ll))
  loo_obj <- loo(ll, r_eff = r_eff, cores = 1)
  saveRDS(loo_obj, file.path(OUT_DIR, paste0(tag, ".loo.rds")))
  cat(sprintf("Wrote LOO (elpd = %.1f +- %.1f, pareto_k > 0.7: %d).\n",
              loo_obj$estimates["elpd_loo", "Estimate"],
              loo_obj$estimates["elpd_loo", "SE"],
              sum(loo_obj$diagnostics$pareto_k > 0.7)))

  rm(fit, ll); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}
cat("\nAll done.\n")
