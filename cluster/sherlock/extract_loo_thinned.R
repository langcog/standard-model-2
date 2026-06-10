## Slim LOO extractor for very-large posteriors (6000+ draws, 2M+ obs).
##
## Thins log_lik to N_THIN draws before computing PSIS-LOO, then computes
## PSIS-LOO in column chunks (each chunk is at most CHUNK_OBS observations).
## PSIS-LOO computations are independent per observation, so chunking is
## exact (no approximation) -- aggregation is just elementwise.
##
## Why chunked: a single 35 GB log_lik on a 128 GB machine plus the loo
## package's internal lw_list copies has previously segfaulted inside
## psis_apply on the NO bundle (2.18M obs x 2000 draws). Chunking caps
## peak memory at chunk_size * draws * 8 bytes + the underlying ll matrix.
##
## Usage: Rscript sherlock/extract_loo_thinned.R <tag1> [<tag2> ...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: extract_loo_thinned.R <tag1> [<tag2> ...]")

suppressPackageStartupMessages({
  library(posterior)
  library(loo)
  library(cmdstanr)
})

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = file.path(Sys.getenv("SCRATCH"),
                                          "standard_model_2/fits"))
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_THIN    <- 2000L      # max draws kept after thinning
CHUNK_OBS <- 250000L    # observations per chunk (~4 GB at 2000 draws)

# ---------- chunked PSIS-LOO --------------------------------------------
# Compute PSIS-LOO on column chunks, return a slim list mimicking a loo()
# return value: elpd_loo / p_loo / looic estimates, pointwise table, and
# Pareto-k diagnostics. Skips the full psis_object (which we don't save
# anyway, and which is the largest object in a loo() return).
chunked_loo <- function(ll, chunk_obs = CHUNK_OBS) {
  N <- ncol(ll)
  S <- nrow(ll)

  elpd_pw  <- numeric(N)
  p_pw     <- numeric(N)
  mcse_pw  <- numeric(N)
  k_pw     <- numeric(N)
  n_eff_pw <- numeric(N)

  n_chunks <- ceiling(N / chunk_obs)
  for (ci in seq_len(n_chunks)) {
    a <- (ci - 1L) * chunk_obs + 1L
    b <- min(ci * chunk_obs, N)
    cat(sprintf("  chunk %d/%d: obs %d..%d (%d obs)\n",
                ci, n_chunks, a, b, b - a + 1L))

    # drop=FALSE keeps matrix shape even at the tail chunk
    ll_chunk <- ll[, a:b, drop = FALSE]
    r_eff_chunk <- rep(1.0, ncol(ll_chunk))

    loo_chunk <- loo(ll_chunk, r_eff = r_eff_chunk, cores = 1L,
                     save_psis = FALSE)

    pw <- loo_chunk$pointwise
    elpd_pw[a:b]  <- pw[, "elpd_loo"]
    p_pw[a:b]     <- pw[, "p_loo"]
    mcse_pw[a:b]  <- pw[, "mcse_elpd_loo"]
    k_pw[a:b]     <- loo_chunk$diagnostics$pareto_k
    n_eff_pw[a:b] <- loo_chunk$diagnostics$n_eff

    rm(ll_chunk, loo_chunk, r_eff_chunk)
    gc(verbose = FALSE)
  }

  # Aggregate point estimates + SEs across observations.
  # SE of a sum of independent obs-level contributions is sqrt(N * var).
  elpd_loo <- sum(elpd_pw)
  elpd_se  <- sqrt(N * stats::var(elpd_pw))
  p_loo    <- sum(p_pw)
  p_se     <- sqrt(N * stats::var(p_pw))
  looic    <- -2 * elpd_loo
  looic_se <- 2 * elpd_se

  estimates <- matrix(
    c(elpd_loo, elpd_se,
      p_loo,    p_se,
      looic,    looic_se),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("elpd_loo", "p_loo", "looic"),
                    c("Estimate", "SE"))
  )

  pointwise <- cbind(
    elpd_loo      = elpd_pw,
    mcse_elpd_loo = mcse_pw,
    p_loo         = p_pw,
    looic         = -2 * elpd_pw
  )

  out <- list(
    estimates   = estimates,
    pointwise   = pointwise,
    diagnostics = list(pareto_k = k_pw, n_eff = n_eff_pw),
    n_eff       = n_eff_pw,
    elpd_loo    = elpd_loo,
    se_elpd_loo = elpd_se,
    p_loo       = p_loo,
    looic       = looic,
    N           = N,
    S           = S,
    method      = "psis_chunked"
  )
  class(out) <- c("psis_loo", "loo")
  out
}

# ---------- per-tag driver ----------------------------------------------
extract_one <- function(tag) {
  in_path  <- file.path(FITS_DIR, paste0(tag, ".rds"))
  csv_dir  <- file.path(FITS_DIR, paste0("csvs_", tag))
  cat(sprintf("\n== %s ==\n", tag))
  if (file.exists(in_path)) {
    cat(sprintf("Reading %s ...\n", in_path))
    fit <- readRDS(in_path)
  } else {
    # No self-contained RDS (e.g. STAN_SKIP_SAVE_OBJECT runs, or a recovery
    # after a failed save): reconstruct the fit straight from the CSVs.
    csvs <- if (dir.exists(csv_dir))
      list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE) else character(0)
    if (length(csvs) == 0) {
      cat("MISSING:", in_path, "and no CSVs in", csv_dir, "\n"); return(NULL)
    }
    cat(sprintf("RDS absent; reconstructing from %d CSVs in %s ...\n",
                length(csvs), csv_dir))
    fit <- cmdstanr::as_cmdstan_fit(csvs)
  }
  stopifnot(inherits(fit, "CmdStanMCMC"))

  cat("Extracting log_lik draws (memory-heavy step)...\n")
  ll_full <- fit$draws("log_lik", format = "draws_matrix")
  D <- nrow(ll_full)
  N <- ncol(ll_full)
  cat(sprintf("  log_lik shape: %d draws x %d obs (%.2f GB)\n",
              D, N, prod(dim(ll_full)) * 8 / 1e9))

  # Thin draws by even spacing
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

  # Free the fit object before doing the heavy PSIS work; we don't need it
  # anymore and it's holding additional memory.
  rm(fit); gc(verbose = FALSE)

  cat(sprintf("Computing PSIS-LOO in chunks of %d obs ...\n", CHUNK_OBS))
  loo_obj <- chunked_loo(ll, chunk_obs = CHUNK_OBS)

  saveRDS(loo_obj, file.path(OUT_DIR, paste0(tag, ".loo.rds")))
  cat(sprintf("Wrote LOO (elpd = %.1f +- %.1f, pareto_k > 0.7: %d / %d obs).\n",
              loo_obj$estimates["elpd_loo", "Estimate"],
              loo_obj$estimates["elpd_loo", "SE"],
              sum(loo_obj$diagnostics$pareto_k > 0.7),
              length(loo_obj$diagnostics$pareto_k)))

  rm(ll, loo_obj); gc(verbose = FALSE)
}

for (tag in args) {
  tryCatch(extract_one(tag),
           error = function(e) cat(sprintf("ERROR for %s: %s\n", tag, conditionMessage(e))))
}
cat("\nAll done.\n")
