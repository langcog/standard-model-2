## Recover a fit's slim scalar summary + draws (+ delta_j medians) directly
## from the cmdstan CSVs, WITHOUT building a cmdstanr fit object.
##
## Why: cmdstanr::as_cmdstan_fit() parses/indexes the full CSV -- including
## the millions of log_lik columns -- and OOM-kills even a 512 GB machine on
## the big fits (Norwegian, ~210 GB). Here we find the columns we want by
## NAME in the header and stream just those out with `cut`, so RAM stays
## tiny at any fit size. (extract_loo_thinned still needs log_lik and is
## guarded by a size check in run_fit.sh.)
##
## Usage: Rscript sherlock/recover_from_csvs.R <tag> [<csv_dir>]

suppressPackageStartupMessages({ library(dplyr); library(posterior) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript sherlock/recover_from_csvs.R <tag> [<csv_dir>]")
TAG <- args[1]
FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR", unset = "fits")
CSV_DIR  <- if (length(args) >= 2) args[2] else file.path(FITS_DIR, sprintf("csvs_%s", TAG))
SUMM_DIR <- file.path(FITS_DIR, "summaries")
dir.create(SUMM_DIR, recursive = TRUE, showWarnings = FALSE)

SCALAR_PARS <- c("sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
                 "pi_alpha", "delta", "s", "sigma_s",
                 "gamma_in")   # D' input-on-slope coupling (the point of D')

csv_files <- sort(list.files(CSV_DIR, pattern = "\\.csv$", full.names = TRUE))
cat(sprintf("[%s] %d CSV files in %s\n", TAG, length(csv_files), CSV_DIR))
if (length(csv_files) == 0) stop("No CSVs found.")

## Keep only COMPLETE chains: cmdstan appends an "Elapsed Time" footer when a
## chain finishes. A cancelled or OOM'd relaunch can leave empty/partial CSVs
## beside good ones (this is what broke the 2026-08-17 auto-recoveries with
## "no lines available in input"). If several runs share the directory, use
## the most recent complete set (run id = <timestamp>-<hash> in the filename).
complete <- vapply(csv_files, function(f)
  any(grepl("Elapsed Time", system(sprintf("tail -c 4000 %s", shQuote(f)), intern = TRUE))),
  logical(1))
if (any(!complete))
  cat(sprintf("[%s] skipping %d incomplete CSV(s): %s\n", TAG, sum(!complete),
              paste(basename(csv_files[!complete]), collapse = ", ")))
csv_files <- csv_files[complete]
if (length(csv_files) == 0) stop("No complete CSVs found (no 'Elapsed Time' footer).")
run_id <- sub("^.*-(\\d{12})-\\d+-([0-9a-f]+)\\.csv$", "\\1-\\2", basename(csv_files))
if (length(unique(run_id)) > 1) {
  latest <- names(which.max(tapply(file.mtime(csv_files), run_id, max)))
  cat(sprintf("[%s] %d runs share this dir; using latest complete run %s (%d chains)\n",
              TAG, length(unique(run_id)), latest, sum(run_id == latest)))
  csv_files <- csv_files[run_id == latest]
}

## --- header: locate each wanted column by name (no full parse) ---
hdr  <- system(sprintf("grep -v '^#' %s | head -1", shQuote(csv_files[1])), intern = TRUE)
cols <- strsplit(hdr, ",", fixed = TRUE)[[1]]
present     <- SCALAR_PARS[SCALAR_PARS %in% cols]
scalar_idx  <- match(present, cols)
dj_idx      <- grep("^delta_j(\\[|$|\\.)", cols)        # per-word difficulties (optional)
want_idx    <- sort(c(scalar_idx, dj_idx))
want_names  <- cols[want_idx]
cat(sprintf("[%s] scalars: %s | delta_j cols: %d\n", TAG,
            paste(present, collapse = ", "), length(dj_idx)))
if (!length(scalar_idx)) stop("No scalar params found in header.")
fields <- paste(want_idx, collapse = ",")

## --- stream just those columns out of every chain (low RAM) ---
draws_list <- lapply(seq_along(csv_files), function(ch) {
  cmd <- sprintf("grep -v '^#' %s | tail -n +2 | cut -d, -f%s", shQuote(csv_files[ch]), fields)
  con <- pipe(cmd, "r"); on.exit(close(con))
  m <- as.matrix(read.csv(con, header = FALSE))
  colnames(m) <- want_names
  d <- as.data.frame(m); d$.chain <- ch; d$.iteration <- seq_len(nrow(d)); d
})
d <- bind_rows(draws_list)
cat(sprintf("[%s] extracted %d draws across %d chains\n", TAG, nrow(d), length(csv_files)))

## --- scalar summary + draws ---
dd <- as_draws_df(d[, c(".chain", ".iteration", present)])
summ <- summarise_draws(dd, mean, median, sd,
                        q025 = ~quantile(.x, 0.025, names = FALSE),
                        q975 = ~quantile(.x, 0.975, names = FALSE),
                        ess_bulk, rhat)
names(summ)[1] <- "variable"
print(as.data.frame(summ))
saveRDS(as.data.frame(summ), file.path(SUMM_DIR, paste0(TAG, ".summary.rds")))
saveRDS(dd, file.path(SUMM_DIR, paste0(TAG, ".draws.rds")))
cat(sprintf("[%s] wrote .summary.rds + .draws.rds\n", TAG))

## --- delta_j medians (for the exposure figure) ---
if (length(dj_idx)) {
  djcols <- want_names[want_names %in% cols[dj_idx]]
  dj_med <- apply(d[, djcols, drop = FALSE], 2, median)
  write.csv(data.frame(jj = seq_along(dj_med), delta_j_median = unname(dj_med)),
            file.path(SUMM_DIR, paste0(TAG, "_psi.csv")), row.names = FALSE)
  cat(sprintf("[%s] wrote _psi.csv (%d items)\n", TAG, length(dj_med)))
}
cat(sprintf("[%s] DONE\n", TAG))
