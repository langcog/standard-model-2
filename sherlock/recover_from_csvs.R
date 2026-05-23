## Recover a fit's slim summary + scalar draws + delta_j medians directly
## from persistent CSV files, bypassing the save_object() step that
## previously OOM'd on disk-full.
##
## Usage: Rscript sherlock/recover_from_csvs.R <tag> [<csv_dir>]
##   tag      - output tag (e.g. long_no_freq_slopes_norwegian)
##   csv_dir  - directory of cmdstan CSV files (default: fits/csvs_<tag>)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript sherlock/recover_from_csvs.R <tag> [<csv_dir>]")
TAG <- args[1]

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR",
                       unset = "fits")
CSV_DIR  <- if (length(args) >= 2) args[2] else
            file.path(FITS_DIR, sprintf("csvs_%s", TAG))
SUMM_DIR <- file.path(FITS_DIR, "summaries")
dir.create(SUMM_DIR, recursive = TRUE, showWarnings = FALSE)

csv_files <- sort(list.files(CSV_DIR, pattern = "\\.csv$", full.names = TRUE))
cat(sprintf("[%s] %d CSV files found in %s\n", TAG, length(csv_files), CSV_DIR))
if (length(csv_files) == 0) stop("No CSVs found.")

fit <- cmdstanr::as_cmdstan_fit(csv_files)

# Scalar params we care about (same set as extract_summary_table_only.R).
SCALAR_PARS <- c("sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
                 "pi_alpha", "delta", "s", "sigma_s")
present <- SCALAR_PARS[SCALAR_PARS %in% fit$metadata()$stan_variables]
cat(sprintf("[%s] scalar pars present: %s\n", TAG, paste(present, collapse = ", ")))

cat(sprintf("[%s] computing scalar summary...\n", TAG))
slim <- fit$summary(variables = present,
                    mean = mean,
                    lo95 = ~ quantile(.x, 0.025, na.rm = TRUE),
                    median = median,
                    hi95 = ~ quantile(.x, 0.975, na.rm = TRUE),
                    ess_bulk = posterior::ess_bulk,
                    rhat = posterior::rhat)
print(slim)
saveRDS(slim, file.path(SUMM_DIR, paste0(TAG, ".summary.rds")))
cat(sprintf("[%s] wrote summary.rds (%d rows)\n", TAG, nrow(slim)))

cat(sprintf("[%s] extracting scalar draws...\n", TAG))
sc_draws <- fit$draws(variables = present, format = "draws_df")
saveRDS(sc_draws, file.path(SUMM_DIR, paste0(TAG, ".draws.rds")))
cat(sprintf("[%s] wrote draws.rds (%d draws)\n", TAG, nrow(sc_draws)))

cat(sprintf("[%s] extracting delta_j medians...\n", TAG))
dj <- fit$draws(variables = "delta_j", format = "draws_df")
dj_mat <- as.data.frame(dj)
dj_cols <- grep("^delta_j", names(dj_mat), value = TRUE)
dj_med <- apply(dj_mat[, dj_cols, drop = FALSE], 2, median)
out_csv <- data.frame(jj = seq_along(dj_med), delta_j_median = unname(dj_med))
write.csv(out_csv, file.path(SUMM_DIR, paste0(TAG, "_psi.csv")), row.names = FALSE)
cat(sprintf("[%s] wrote _psi.csv (%d items)\n", TAG, nrow(out_csv)))

cat(sprintf("[%s] DONE\n", TAG))
