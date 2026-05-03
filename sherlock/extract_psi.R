## Extract per-word psi posterior summary from a longitudinal fit.
##
## For each psi[j] in the fit's draws, computes mean/median/SD/95% CrI
## across MCMC draws, writes a small CSV (~30 KB at J=200). Used by
## the Phase C psi_j ~ log p_j + class regression. Lives separately
## from extract_summaries.R because we don't always want this for
## every fit (psi is large; only m1_time_only / m1 / baseline want it).
##
## Usage on Sherlock:
##   Rscript sherlock/extract_psi.R <tag>
##
## Output:
##   $SCRATCH/standard_model_2/summaries/<tag>_psi.csv

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: extract_psi.R <tag>")
tag <- args[1]

FITS_DIR <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/fits")
OUT_DIR  <- file.path(Sys.getenv("SCRATCH"), "standard_model_2/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

in_path <- file.path(FITS_DIR, paste0(tag, ".rds"))
if (!file.exists(in_path)) stop("Fit not found: ", in_path)

cat(sprintf("Reading %s ...\n", in_path))
fit <- readRDS(in_path)
d <- as_draws_df(fit)

psi_cols <- grep("^psi\\[", names(d), value = TRUE)
if (length(psi_cols) == 0) stop("No psi[] columns in fit draws.")
cat(sprintf("Found %d psi columns.\n", length(psi_cols)))

psi_summary <- tibble(
  jj         = as.integer(sub("psi\\[(\\d+)\\]", "\\1", psi_cols)),
  psi_mean   = sapply(psi_cols, function(p) mean(d[[p]])),
  psi_median = sapply(psi_cols, function(p) median(d[[p]])),
  psi_sd     = sapply(psi_cols, function(p) sd(d[[p]])),
  psi_q025   = sapply(psi_cols, function(p) quantile(d[[p]], 0.025, names = FALSE)),
  psi_q975   = sapply(psi_cols, function(p) quantile(d[[p]], 0.975, names = FALSE))
) %>% arrange(jj)

out_path <- file.path(OUT_DIR, paste0(tag, "_psi.csv"))
write.csv(psi_summary, out_path, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", out_path, nrow(psi_summary)))
