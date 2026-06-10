## Extract per-item delta_j medians from a kept cmdstan CSV dir (cut-streaming,
## low RAM -- skips the huge log_lik columns). Used to recover the proc_dp
## item difficulties for the RT-quantile fan figure.
## Usage: STANDARD_MODEL_FITS_DIR=<dir> Rscript sherlock/extract_proc_deltaj.R <tag>
suppressPackageStartupMessages(library(stats))
args <- commandArgs(trailingOnly = TRUE)
TAG  <- args[1]
FITS <- Sys.getenv("STANDARD_MODEL_FITS_DIR", unset = "fits")
CSVDIR <- file.path(FITS, paste0("csvs_", TAG))
files <- sort(list.files(CSVDIR, pattern = "\\.csv$", full.names = TRUE))
stopifnot(length(files) > 0)
hdr  <- system(sprintf("grep -v '^#' %s | head -1", shQuote(files[1])), intern = TRUE)
cols <- strsplit(hdr, ",", fixed = TRUE)[[1]]
dj   <- grep("^delta_j\\.[0-9]+$", cols)   # transformed delta_j (dot notation), NOT delta_j_raw
stopifnot(length(dj) > 0)
fields <- paste(dj, collapse = ",")
allv <- do.call(rbind, lapply(files, function(f) {
  cmd <- sprintf("grep -v '^#' %s | tail -n +2 | cut -d, -f%s", shQuote(f), fields)
  as.matrix(read.csv(pipe(cmd), header = FALSE))
}))
med <- apply(allv, 2, median)
out <- file.path(FITS, "summaries", paste0(TAG, "_psi.csv"))
write.csv(data.frame(jj = seq_along(med), delta_j = med), out, row.names = FALSE)
cat(sprintf("wrote %s (%d delta_j, %d draws)\n", out, length(med), nrow(allv)))
