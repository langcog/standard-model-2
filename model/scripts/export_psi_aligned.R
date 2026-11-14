## Write fits/summaries/long_no_freq_slopes_psi.csv keyed to the CURRENT longitudinal bundle by item name.
##
## Two modes.
##
## (A) --refit <tag>_psi.csv   [preferred]
##   Import delta_j medians from a fit on the current bundle (the Sherlock `_2k` refits of 2026-09-03;
##   `cluster/sherlock/recover_from_csvs.R` writes `<tag>_psi.csv` with columns jj, delta_j_median).
##   jj must be the bundle's row index (1..J); the script checks J and the correlation with production.
##   Rscript model/scripts/export_psi_aligned.R --refit fits/summaries/no_freq_slopes_sigmaR_0p44_2k_psi.csv
##
## (B) default: re-key the OLD June fit (fit on the 682-item June-2026 bundle; git 446820b / psi 4d514ed)
##   by item_definition. This is an alignment, not a refit; used as the stopgap on 2026-11-04.
##   Rscript model/scripts/export_psi_aligned.R [bundle_rev] [psi_rev]
##   Old `in` (jj 301) and `inside/in` (jj 303) are one item (ul:in) now; the frequency-matched `in` is used.
##
## Why: the 2026-08-15 bundle rebuild (c97256d) renumbered jj and renamed items ("id:a", "ul:in"), so the
## old jj-keyed file no longer lined up (corr with production ~0 vs 0.99 on the old bundle).
## Optional: --out <path> (default fits/summaries/long_no_freq_slopes_psi.csv).
suppressPackageStartupMessages({library(dplyr); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
take <- function(flag) { i <- match(flag, args); if (is.na(i)) return(NULL)
  v <- args[i + 1]; args <<- args[-c(i, i + 1)]; v }
refit_path <- take("--refit")
out <- take("--out"); if (is.null(out)) out <- "fits/summaries/long_no_freq_slopes_psi.csv"
bundle_rev <- if (length(args) >= 1) args[1] else "446820b"   # June-7 dedup'd 682-item bundle
psi_rev    <- if (length(args) >= 2) args[2] else "4d514ed"   # psi committed 2026-06-10

bundle <- readRDS("fits/long_subset_data.rds")
cur    <- bundle$word_info

if (!is.null(refit_path)) {
  rf <- read_csv(refit_path, show_col_types = FALSE)
  if (!all(c("jj", "delta_j_median") %in% names(rf))) stop("refit csv needs columns jj, delta_j_median")
  if (nrow(rf) != nrow(cur) || !identical(as.integer(rf$jj), seq_len(nrow(cur))))
    stop(sprintf("refit csv has %d rows (jj %s..%s) but the bundle has %d items keyed 1..J; ",
                 nrow(rf), min(rf$jj), max(rf$jj), nrow(cur)),
         "it was fit on a different bundle. Do not import.")
  if (!"delta_j_mean" %in% names(rf)) rf$delta_j_mean <- NA_real_
  aligned <- cur |> select(jj, item, item_definition) |>
    left_join(rf |> select(jj, delta_j_mean, delta_j_median) |> mutate(source_jj = jj), by = "jj")
  note <- c(sprintf("Values: refit on the current bundle, imported from %s (delta_j posterior medians; no means).",
                    basename(refit_path)),
            "jj in the refit csv is the bundle's row index; checked J and the production-rate correlation.")
} else {
  from_git <- function(rev, path, dest) {
    if (system(sprintf("git show %s:%s > %s", rev, path, shQuote(dest))) != 0)
      stop("git show failed for ", rev, ":", path)
    dest
  }
  old_b   <- readRDS(from_git(bundle_rev, "fits/long_subset_data.rds", tempfile(fileext = ".rds")))
  old_psi <- read_csv(from_git(psi_rev, "fits/summaries/long_no_freq_slopes_psi.csv", tempfile(fileext = ".csv")),
                      show_col_types = FALSE)
  stopifnot(nrow(old_b$word_info) == nrow(old_psi))
  old <- old_b$word_info |> select(source_jj = jj, old_item = item) |>
    inner_join(old_psi |> rename(source_jj = jj), by = "source_jj")
  aligned <- cur |> select(jj, item, item_definition) |>
    left_join(old |> select(old_item, source_jj, delta_j_mean, delta_j_median),
              by = c("item_definition" = "old_item"))
  cat(sprintf("old items not used: %s\n",
              paste(setdiff(old$old_item, cur$item_definition), collapse = ", ")))
  note <- c(sprintf("Values: fit on the 682-item June-2026 bundle (git %s), psi as committed in %s.", bundle_rev, psi_rev),
            "Re-keyed by item_definition (source_jj = jj in the old bundle). Alignment only, NOT a refit.")
}
stopifnot(nrow(aligned) == nrow(cur), !anyNA(aligned$delta_j_median), !anyDuplicated(aligned$item))
cat(sprintf("%d / %d current items have a delta_j\n", nrow(aligned), nrow(cur)))

## sanity: delta_j should track production difficulty in the CURRENT bundle
d <- bundle$df |> filter(age >= 22, age <= 28) |>
  group_by(item) |> summarise(rate = mean(produces), .groups = "drop") |>
  inner_join(aligned, by = "item") |> filter(rate > 0, rate < 1)
r <- cor(d$delta_j_median, -qlogis(d$rate))
cat(sprintf("cor(delta_j, -logit production rate, 22-28 mo) = %.3f over %d items\n", r, nrow(d)))
if (r < 0.9) stop("failed the production-rate sanity check (r < 0.9); not written")

write_csv(aligned |> select(jj, item, item_definition, delta_j_mean, delta_j_median, source_jj), out)
writeLines(c(
  "long_no_freq_slopes_psi.csv: EN longitudinal per-word difficulty, keyed to fits/long_subset_data.rds by item.",
  note, "Written by model/scripts/export_psi_aligned.R. Read it via model/R/psi_anchor.R::load_en_psi_by_name()."),
  sub("\\.csv$", ".README.txt", out))
cat("wrote", out, "\n")
