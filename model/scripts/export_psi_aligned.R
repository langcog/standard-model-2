## Re-key the EN longitudinal fit's per-word difficulty to the CURRENT longitudinal bundle by item name.
##
## Why: fits/summaries/long_no_freq_slopes_psi.csv came from a fit on the June-2026 bundle (682 items,
## jj-keyed, no item column). The 2026-08-15 bundle rebuild (c97256d) renumbered jj and renamed items
## ("id:a", "ul:in"), so the old file no longer lines up (corr with production ~0 vs 0.99 on the old
## bundle). This script recovers the old bundle from git, maps psi to it by position, then to the current
## bundle by item_definition, and writes a psi file WITH an item column. It is an alignment of the existing
## fit, NOT a refit: the values still come from the June data. A refit on the current bundle is the
## real fix (journal/notes/TODO_psi_681_refit.md).
##
## Old bundle's `in` (jj 301) and `inside/in` (jj 303) are one item (ul:in) in the current bundle; the
## frequency-matched `in` is used.
##
## Usage (repo root): Rscript model/scripts/export_psi_aligned.R [bundle_rev=446820b] [psi_rev=4d514ed]
suppressPackageStartupMessages({library(dplyr); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
bundle_rev <- if (length(args) >= 1) args[1] else "446820b"   # June-7 dedup'd 682-item bundle
psi_rev    <- if (length(args) >= 2) args[2] else "4d514ed"   # psi committed 2026-06-10
out        <- "fits/summaries/long_no_freq_slopes_psi.csv"

from_git <- function(rev, path, dest) {
  st <- system(sprintf("git show %s:%s > %s", rev, path, shQuote(dest)))
  if (st != 0) stop("git show failed for ", rev, ":", path)
  dest
}
old_b   <- readRDS(from_git(bundle_rev, "fits/long_subset_data.rds", tempfile(fileext = ".rds")))
old_psi <- read_csv(from_git(psi_rev, out, tempfile(fileext = ".csv")), show_col_types = FALSE)
cur     <- readRDS("fits/long_subset_data.rds")$word_info

old <- old_b$word_info |> select(source_jj = jj, old_item = item) |>
  inner_join(old_psi |> rename(source_jj = jj), by = "source_jj")
stopifnot(nrow(old_b$word_info) == nrow(old_psi))

aligned <- cur |> select(jj, item, item_definition) |>
  left_join(old |> select(old_item, source_jj, delta_j_mean, delta_j_median),
            by = c("item_definition" = "old_item"))
stopifnot(nrow(aligned) == nrow(cur), !anyNA(aligned$delta_j_median),
          !anyDuplicated(aligned$item))
unused <- setdiff(old$old_item, cur$item_definition)
cat(sprintf("aligned %d / %d current items; old items not used: %s\n",
            nrow(aligned), nrow(cur), paste(unused, collapse = ", ")))

## sanity: delta_j should track production difficulty in the CURRENT bundle
d <- readRDS("fits/long_subset_data.rds")$df |> filter(age >= 22, age <= 28) |>
  group_by(item) |> summarise(rate = mean(produces), .groups = "drop") |>
  inner_join(aligned, by = "item") |> filter(rate > 0, rate < 1)
r <- cor(d$delta_j_median, -qlogis(d$rate))
cat(sprintf("cor(delta_j, -logit production rate, 22-28 mo) = %.3f over %d items\n", r, nrow(d)))
if (r < 0.9) stop("alignment failed the production-rate sanity check")

write_csv(aligned |> select(jj, item, item_definition, delta_j_mean, delta_j_median, source_jj), out)
writeLines(c(
  "long_no_freq_slopes_psi.csv: EN longitudinal per-word difficulty, keyed to fits/long_subset_data.rds by item.",
  sprintf("Values: fit on the 682-item June-2026 bundle (git %s), psi as committed in %s.", bundle_rev, psi_rev),
  "Re-keyed by model/scripts/export_psi_aligned.R (source_jj = jj in the old bundle). Not a refit.",
  "Read it through model/R/psi_anchor.R::load_en_psi_by_name()."),
  sub("\\.csv$", ".README.txt", out))
cat("wrote", out, "\n")
