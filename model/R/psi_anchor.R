## Name-keyed access to the EN longitudinal fit's per-word difficulty (delta_j).
##
## fits/summaries/long_no_freq_slopes_psi.csv must carry an `item` column that matches
## fits/long_subset_data.rds$word_info$item. A bare jj-keyed psi silently misaligns when the bundle
## is rebuilt (June 2026: 682 items; Aug 2026 rebuild: 681 items with new "id:"/"ul:" item names), so
## every consumer reads it through this function instead of joining on jj.
##
## Returns tibble(item, item_definition, delta_j_median) in word_info order. Callers key it to their
## own item strings through `item_definition` (the plain CDI word, e.g. "a lot"), not `item`.
load_en_psi_by_name <- function(fits_dir) {
  suppressPackageStartupMessages({library(dplyr); library(readr)})
  b   <- readRDS(file.path(fits_dir, "long_subset_data.rds"))
  psi <- read_csv(file.path(fits_dir, "summaries", "long_no_freq_slopes_psi.csv"),
                  show_col_types = FALSE)
  wi  <- b$word_info
  if (!"item" %in% names(psi))
    stop("long_no_freq_slopes_psi.csv has no `item` column; a jj join would misalign delta_j ",
         "with words. Regenerate it: Rscript model/scripts/export_psi_aligned.R")
  if (nrow(psi) != nrow(wi) || !setequal(psi$item, wi$item))
    stop(sprintf("psi (%d items) does not match long_subset_data word_info (%d items) by item name; ",
                 nrow(psi), nrow(wi)),
         "refit or re-export psi for the current bundle (model/scripts/export_psi_aligned.R).")
  wi |> select(item, item_definition) |>
    inner_join(psi |> select(item, delta_j_median), by = "item")
}

## Fail loudly when an anchor join matched (almost) nothing, e.g. after an item-naming change.
check_anchor_coverage <- function(anchored, label = "delta_j anchor", min_frac = 0.25) {
  frac <- mean(anchored)
  if (sum(anchored) == 0)
    stop(label, ": 0 items anchored; the item-name join failed.")
  if (frac < min_frac)
    warning(sprintf("%s: only %.0f%% of items anchored (%d / %d)", label, 100 * frac,
                    sum(anchored), length(anchored)))
  invisible(frac)
}
