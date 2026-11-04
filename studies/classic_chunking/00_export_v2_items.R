# Per-item table of v2 word difficulty for the CLASSIC comparison.
#   delta198 : posterior-median per-word difficulty (psi_j) from the named per-word fit
#              (figs/longitudinal/psi_freq_regression_per_word.csv; 197 items)
#   delta681 : v2 delta_j for all 681 items, from fits/summaries/long_no_freq_slopes_psi.csv once it
#              is keyed to the bundle by item (model/scripts/export_psi_aligned.R); June-fit values
#   emp_diff : empirical difficulty, -logit(production rate) at 22-28 months (all items)
# NB: before 2026-11-04 fits/summaries/long_no_freq_slopes_psi.csv was jj-keyed to an older 682-item
# bundle and did NOT align with this one (corr with production ~0); it now carries an item column.
# usage: Rscript 00_export_v2_items.R <out_dir>   (run from the repo root)
suppressPackageStartupMessages({library(dplyr); library(readr)})
out <- commandArgs(TRUE)[1]
b <- readRDS("fits/long_subset_data.rds")
rates <- b$df %>%
  filter(age >= 22, age <= 28) %>%
  group_by(jj) %>%
  summarise(rate = mean(produces), n_obs = n(), .groups = "drop")
items <- b$word_info %>%
  left_join(rates, by = "jj") %>%
  mutate(emp_diff = log((1 - rate) / rate))
psi <- read_csv("figs/longitudinal/psi_freq_regression_per_word.csv", show_col_types = FALSE) %>%
  select(item_definition = item, delta198 = psi_median)   # its jj is a private 1..198 index: join by name
psi681 <- read_csv("fits/summaries/long_no_freq_slopes_psi.csv", show_col_types = FALSE)
if (!"item" %in% names(psi681)) stop("psi file has no item column; run model/scripts/export_psi_aligned.R")
items <- items %>% left_join(psi, by = "item_definition") %>%
  left_join(psi681 %>% select(item, delta681 = delta_j_median), by = "item") %>%
  select(jj, item, item_definition, prob, cc, emp_diff, rate, n_obs, delta198, delta681)
write_csv(items, file.path(out, "v2_items_delta.csv"))
cat(nrow(items), "items;", sum(!is.na(items$delta198)), "with delta198;", sum(!is.na(items$delta681)), "with delta681\n")
