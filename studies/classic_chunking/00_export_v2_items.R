# Per-item table of v2 word difficulty for the CLASSIC comparison.
#   delta198 : posterior-median per-word difficulty (psi_j) from the named per-word fit
#              (figs/longitudinal/psi_freq_regression_per_word.csv; 197 items)
#   emp_diff : empirical difficulty, -logit(production rate) at 22-28 months (all items)
# NB: fits/summaries/long_no_freq_slopes_psi.csv (682 rows) does NOT align with this
# bundle's 681 items (corr with production ~0); do not use it for item-level work.
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
items <- items %>% left_join(psi, by = "item_definition") %>%
  select(jj, item, item_definition, prob, cc, emp_diff, rate, n_obs, delta198)
write_csv(items, file.path(out, "v2_items_delta.csv"))
cat(nrow(items), "items;", sum(!is.na(items$delta198)), "with delta198\n")
