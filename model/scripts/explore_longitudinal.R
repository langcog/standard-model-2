## Explore the Wordbook admins feather file to find longitudinal
## English (American) WS and Norwegian WS administrations.
##
## Adapted from the snippet the user provided.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

ADMINS_PATH <- "/Users/mcfrank/Projects/wordbank/wordbank-book/data/_common/admins.feather"

admins <- read_feather(ADMINS_PATH)
cat("Total admins:", nrow(admins), "\n")
cat("Columns:", paste(colnames(admins), collapse = ", "), "\n")
cat("Languages (first 20):", paste(head(sort(unique(admins$language)), 20),
                                    collapse = ", "), "\n")
cat("Forms:", paste(unique(admins$form), collapse = ", "), "\n\n")

# Longitudinal identifier. Some contributors use the same child_id for
# one child across administrations; the langform pair disambiguates when
# a child appears in WS and WG both.
.inst_sep <- "|"
longitudinal_admins <- admins %>%
  mutate(langform = paste(language, form, sep = .inst_sep)) %>%
  group_by(langform, child_id) %>%
  count() %>%
  filter(n > 1)
cat(sprintf("Longitudinal (same langform, child_id) rows: %d\n",
            nrow(longitudinal_admins)))

# English WS + Norwegian WS
n_long_ws <- admins %>%
  filter(child_id %in% longitudinal_admins$child_id,
         language %in% c("Norwegian", "English (American)"),
         form == "WS") %>%
  group_by(child_id, language, dataset_name) %>%
  mutate(n_admins = n()) %>%
  filter(n_admins > 1) %>%
  ungroup()

cat("\n=== Longitudinal WS (English + Norwegian) ===\n")
cat(sprintf("  admins: %d\n", nrow(n_long_ws)))
cat(sprintf("  unique children: %d\n",
            length(unique(paste(n_long_ws$child_id,
                                 n_long_ws$language,
                                 n_long_ws$dataset_name)))))
cat("  by language:\n")
print(n_long_ws %>% distinct(child_id, language, dataset_name) %>%
      count(language, dataset_name))

# English WG
n_long_wg <- admins %>%
  filter(child_id %in% longitudinal_admins$child_id,
         language %in% c("Norwegian", "English (American)"),
         form == "WG") %>%
  group_by(child_id, language, dataset_name) %>%
  mutate(n_admins = n()) %>%
  filter(n_admins > 1) %>%
  ungroup()

cat("\n=== Longitudinal WG (English + Norwegian) ===\n")
cat(sprintf("  admins: %d\n", nrow(n_long_wg)))
cat(sprintf("  unique children: %d\n",
            length(unique(paste(n_long_wg$child_id,
                                 n_long_wg$language,
                                 n_long_wg$dataset_name)))))
print(n_long_wg %>% distinct(child_id, language, dataset_name) %>%
      count(language, dataset_name))

# Age spans
cat("\n=== Age spans (WS longitudinal) ===\n")
span_tbl <- n_long_ws %>%
  group_by(child_id, language, dataset_name) %>%
  summarise(n_admins = n(),
            min_age = min(age), max_age = max(age),
            span = max_age - min_age,
            .groups = "drop")
cat(sprintf("  median #admins: %.0f\n", median(span_tbl$n_admins)))
cat(sprintf("  median span: %.1f mo (range %.1f–%.1f)\n",
            median(span_tbl$span), min(span_tbl$span), max(span_tbl$span)))

# Save the filtered longitudinal admin info so downstream scripts can use it
saveRDS(list(admins_ws = n_long_ws, admins_wg = n_long_wg),
        file.path("/Users/mcfrank/Projects/standard_model_2/fits",
                  "longitudinal_admins.rds"))
cat("\nSaved: fits/longitudinal_admins.rds\n")
