## Pull English + Norwegian longitudinal item-level data from Wordbank,
## across BOTH WG and WS forms. Combining the two extends coverage from
## ages 8-30 (instead of just WS's 16-30) and gives every child the
## maximal item set seen across their admin history.
##
## Usage:   Rscript model/scripts/pull_longitudinal.R
## Outputs: fits/long_items.rds  (renamed from old long_ws_items.rds)
##          For backwards compatibility a symlink long_ws_items.rds is
##          also written, pointing at the new file.
##
## Filter rule (longitudinal): keep child if they have >= 2 admins on
## the same form (WG or WS). Some children have admins on both forms;
## both contribute, indexed by (child, age, form).

suppressPackageStartupMessages({
  library(wordbankr)
  library(dplyr)
  library(tidyr)
})

source("model/R/config.R")
source("model/R/helpers.R")

LANGUAGES <- c("English (American)", "Norwegian")
FORMS     <- c("WG", "WS")

all_rows <- list()
for (lang in LANGUAGES) {
  for (form in FORMS) {
    message(sprintf("Fetching %s / %s ...", lang, form))
    d <- tryCatch(
      get_instrument_data(language = lang, form = form,
                          administration_info = TRUE,
                          item_info = TRUE),
      error = function(e) {
        message(sprintf("    skipped: %s", conditionMessage(e)))
        NULL
      })
    if (is.null(d) || nrow(d) == 0) next
    message(sprintf("    got %d item rows", nrow(d)))
    all_rows[[length(all_rows) + 1]] <- d %>%
      filter(item_kind == "word") %>%
      mutate(language = lang, form = form,
             produces = as.integer(value == "produces")) %>%
      select(language, form, child_id, age, item = item_definition,
             lexical_category, produces)
  }
}

d_long <- bind_rows(all_rows)
message(sprintf("\nCombined %d item rows across %d children",
                nrow(d_long), n_distinct(paste(d_long$language, d_long$child_id))))

# Identify longitudinal children: >= 2 admins on the same form
admin_counts <- d_long %>%
  distinct(language, form, child_id, age) %>%
  count(language, form, child_id, name = "n_admins")

long_keys <- admin_counts %>%
  filter(n_admins >= 2) %>%
  distinct(language, form, child_id)

# A child is "longitudinal" if they're long on EITHER form. Keep all
# their data across both forms.
long_child_keys <- long_keys %>% distinct(language, child_id)
d_long <- d_long %>% inner_join(long_child_keys, by = c("language", "child_id"))

# Attach English CHILDES p_j (only for English rows). Norwegian rows
# still have prob = NA and get freq attached by prepare_longitudinal_norwegian.R.
e <- new.env(); load(PATHS$wordbank, envir = e)
prob_tbl <- e$d_wf %>% distinct(item, prob)
d_long <- d_long %>% left_join(prob_tbl, by = "item")

# ---- Reporting ----
cat("\n=== Per-language summary ===\n")
print(d_long %>%
        group_by(language) %>%
        summarise(rows = n(),
                  children = n_distinct(child_id),
                  admins = n_distinct(paste(child_id, age, form)),
                  items = n_distinct(item),
                  forms = paste(sort(unique(form)), collapse = "+"),
                  age_range = sprintf("%.0f-%.0f", min(age), max(age)),
                  with_prob = sprintf("%.0f%%", 100 * mean(!is.na(prob))),
                  .groups = "drop"))

cat("\n=== Per-language admin counts per child ===\n")
for (lang in unique(d_long$language)) {
  cat(sprintf("\n%s:\n", lang))
  print(d_long %>% filter(language == lang) %>%
          distinct(child_id, age, form) %>% count(child_id) %>%
          count(n, name = "n_children") %>% arrange(n))
}

# ---- Save ----
out <- file.path(PATHS$fits_dir, "long_items.rds")
saveRDS(d_long, out)
cat(sprintf("\nSaved %s (%.1f MB)\n", out, file.info(out)$size / 1e6))

# (Old long_ws_items.rds dropped now that all scripts read long_items.rds.)
