## Pull English + Norwegian longitudinal WS item-level data from Wordbank.
##
## Usage:   Rscript model/scripts/pull_longitudinal.R
## Outputs: model/fits/long_ws_items.rds
##
## Filters to children with >=2 admins per (language, form) combination and
## keeps only items that appear on our CHILDES frequency table (so we can
## model them with the same prob scale as the cross-sectional subset).

suppressPackageStartupMessages({
  library(wordbankr)
  library(dplyr)
  library(tidyr)
})

source("model/R/config.R")
source("model/R/helpers.R")

LANGUAGES <- c("English (American)", "Norwegian")
FORM      <- "WS"

all_d <- list()
for (lang in LANGUAGES) {
  message(sprintf("Fetching %s %s item data from Wordbank...", lang, FORM))
  d <- get_instrument_data(language = lang, form = FORM,
                           administration_info = TRUE,
                           item_info = TRUE)
  message(sprintf("  got %d item rows", nrow(d)))
  all_d[[lang]] <- d
}

# Identify longitudinal children per (language) in the admin metadata
# Each row of `d` is (admin x item); collapse to admins first.
pick_longitudinal <- function(d, lang) {
  d %>%
    distinct(child_id, age) %>%
    count(child_id, name = "n_admins") %>%
    filter(n_admins >= 2)
}

assemble <- function(d, lang) {
  long_ids <- pick_longitudinal(d, lang)$child_id
  d %>%
    filter(child_id %in% long_ids) %>%
    filter(item_kind == "word") %>%
    select(child_id, age, item = item_definition,
           lexical_category, value) %>%
    mutate(language = lang,
           produces = as.integer(value == "produces"))
}

eng <- assemble(all_d[["English (American)"]], "English (American)")
nor <- assemble(all_d[["Norwegian"]], "Norwegian")
d_long <- bind_rows(eng, nor)

# For English, attach CHILDES frequencies from the preprocessed file we
# used for the cross-sectional fits. (Norwegian will need separate freq
# data — save language for later.)
e <- new.env(); load(PATHS$wordbank, envir = e)
prob_tbl <- e$d_wf %>% distinct(item, prob)

d_long <- d_long %>% left_join(prob_tbl, by = "item")

# For Norwegian + items not in English freq table, leave prob NA for now
cat("\nSummary:\n")
cat(sprintf("  English rows: %d\n", nrow(eng)))
cat(sprintf("  Norwegian rows: %d\n", nrow(nor)))
cat(sprintf("  English unique children: %d\n", length(unique(eng$child_id))))
cat(sprintf("  Norwegian unique children: %d\n", length(unique(nor$child_id))))
cat(sprintf("  English with freq: %d / %d items\n",
            sum(d_long$language == "English (American)" & !is.na(d_long$prob)) /
            sum(d_long$language == "English (American)") * 100,
            n_distinct(d_long$item[d_long$language == "English (American)"])))

# Admin count distributions
cat("\n#admins per child (English):\n")
print(eng %>% distinct(child_id, age) %>% count(child_id) %>%
      count(n, name = "n_children"))

cat("\n#admins per child (Norwegian):\n")
print(nor %>% distinct(child_id, age) %>% count(child_id) %>%
      count(n, name = "n_children"))

# Age-span summary
cat("\nAge span per child (months):\n")
span_tbl <- d_long %>%
  distinct(language, child_id, age) %>%
  group_by(language, child_id) %>%
  summarise(n = n(), span = max(age) - min(age), .groups = "drop")
print(span_tbl %>% group_by(language) %>%
      summarise(n_children = n_distinct(child_id),
                median_n = median(n),
                median_span = median(span),
                max_span = max(span)))

saveRDS(d_long, file.path(PATHS$fits_dir, "long_ws_items.rds"))
cat(sprintf("\nSaved model/fits/long_ws_items.rds (%.1f MB)\n",
            file.info(file.path(PATHS$fits_dir, "long_ws_items.rds"))$size / 1e6))
