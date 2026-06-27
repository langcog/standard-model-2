## Build the harmonized item-level CDI table: union every dataset's long CDI into one
## clearly-formatted CSV (for hand-checking + eventual return to Wordbank/Peekbank).
## Inputs (per-dataset long CSVs, each produced by a parse_*_cdi.R):
##   data/intermediates/stanford_cdi_items_long.csv   (studies tlo/totlot2/totlot3)
##   data/intermediates/totlot_cdi_items_long.csv     (totlot)
##   data/intermediates/elena_cdi_items_long.csv      (elena)
##   data/raw/seedlings/cdi_items_long.csv           (seedlings; has comprehends)
##   data/raw/babyview/cdi_items_long.csv            (babyview; raw labels -> canonicalized here)
## Output: data/harmonized/cdi_item_level.csv
##   schema: dataset, paper_code, cohort, child_id, age, form, item, produces, comprehends
## SLENA / WF2013 (Spanish) is pending its parser; not yet included.
## RUN LOCALLY.
suppressPackageStartupMessages({library(here); library(dplyr); library(tidyr); library(readr)})

## study (granular cohort) -> dataset_name + paper code
MAP <- tribble(
  ~cohort,     ~dataset,                  ~paper_code,
  "totlot3",   "adams_marchman_2018",     "AM2018",
  "totlot2",   "fernald_marchman_2012",   "FM2012",
  "tlo",       "fmw_2013",                "FMW2013",
  "elena",     "fmw_2013",                "FMW2013",
  "totlot",    "fernald_totlot",          "FPM2006",
  "seedlings", "seedlings",               "SEEDLingS",
  "babyview",  "babyview",                "babyview")
add_meta <- function(d) left_join(d, MAP, by = "cohort")

## --- peekbank Stanford cohorts (item already canonical Wordbank name) ---
stanford <- read_csv(here("data/intermediates/stanford_cdi_items_long.csv"), show_col_types = FALSE) %>%
  transmute(cohort = study, child_id = as.character(lab_subject_id), age, form, item,
            produces = as.integer(produces), comprehends = as.integer(comprehends), item_canonical = TRUE)
totlot <- read_csv(here("data/intermediates/totlot_cdi_items_long.csv"), show_col_types = FALSE) %>%
  transmute(cohort = "totlot", child_id = as.character(lab_subject_id), age, form, item,
            produces = as.integer(produces), comprehends = NA_integer_, item_canonical = TRUE)
elena <- read_csv(here("data/intermediates/elena_cdi_items_long.csv"), show_col_types = FALSE) %>%
  transmute(cohort = "elena", child_id = as.character(lab_subject_id), age, form, item,
            produces = as.integer(produces), comprehends = as.integer(comprehends), item_canonical = TRUE)

## --- seedlings (item canonical; carries comprehends) ---
seed <- read_csv(here("data/raw/seedlings/cdi_items_long.csv"), show_col_types = FALSE) %>%
  transmute(cohort = "seedlings", child_id = as.character(subject_id), age, form, item,
            produces = as.integer(produces), comprehends = as.integer(comprehends), item_canonical = TRUE)

## --- babyview (webCDI item-keys already mapped to canonical item in parse_babyview_cdi.R) ---
## Keep production-vocab items only, matching the other cohorts (which are vocab by construction).
bv <- read_csv(here("data/raw/babyview/cdi_items_long.csv"), show_col_types = FALSE) %>%
  filter(is_production_vocab) %>%
  transmute(cohort = "babyview", child_id = as.character(subject_id), age, form, item,
            produces = as.integer(produces), comprehends = as.integer(comprehends), item_canonical)

harm <- bind_rows(stanford, totlot, elena, seed, bv) %>% add_meta() %>%
  select(dataset, paper_code, cohort, child_id, age, form, item, produces, comprehends, item_canonical) %>%
  filter(!is.na(produces))

dir.create(here("data/harmonized"), showWarnings = FALSE)
out <- here("data/harmonized/cdi_item_level.csv")
write_csv(harm, out)

cat(sprintf("\nwrote %s  (%d rows)\n\n", out, nrow(harm)))
cat("=== per-dataset CHECK ===\n")
harm %>% group_by(paper_code, cohort) %>%
  summarise(kids = n_distinct(child_id),
            admins = n_distinct(paste(child_id, age, form)),
            items = n_distinct(item),
            age_lo = round(min(age),1), age_hi = round(max(age),1),
            forms = paste(sort(unique(form)), collapse="/"),
            comp = if (all(is.na(comprehends))) "" else "yes", .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
