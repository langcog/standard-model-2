## Build the hand-checkable wide CDI item crosswalk: one row per canonical item, one
## column per dataset showing the EXACT raw short-code / item-key that dataset used for it,
## plus a review flag for the non-obvious merges. Lets a human verify every cross-dataset
## item merge at a glance (sort by `flag`). Output: data/harmonized/cdi_item_crosswalk.csv
##   item, form, AM2018, FM2012, FMW2013, FPM2006, SEEDLingS, babyview, n_datasets, flag
## Raw-code sources: stanford/totlot long `short`; seedlings Talk_ map; babyview webCDI key.
## RUN LOCALLY.
suppressPackageStartupMessages({library(here); library(dplyr); library(tidyr); library(readr); library(stringr)})

collapse_codes <- function(x) paste(sort(unique(x[!is.na(x)])), collapse = ";")

## stanford cohorts (AM2018/FM2012/FMW2013): item -> raw short, + manual-disambig flag
st <- read_csv(here("data/intermediates/stanford_cdi_items_long.csv"), show_col_types = FALSE) %>%
  mutate(ds = recode(study, totlot3 = "AM2018", totlot2 = "FM2012", tlo = "FMW2013")) %>%
  group_by(ds, item) %>%
  summarise(code = collapse_codes(short), manual = any(mapping_status == "manual_disambig"), .groups = "drop")

## totlot (FPM2006)
tl <- read_csv(here("data/intermediates/totlot_cdi_items_long.csv"), show_col_types = FALSE) %>%
  group_by(item) %>%
  summarise(code = collapse_codes(short), manual = any(mapping_status == "manual_disambig"), .groups = "drop") %>%
  mutate(ds = "FPM2006")

## seedlings: Talk_ short per canonical item (from the seedlings short-code map)
se <- read_csv(here("data/raw/seedlings/cdi_seedlings_short_code_map.csv"), show_col_types = FALSE) %>%
  transmute(ds = "SEEDLingS", item = item_definition, code = short, manual = status == "manual_disambig")

## babyview: webCDI item-key per canonical production-vocab item (+ homonym flag)
bv <- read_csv(here("data/raw/babyview/cdi_items_long.csv"), show_col_types = FALSE) %>%
  filter(is_production_vocab) %>%
  group_by(item) %>%
  summarise(code = collapse_codes(item_key), homonym = !any(item_canonical), .groups = "drop") %>%
  mutate(ds = "babyview", manual = FALSE)

long <- bind_rows(st, tl, se, bv %>% select(ds, item, code, manual))
bv_homonym <- bv %>% filter(homonym) %>% pull(item)

## form lookup (which CDI form(s) the item lives on) from the harmonized table
forms <- read_csv(here("data/harmonized/cdi_item_level.csv"), show_col_types = FALSE) %>%
  filter(item_canonical) %>% distinct(item, form) %>%
  group_by(item) %>% summarise(form = paste(sort(unique(form)), collapse = "/"), .groups = "drop")

DS <- c("AM2018", "FM2012", "FMW2013", "FPM2006", "SEEDLingS", "babyview")
wide <- long %>%
  select(ds, item, code) %>%
  pivot_wider(names_from = ds, values_from = code) %>%
  left_join(forms, by = "item") %>%
  mutate(n_datasets = rowSums(!is.na(across(all_of(DS)))),
         any_manual = item %in% (long %>% filter(manual) %>% pull(item)),
         flag = case_when(
           item %in% bv_homonym                 ~ "babyview_homonym_check",
           !is.na(babyview) & n_datasets == 1    ~ "babyview_only",
           n_datasets == 1                       ~ "single_dataset",
           any_manual                            ~ "manual_disambig",
           TRUE                                  ~ "")) %>%
  select(item, form, all_of(DS), n_datasets, flag) %>%
  arrange(desc(flag != ""), desc(n_datasets), item)

out <- here("data/harmonized/cdi_item_crosswalk.csv")
write_csv(wide, out)
cat(sprintf("wrote %s  (%d canonical items)\n\n", out, nrow(wide)))
cat("=== flag summary (rows to hand-check) ===\n")
print(wide %>% count(flag, sort = TRUE) %>% as.data.frame(), row.names = FALSE)
cat("\n=== babyview homonyms needing disambiguation ===\n")
print(wide %>% filter(flag == "babyview_homonym_check") %>% select(item, form, babyview) %>% as.data.frame(), row.names = FALSE)
