## Parse BabyView wide CDI exports -> long item-level CSV (CDI-only; fast, no transcripts).
## BabyView is collected via webCDI, so the wide columns are webCDI item-KEYS (e.g. "becarefl",
## "bshteeth"), NOT Wordbank definitions. We map them through the webCDI item dictionary
## (data/raw/babyview/webcdi_items_{WG,WS}.csv, pulled from langcog/web-cdi) to get the canonical
## definition + a production-vocab classifier (choices contains "produces"), then match the
## definition to the Wordbank item list for the shared canonical `item` name.
## Mirrors the wide->long item-column logic in prepare_babyview.R (~lines 86-137).
## Output: data/raw/babyview/cdi_items_long.csv
##   subject_id, age, form, item_key, item, definition, produces, is_production_vocab, item_canonical
## RUN LOCALLY.
suppressPackageStartupMessages({library(here); library(dplyr); library(tidyr); library(readr)})
CDI_DIR <- here("data/raw/babyview/data_june_2026")

META_COLS <- c("opt_out","study_name","subject_id","local_lab_id","repeat_num",
  "administration_id","link","completed","completedBackgroundInfo","due_date",
  "last_modified","created_date","completed_date","event_id","age","sex","country",
  "zip_code","birth_order","birth_weight_lb","birth_weight_kg","multi_birth_boolean",
  "multi_birth","sibling_boolean","sibling_count","sibling_data","born_on_due_date",
  "early_or_late","due_date_diff","form_filler","form_filler_other","primary_caregiver",
  "primary_caregiver_other","mother_yob","mother_education","secondary_caregiver",
  "secondary_caregiver_other","father_yob","father_education","annual_income",
  "child_hispanic_latino","child_ethnicity","caregiver_info","caregiver_other",
  "other_languages_boolean","other_languages","language_from","language_days_per_week",
  "language_hours_per_day","ear_infections_boolean","ear_infections","hearing_loss_boolean",
  "hearing_loss","vision_problems_boolean","vision_problems","illnesses_boolean","illnesses",
  "services_boolean","services","worried_boolean","worried","learning_disability_boolean",
  "learning_disability","children_comforted","show_respect","close_bonds","parents_help_learn",
  "play_learning","explore_experiment","do_as_told","read_at_home","teach_alphbet",
  "rhyming_games","read_for_pleasure","child_asks_for_reading","child_self_reads",
  "child_asks_words_say","place_of_residence","primary_caregiver_occupation",
  "primary_caregiver_occupation_description","secondary_caregiver_occupation",
  "secondary_caregiver_occupation_description","kindergarten_since_when",
  "kindergarten_hpd","kindergarten_dpw")

to_long <- function(df, form_label) {
  meta_present <- intersect(META_COLS, colnames(df))
  item_cols <- setdiff(colnames(df), meta_present)
  item_cols <- item_cols[!grepl("^Total|^Word|^Combining|^Complexity|benchmark|^How |^Combination |Percentile|Combining",
                                item_cols)]
  df %>%
    filter(completed, !is.na(age), age != 999, age >= 8, age <= 36) %>%
    select(subject_id, age, all_of(item_cols)) %>%
    mutate(across(-c(subject_id, age), as.character)) %>%
    pivot_longer(-c(subject_id, age), names_to = "item_key", values_to = "raw") %>%
    mutate(produces = as.integer(!is.na(raw) & raw == "produces"),
           # WG records comprehension ("understands" or "produces"); WS is production-only.
           comprehends = if (form_label == "WG")
             as.integer(!is.na(raw) & raw %in% c("understands", "produces")) else NA_integer_,
           form = form_label)
}

## webCDI item dictionary: item-key -> definition + choices (the vocab classifier)
norm <- function(x) x %>% tolower() %>% gsub("[\\.\\(\\)]","_",.,perl=TRUE) %>% gsub("\\s+|/","_",.) %>%
  gsub("__+","_",.) %>% gsub("[^a-z0-9_]","",.) %>% sub("_+$","",.)
## is_production_vocab = a "produces"-type item in a real vocabulary category -- EXCLUDES the
## word_forms_/word_endings_ morphology sections (irregular pasts/plurals: ate, blew, children,
## feet), which are "produces"-type but not part of the production-vocabulary checklist.
dict <- bind_rows(
  read_csv(here("data/raw/babyview/webcdi_items_WG.csv"), show_col_types = FALSE) %>% mutate(form = "WG"),
  read_csv(here("data/raw/babyview/webcdi_items_WS.csv"), show_col_types = FALSE) %>% mutate(form = "WS")) %>%
  transmute(form, item_key = item, definition, category,
            is_production_vocab = grepl("produces", coalesce(choices, "")) &
                                  !grepl("^word_(forms|endings)", coalesce(category, ""))) %>%
  distinct(form, item_key, .keep_all = TRUE)

## Wordbank canonical item list (match webCDI definition -> shared canonical name)
wb <- readRDS(here("fits/long_items.rds")) %>% filter(language == "English (American)") %>%
  distinct(item) %>% mutate(def_norm = norm(item)) %>%
  group_by(def_norm) %>% slice(1) %>% ungroup() %>% rename(item_canon = item)

ws <- read_csv(file.path(CDI_DIR, "babyview-english-ws_items.csv"), show_col_types = FALSE, progress = FALSE)
wg <- read_csv(file.path(CDI_DIR, "babyview-english-wg_items.csv"), show_col_types = FALSE, progress = FALSE)
cdi <- bind_rows(to_long(ws, "WS"), to_long(wg, "WG")) %>% filter(!is.na(produces)) %>%
  left_join(dict, by = c("form", "item_key")) %>%
  mutate(is_production_vocab = coalesce(is_production_vocab, FALSE),
         def_norm = norm(definition)) %>%
  left_join(wb, by = "def_norm") %>%
  mutate(item = coalesce(item_canon, definition, item_key),
         item_canonical = !is.na(item_canon))

## hand-checked overrides for the ~14 items where webCDI labels diverge from the Wordbank
## canonical strings used by the stanford cohorts (spelling/plural/homonym). See
## data/raw/babyview/babyview_item_overrides.csv (MCF-reviewable). Applied by webCDI item_key.
ov <- read_csv(here("data/raw/babyview/babyview_item_overrides.csv"), show_col_types = FALSE)
cdi <- cdi %>% left_join(ov, by = "item_key") %>%
  mutate(item = coalesce(canonical_item, item),
         item_canonical = item_canonical | !is.na(canonical_item)) %>%
  select(subject_id, age, form, item_key, item, definition, category, produces, comprehends, is_production_vocab, item_canonical)

out <- here("data/raw/babyview/cdi_items_long.csv")
write_csv(cdi, out)
v <- cdi %>% filter(is_production_vocab)
cat(sprintf("wrote %s\n  rows=%d  admins=%d  kids=%d\n", out, nrow(cdi),
            n_distinct(paste(cdi$subject_id, cdi$age, cdi$form)), n_distinct(cdi$subject_id)))
cat(sprintf("  production-vocab items: %d (canonical-matched: %d) | non-vocab item rows kept too\n",
            n_distinct(v$item_key), n_distinct(v$item_key[v$item_canonical])))
v %>% group_by(form) %>% summarise(prod_items = n_distinct(item_key),
  canonical = n_distinct(item_key[item_canonical]), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
