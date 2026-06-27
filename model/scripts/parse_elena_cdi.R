## Parse ELENA (FMW2013 Batch 1, Stanford) WG CDI from the CORRECTED wide export
## (data/raw/FMW2013/elena/ELENA_WG_items.xlsx). One WG administration per kid at
## 16-18mo. The ids (4943,6117,...) now reconcile with the ELENA LENA SubjectID1
## (all 26 overlap), unblocking the study. Vocabulary items are short-coded columns
## matching cdi_short_code_map_wg.csv; WG production = value 2 ("understands and
## says") (1 = understands only, blank = neither). Emits the long schema used by
## prepare_io_dataset.R's fmw2013 hook (study == "elena").
## NOTE: supersedes the earlier txt-based WS parse, whose ids (6143..) did NOT
## reconcile with the RT/LENA -- that export was on a different id scheme.
## RUN LOCALLY.  Output: data/intermediates/elena_cdi_items_long.csv
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(readxl) })

e  <- read_excel("data/raw/FMW2013/elena/ELENA_WG_items.xlsx")
wg <- read_csv("data/intermediates/cdi_short_code_map_wg.csv", show_col_types = FALSE)
itemcols <- intersect(names(e), wg$short)                  # vocab short-codes (395)
map <- setNames(wg$item_definition, wg$short)

long <- e %>%
  select(id, CDIAge, all_of(itemcols)) %>%
  pivot_longer(all_of(itemcols), names_to = "short", values_to = "raw") %>%
  transmute(lab_subject_id = as.character(id), study = "elena",
            age = as.integer(CDIAge), form = "WG", item = unname(map[short]),
            produces = as.integer(!is.na(raw) & suppressWarnings(as.numeric(raw)) == 2),
            # WG: 1 = understands, 2 = understands+says -> comprehension is either.
            comprehends = as.integer(!is.na(raw) & suppressWarnings(as.numeric(raw)) %in% c(1, 2)))

pk <- long %>% group_by(lab_subject_id) %>% summarise(n_prod = sum(produces), .groups = "drop")
cat(sprintf("ELENA WG CDI: %d kids, ages %s, %d items, %d rows\n",
            n_distinct(long$lab_subject_id), paste(sort(unique(long$age)), collapse = "/"),
            n_distinct(long$item), nrow(long)))
cat(sprintf("  words produced / kid: median=%.0f  range=%d-%d  (16-18mo WG)\n",
            median(pk$n_prod), min(pk$n_prod), max(pk$n_prod)))
write_csv(long, "data/intermediates/elena_cdi_items_long.csv")
cat("wrote data/intermediates/elena_cdi_items_long.csv\n")
