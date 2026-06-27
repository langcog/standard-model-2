## Build the single human-checkable CDI item master key across all io-proc
## datasets: one row per CHILDES-prob item, with its short-code, CHILDES
## frequency, lexical class, per-dataset presence, and whether it's in the
## model's 200-item subsample. Lets a human verify the short-code -> item ->
## CHILDES merge and see item coverage/loss at a glance.
## RUN LOCALLY.  Output: data/intermediates/cdi_master_item_key.csv
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

iop     <- readRDS("fits/io_pooled_subset_data.rds")$df
childes <- iop %>% distinct(item, prob, lexical_category) %>% filter(!is.na(prob), prob > 0)
stanford <- read_csv("data/intermediates/stanford_cdi_items_long.csv", show_col_types = FALSE)
totlot   <- read_csv("data/intermediates/totlot_cdi_items_long.csv",   show_col_types = FALSE)
ws_map   <- read.csv("data/intermediates/cdi_short_code_map_ws.csv")     # short -> item_definition
DSCOLS   <- c("AM2018", "FM2012", "FMW2013", "fernald_totlot", "BabyView", "SEEDLingS")

## short-code per item (WS form; first match)
sc <- ws_map %>% transmute(item = item_definition, short) %>% distinct(item, .keep_all = TRUE)

## per-dataset presence (after the CHILDES-prob gate)
pres <- bind_rows(
  stanford %>% transmute(item, ds = recode(study, totlot2 = "FM2012", totlot3 = "AM2018", tlo = "FMW2013")),
  totlot   %>% transmute(item, ds = "fernald_totlot"),
  iop %>% filter(study %in% c("BabyView", "SEEDLingS")) %>% transmute(item, ds = study)) %>%
  filter(item %in% childes$item) %>% distinct(item, ds) %>% mutate(v = 1L) %>%
  pivot_wider(names_from = ds, values_from = v, values_fill = 0L)

chosen <- readRDS("fits/joint_io_proc_subset_data.rds")$word_info$item
master <- childes %>%
  left_join(sc, by = "item") %>%
  mutate(in_model_200 = item %in% chosen) %>%
  left_join(pres, by = "item")
for (cc in DSCOLS) if (is.null(master[[cc]])) master[[cc]] <- 0L
master <- master %>%
  mutate(across(all_of(DSCOLS), ~ replace_na(., 0L)),
         n_datasets = rowSums(across(all_of(DSCOLS)))) %>%
  select(item, short, lexical_category, prob, in_model_200, n_datasets, all_of(DSCOLS)) %>%
  arrange(desc(in_model_200), lexical_category, desc(prob))

write.csv(master, "data/intermediates/cdi_master_item_key.csv", row.names = FALSE)
cat(sprintf("wrote data/intermediates/cdi_master_item_key.csv: %d items (%d in model-200; %d in all 6 datasets)\n",
            nrow(master), sum(master$in_model_200), sum(master$n_datasets == 6)))
cat(sprintf("  model-200 coverage: AM2018/FM2012/FMW2013/fernald_totlot=200; BabyView=%d; SEEDLingS=%d\n",
            sum(master$BabyView[master$in_model_200]), sum(master$SEEDLingS[master$in_model_200])))
cat(sprintf("  items missing a WS short-code (check): %d\n", sum(is.na(master$short))))
