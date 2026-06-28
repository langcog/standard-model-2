## Parse SLENA (Weisleder/Fernald Spanish, WF2013) wide CDI -> long item-level.
## WG (SLENA_WGCDIS): each item has ITEM/P (produces) + ITEM/U (understands) columns.
## WS (SLENA_CDIS):  plain item columns, production only.
## The Spanish item set is NON-OVERLAPPING with English (its own bank in the bilingual model);
## item codes here are the SLENA Spanish short-codes (canonical Wordbank-Spanish mapping = TODO).
## ParticipantId == peekbank weisleder_stl lab_subject_id (links to the LWL).
## Output: data/intermediates/slena_cdi_items_long.csv
##   subject_id, age, form, item, produces, comprehends
## RUN LOCALLY.
suppressPackageStartupMessages({library(here); library(readxl); library(dplyr); library(tidyr)})
WF <- here("data/raw/WF2013")
item_cols <- function(nm) {
  ## vocabulary items run from after CDIDate up to the VOCAB total; everything past VOCAB
  ## (word-forms, how-use, sentence complexity) is non-vocab -> excluded (MCF).
  start <- which(nm == "CDIDate") + 1
  stop  <- if (any(nm == "VOCAB")) which(nm == "VOCAB")[1] - 1 else length(nm)
  it <- nm[start:stop]
  it[!grepl("VOCAB|VOCPER|complex|combin|Total|gestos|frases", it, ignore.case = TRUE) & !startsWith(it, "...")]
}
b01 <- function(v) as.integer(!is.na(v) & suppressWarnings(as.numeric(v)) == 1)

## ---- WG: /P -> produces, /U -> understands ----
wg <- read_excel(file.path(WF, "SLENA_WGCDIS_toWordbank.xlsx"), .name_repair = "unique_quiet")
ic <- item_cols(names(wg)); pcols <- ic[grepl("/P$", ic)]; ucols <- ic[grepl("/U$", ic)]
prod <- wg %>% transmute(subject_id = as.character(ParticipantId), age = as.numeric(CDIAge), across(all_of(pcols))) %>%
  pivot_longer(all_of(pcols), names_to = "item", values_to = "p") %>% mutate(item = sub("/P$", "", item), produces = b01(p))
comp <- wg %>% transmute(subject_id = as.character(ParticipantId), age = as.numeric(CDIAge), across(all_of(ucols))) %>%
  pivot_longer(all_of(ucols), names_to = "item", values_to = "u") %>% mutate(item = sub("/U$", "", item), u = b01(u))
wg_long <- prod %>% select(subject_id, age, item, produces) %>%
  full_join(comp %>% select(subject_id, age, item, u), by = c("subject_id", "age", "item")) %>%
  transmute(subject_id, age, form = "WG", item,
            produces = coalesce(produces, 0L), comprehends = as.integer(coalesce(u, 0L) == 1 | coalesce(produces, 0L) == 1))

## ---- WS: production only ----
ws <- read_excel(file.path(WF, "SLENA_CDIS_toWordbank.xlsx"), .name_repair = "unique_quiet")
wsc <- item_cols(names(ws))
ws_long <- ws %>% transmute(subject_id = as.character(ParticipantId), age = as.numeric(CDIAge), across(all_of(wsc))) %>%
  pivot_longer(all_of(wsc), names_to = "item", values_to = "p") %>%
  transmute(subject_id, age, form = "WS", item, produces = b01(p), comprehends = NA_integer_)

long <- bind_rows(wg_long, ws_long) %>% filter(!is.na(age), !is.na(produces))
out <- here("data/intermediates/slena_cdi_items_long.csv")
readr::write_csv(long, out)
cat(sprintf("wrote %s  (%d rows)\n", out, nrow(long)))
long %>% group_by(form) %>% summarise(kids = n_distinct(subject_id),
  admins = n_distinct(paste(subject_id, age)), items = n_distinct(item),
  age_lo = min(age), age_hi = max(age), comp = any(!is.na(comprehends)), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
# sanity: comprehends >= produces on WG
chk <- long %>% filter(form == "WG") %>% group_by(subject_id, age) %>% summarise(p = sum(produces), c = sum(comprehends), .groups = "drop")
cat(sprintf("WG admins with comprehends >= produces: %d / %d\n", sum(chk$c >= chk$p), nrow(chk)))
