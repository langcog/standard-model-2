## Parse the fernald_totlot ("totlot"/TotLot original) item-level WS CDI from the
## three age-band Excel files (data/peekbank/totlot/TL_{18,21,25}m_WS.xlsx) into
## the long format used by data/peekbank/stanford_cdi_items_long.csv, so it can
## join the proc_dp / joint bundle. These 58 kids have longitudinal LWL RT in
## Peekbank (dataset fernald_totlot, lab_subject_id == TL `id`) but were dropped
## before for lack of item-level CDI -- now recovered.
##
## RUN LOCALLY. Output: data/peekbank/totlot_cdi_items_long.csv
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })

TL_DIR <- "data/peekbank/totlot"
FILES  <- c(`18` = "TL_18m_WS.xlsx", `21` = "TL_21m_WS.xlsx", `25` = "TL_25m_WS.xlsx")
ws_map <- read.csv("data/peekbank/cdi_short_code_map_ws.csv")   # short -> item_definition
strip  <- function(x) sub("\\.\\.\\.[0-9]+$", "", x)            # undo name_repair="unique" suffix

long <- bind_rows(lapply(names(FILES), function(band) {
  d  <- read_excel(file.path(TL_DIR, FILES[band]), .name_repair = "unique")
  sc <- strip(names(d))                                          # short-code per column
  id <- as.character(d[[which(sc == "id")[1]]])
  item_idx <- which(sc %in% ws_map$short)                        # only mapped CDI items
  items <- d[, item_idx, drop = FALSE]; names(items) <- names(d)[item_idx]
  bind_cols(lab_subject_id = id, items) %>%
    pivot_longer(-lab_subject_id, names_to = "col", values_to = "raw") %>%
    mutate(short = strip(col),
           produces = as.integer(!is.na(raw) & raw == 1)) %>%    # WS: 1 = says, blank = no
    group_by(lab_subject_id, short) %>%                          # collapse dup item ("feet")
    summarise(produces = max(produces), .groups = "drop") %>%
    mutate(age = as.integer(band))
}))

long <- long %>%
  inner_join(ws_map %>% select(short, item = item_definition, mapping_status = status), by = "short") %>%
  mutate(study = "totlot", form = "WS", source_file = "TL_{18,21,25}m_WS.xlsx") %>%
  select(lab_subject_id, study, age, form, item, produces, short, mapping_status, source_file)

out <- "data/peekbank/totlot_cdi_items_long.csv"
write.csv(long, out, row.names = FALSE)
cat(sprintf("Wrote %s\n  %d rows | %d kids | %d items | ages %s | admins(kid x age) %d | mean produces %.2f\n",
            out, nrow(long), n_distinct(long$lab_subject_id), n_distinct(long$item),
            paste(sort(unique(long$age)), collapse="/"),
            n_distinct(paste(long$lab_subject_id, long$age)), mean(long$produces)))
