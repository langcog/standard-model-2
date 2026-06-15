## Parse the ELENA (FMW2013 Batch 1, Stanford, 4-digit ids) item-level CDI at
## ~24mo from the wide Wordbank export (data/peekbank/fmw_2013/FernaldCDIstoWordbank.txt)
## into the long format used by stanford_cdi_items_long.csv, so it joins the
## FMW2013 (fmw_2013) datasets. The file holds both batches at 24mo; we take ONLY
## ELENA (4-digit) here -- TLO (20xxx) already comes from the stanford parse at
## 18/24/30mo (more ages). Vocab = cols 13-692 (680 WS items, 0/1 production),
## in canonical WS form order -> mapped by POSITION to cdi_short_code_map_ws.csv's
## item_definition (avoids the spelling/sense-split ambiguity in the file's names).
## Cols 693+ are Word Forms/Complexity sections (not vocabulary) -> dropped.
## RUN LOCALLY.  Output: data/peekbank/elena_cdi_items_long.csv
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

FILE <- "data/peekbank/fmw_2013/FernaldCDIstoWordbank.txt"
ws   <- read.csv("data/peekbank/cdi_short_code_map_ws.csv")$item_definition   # 680, form order
stopifnot(length(ws) == 680)

## data rows: skip the 2-row header (row1 = item names, row2 = labels)
raw <- read.delim(FILE, skip = 1, header = TRUE, check.names = FALSE, colClasses = "character")
id  <- trimws(raw[[1]]); age <- suppressWarnings(as.integer(round(as.numeric(raw[[2]]))))
vocab <- raw[, 13:692]                                  # the 680 vocabulary columns
stopifnot(ncol(vocab) == 680)
names(vocab) <- ws                                      # POSITION -> canonical item_definition

long <- as_tibble(vocab) %>%
  mutate(lab_subject_id = id, age = age) %>%
  filter(grepl("^[0-9]{4}$", lab_subject_id)) %>%        # ELENA = 4-digit ids only
  pivot_longer(-c(lab_subject_id, age), names_to = "item", values_to = "raw") %>%
  transmute(lab_subject_id, study = "elena", age, form = "WS", item,
            produces = as.integer(trimws(raw) == "1"))

cat(sprintf("ELENA CDI: %d kids, ages %s, %d items, %d rows (mean produces %.2f)\n",
            n_distinct(long$lab_subject_id), paste(sort(unique(long$age)), collapse="/"),
            n_distinct(long$item), nrow(long), mean(long$produces, na.rm = TRUE)))
write_csv(long, "data/peekbank/elena_cdi_items_long.csv")
cat("wrote data/peekbank/elena_cdi_items_long.csv\n")
