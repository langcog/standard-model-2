## Phase-1 modeling inputs for the Spanish + sumscore extension.
## (1) Spanish short-code map DRAFT: SLENA lab codes -> Wordbank-Spanish item_definitions by
##     accent-stripped exact match; the rest are flagged needs_review (lab has the short-code key).
##     -> data/intermediates/cdi_short_code_map_es_{wg,ws}.csv  (short, item_definition, status)
## (2) Count-likelihood admin table for the SUMSCORE cohorts (ELENA-WS English, HABLA Spanish):
##     one row per (child, age, form) with k = words produced, n = form size. The count likelihood
##     (Poisson-binomial over the form's items) attaches these to the latent theta.
##     -> data/intermediates/sumscore_count_admins.csv
## RUN LOCALLY.
suppressPackageStartupMessages({library(here); library(dplyr); library(readr); library(readxl)})
fp <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(trimws(x), to = "ASCII//TRANSLIT")))   # accent-stripped fingerprint

## ---------- (1) Spanish short-code map: SLENA codes -> Wordbank-Spanish definitions ----------
## The Wordbank raw instrument CSVs (data/raw/WF2013/wordbank_spanish_{WG,WS}.csv) are the
## authoritative dictionary: their `item` column IS the SLENA lab short-code, `definition` the
## canonical item. So this maps the codes directly -- the residual `needs_review` rows are mostly
## the non-vocab word-form/complexity codes (leave item_definition blank to exclude).
sl <- read_csv(here("data/intermediates/slena_cdi_items_long.csv"), show_col_types = FALSE) %>% distinct(form, item)
dict_of <- function(frm) read_csv(here(sprintf("data/raw/WF2013/wordbank_spanish_%s.csv", frm)), show_col_types = FALSE) %>%
  filter(type == "word") %>% transmute(key = tolower(item), item_definition = definition)
DWG <- dict_of("WG"); DWS <- dict_of("WS")
DALL <- bind_rows(DWG, DWS) %>% distinct(key, .keep_all = TRUE)         # union -> cross-form coverage
es_sizes <- list(WG = nrow(DWG), WS = nrow(DWS))
for (frm in c("WG", "WS")) {
  m <- sl %>% filter(form == frm) %>% mutate(key = tolower(item)) %>%
    left_join(DALL, by = "key") %>%
    mutate(status = ifelse(!is.na(item_definition), "wordbank", "needs_review"),
           fuzzy_candidate = NA_character_, fuzzy_dist = NA_integer_)
  nr <- which(m$status == "needs_review")               # residual: nearest dict code by Levenshtein
  if (length(nr) > 0) {
    D  <- adist(m$key[nr], DALL$key); bi <- max.col(-D, ties.method = "first")
    m$fuzzy_candidate[nr] <- DALL$item_definition[bi]
    m$fuzzy_dist[nr]      <- D[cbind(seq_along(nr), bi)]
  }
  m <- m %>% transmute(short = item, item_definition, status, fuzzy_candidate, fuzzy_dist) %>%
    arrange(status, fuzzy_dist, short)
  out <- here(sprintf("data/intermediates/cdi_short_code_map_es_%s.csv", tolower(frm)))
  write_csv(m, out)
  cat(sprintf("(1) %s: %d codes -> wordbank %d, needs_review %d  [%s]\n",
              frm, nrow(m), sum(m$status == "wordbank"), sum(m$status == "needs_review"), basename(out)))
}
N_ES_WG <- es_sizes[["WG"]]; N_ES_WS <- es_sizes[["WS"]]; N_EN_WS <- 680

## ---------- (2) count-likelihood admin table ----------
elena <- read_excel(here("data/raw/FMW2013/elena/ELENA_WS_SumScores.xlsx")) %>%
  transmute(cohort = "ELENA", paper_code = "FMW2013", language = "English",
            child_id = as.character(ParticipantId), age = as.numeric(CDIAge), form = "WS",
            k = as.numeric(VOCAB), n = N_EN_WS)
hb <- read_csv(here("data/raw/Bang2025/Habla1.0_LENALWLCDISumScores.csv"), show_col_types = FALSE)
habla <- bind_rows(
  hb %>% transmute(child_id = as.character(ID), age = as.numeric(CDIAge18m),   form = "WG", k = as.numeric(WordsProd18m), n = N_ES_WG),
  hb %>% transmute(child_id = as.character(ID), age = as.numeric(CDIWSAge),     form = "WS", k = as.numeric(CDIVocPost21), n = N_ES_WS),
  hb %>% transmute(child_id = as.character(ID), age = as.numeric(CDIAgePost25), form = "WS", k = as.numeric(CDIVocPost25), n = N_ES_WS)) %>%
  mutate(cohort = "HABLA", paper_code = "Bang2025", language = "Spanish")
adm <- bind_rows(elena, habla) %>%
  filter(!is.na(k), !is.na(age), k >= 0, k <= n) %>%
  select(cohort, paper_code, language, child_id, age, form, k, n) %>% arrange(cohort, child_id, age)
write_csv(adm, here("data/intermediates/sumscore_count_admins.csv"))
cat(sprintf("\n(2) sumscore_count_admins.csv: %d admins\n", nrow(adm)))
adm %>% group_by(cohort, language, form) %>%
  summarise(kids = n_distinct(child_id), admins = n(), med_k = median(k), n_form = first(n),
            age_lo = min(age), age_hi = max(age), .groups = "drop") %>% as.data.frame() %>% print(row.names = FALSE)
