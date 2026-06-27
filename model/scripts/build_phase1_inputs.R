## Phase-1 modeling inputs for the Spanish + sumscore extension.
## (1) Spanish short-code map DRAFT: SLENA lab codes -> Wordbank-Spanish item_definitions by
##     accent-stripped exact match; the rest are flagged needs_review (lab has the short-code key).
##     -> data/intermediates/cdi_short_code_map_es_{wg,ws}.csv  (short, item_definition, status)
## (2) Count-likelihood admin table for the SUMSCORE cohorts (ELENA-WS English, HABLA Spanish):
##     one row per (child, age, form) with k = words produced, n = form size. The count likelihood
##     (Poisson-binomial over the form's items) attaches these to the latent theta.
##     -> data/intermediates/sumscore_count_admins.csv
## RUN LOCALLY (needs wordbankr).
suppressPackageStartupMessages({library(here); library(dplyr); library(readr); library(readxl); library(wordbankr)})
fp <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(trimws(x), to = "ASCII//TRANSLIT")))   # accent-stripped fingerprint

## ---------- (1) Spanish short-code map draft ----------
sl <- read_csv(here("data/intermediates/slena_cdi_items_long.csv"), show_col_types = FALSE) %>% distinct(form, item)
es_sizes <- list()
for (frm in c("WG", "WS")) {
  wb <- get_item_data(language = "Spanish (Mexican)", form = frm) %>% filter(item_kind == "word") %>%
    distinct(item_definition) %>% mutate(fp = fp(item_definition)) %>% group_by(fp) %>% slice(1) %>% ungroup()
  es_sizes[[frm]] <- nrow(wb)
  m <- sl %>% filter(form == frm) %>% mutate(fp = fp(item)) %>%
    left_join(wb, by = "fp") %>%
    transmute(short = item, item_definition, status = ifelse(!is.na(item_definition), "auto_exact", "needs_review")) %>%
    arrange(status, short)
  out <- here(sprintf("data/intermediates/cdi_short_code_map_es_%s.csv", tolower(frm)))
  write_csv(m, out)
  cat(sprintf("(1) %s: %d codes -> auto_exact %d, needs_review %d  [%s]\n",
              frm, nrow(m), sum(m$status == "auto_exact"), sum(m$status == "needs_review"), basename(out)))
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
