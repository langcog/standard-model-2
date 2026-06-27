## Build the harmonized INPUT and LWL-SESSION tables: one clearly-formatted CSV each,
## merged across datasets (for hand-checking + eventual return to Peekbank). Trial-level
## LWL is a separate (heavier) build; this is the per-session rollup the model uses.
## Outputs:
##   data/harmonized/input_level.csv      dataset, paper_code, cohort, child_id, age, instrument, input_rate, log_input, units
##   data/harmonized/lwl_session_level.csv dataset, paper_code, child_id, age, log_rt, rt_ms, accuracy, n_trials_rt
## RUN LOCALLY.
suppressPackageStartupMessages({library(here); library(dplyr); library(tidyr); library(readr)})

PC <- c(adams_marchman_2018 = "AM2018", fernald_marchman_2012 = "FM2012", fmw_2013 = "FMW2013",
        fernald_totlot = "FPM2006", seedlings_zhu = "SEEDLingS")

## ---------- INPUT ----------
## AM2018 (TL3): LENA AWC/hr at 16 + 18 mo
am <- read_csv(here("data/raw/AM2018/lena_am2018_fmw2013.csv"), show_col_types = FALSE) %>%
  filter(Study == "TL3") %>%
  transmute(child_id = as.character(SubjectID1), a16 = AGE16M, r16 = AWCHr16M, a18 = AGE18M, r18 = AWCHr18M) %>%
  pivot_longer(-child_id, names_to = c(".value", "tp"), names_pattern = "([ar])(16|18)") %>%
  transmute(paper_code = "AM2018", cohort = "totlot3", child_id, age = a, input_rate = r)
## FMW2013 (TLO + ELENA): LENA AWC/hr at 18 + 24 mo
fmw <- read_csv(here("data/raw/FMW2013/TLOELENA_LENA_1824.csv"), show_col_types = FALSE) %>%
  filter(Study %in% c("TLO", "ELENA")) %>%
  transmute(cohort = ifelse(Study == "ELENA", "elena", "tlo"), child_id = as.character(SubjectID1),
            a18 = AGE18M, r18 = AWCHr18M, a24 = AGE24M, r24 = AWCHr24M) %>%
  pivot_longer(c(a18, r18, a24, r24), names_to = c(".value", "tp"), names_pattern = "([ar])(18|24)") %>%
  transmute(paper_code = "FMW2013", cohort, child_id, age = a, input_rate = r)
## SEEDLingS: monthly LENA AWC/hr
seed <- read_csv(here("data/raw/seedlings/lena_data.csv"), show_col_types = FALSE) %>%
  transmute(paper_code = "SEEDLingS", cohort = "seedlings", child_id = as.character(subj),
            age = month, input_rate = awc_perhr)
lena <- bind_rows(am, fmw, seed) %>% filter(!is.na(input_rate), input_rate > 0, !is.na(age)) %>%
  mutate(instrument = "LENA", log_input = log(input_rate), units = "adult words / hour")
## BabyView: per-video head-cam adult-token log-rate (from the prepared bundle's video table)
bv <- readRDS(here("fits/babyview_subset_data.rds"))$videos %>%
  transmute(paper_code = "babyview", cohort = "babyview", child_id = as.character(subject_id),
            age = age_mo, log_input = log_r_obs, input_rate = exp(log_r_obs),
            instrument = "head-cam transcript", units = "adult tokens (study-specific log-rate)")
input <- bind_rows(lena, bv) %>%
  mutate(dataset = recode(paper_code, AM2018 = "adams_marchman_2018", FM2012 = "fernald_marchman_2012",
                          FMW2013 = "fmw_2013", FPM2006 = "fernald_totlot", SEEDLingS = "seedlings",
                          babyview = "babyview")) %>%
  select(dataset, paper_code, cohort, child_id, age, instrument, input_rate, log_input, units) %>%
  arrange(paper_code, child_id, age)
write_csv(input, here("data/harmonized/input_level.csv"))

## ---------- LWL (session level) ----------
OURS <- c("adams_marchman_2018", "fernald_marchman_2012", "fmw_2013", "fernald_totlot")
admins <- readRDS(here("data/raw/peekbank/_pb2026_admins.rds")) %>%
  distinct(dataset_name, subject_id, administration_id, lab_subject_id)
pk <- readRDS(here("data/raw/peekbank/1_d_sub.Rds")) %>% filter(dataset_name %in% OURS) %>%
  left_join(admins, by = c("dataset_name", "subject_id", "administration_id")) %>%
  transmute(dataset = dataset_name, child_id = as.character(lab_subject_id), age,
            log_rt, rt_ms = rt, accuracy = long_window_accuracy, n_trials_rt)
sd <- read_csv(here("data/raw/seedlings/seedlings_lwl_rt.csv"), show_col_types = FALSE) %>%
  transmute(dataset = "seedlings", child_id = as.character(lab_subject_id), age = lwl_age,
            log_rt = lwl_log_rt, rt_ms = exp(lwl_log_rt), accuracy = NA_real_, n_trials_rt = NA_integer_)
lwl <- bind_rows(pk, sd) %>%
  mutate(paper_code = recode(dataset, adams_marchman_2018 = "AM2018", fernald_marchman_2012 = "FM2012",
                             fmw_2013 = "FMW2013", fernald_totlot = "FPM2006", seedlings = "SEEDLingS")) %>%
  select(dataset, paper_code, child_id, age, log_rt, rt_ms, accuracy, n_trials_rt) %>%
  filter(!is.na(log_rt)) %>% arrange(paper_code, child_id, age)
write_csv(lwl, here("data/harmonized/lwl_session_level.csv"))

cat(sprintf("wrote input_level.csv (%d rows) + lwl_session_level.csv (%d rows)\n\n", nrow(input), nrow(lwl)))
cat("=== INPUT by dataset ===\n")
input %>% group_by(paper_code, instrument) %>%
  summarise(kids = n_distinct(child_id), recs = n(), age_lo = min(age, na.rm=TRUE), age_hi = max(age, na.rm=TRUE), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n=== LWL (session) by dataset ===\n")
lwl %>% group_by(paper_code) %>%
  summarise(kids = n_distinct(child_id), sessions = n(), has_acc = sum(!is.na(accuracy)) > 0,
            age_lo = min(age, na.rm=TRUE), age_hi = max(age, na.rm=TRUE), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
