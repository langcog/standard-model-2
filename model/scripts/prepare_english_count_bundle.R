## English-only bundle WITH the sumscore count branch (no Spanish): the full English io-proc MM
## bundle + the ELENA-WS English sumscore admins (the "no-item-level" English kids). Isolates
## whether the English count cohort alone gives stable estimates, separate from Spanish/HABLA.
## Same count machinery as the bilingual bundle, one form (English-WS = all English items).
## RUN LOCALLY. -> fits/joint_io_proc_english_count_subset_data.rds
suppressPackageStartupMessages({library(here); library(dplyr); library(readxl)})

b  <- readRDS(here("fits/joint_io_proc_mm_subset_data.rds"))
sd <- b$stan_data; I0 <- sd$I; A0 <- sd$A; J0 <- sd$J; S0 <- sd$S; a0 <- sd$a0
ci <- b$child_info %>% mutate(subject_id = as.character(subject_id))
cat(sprintf("English base: I=%d J=%d S=%d N=%d N_lwl=%d V_obs=%d\n", I0, J0, S0, sd$N, sd$N_lwl, sd$V_obs))

## ELENA-WS sumscores -> attach to existing fmw kids where ids match, else new English (study 3)
el <- read_excel(here("data/raw/FMW2013/elena/ELENA_WS_SumScores.xlsx")) %>%
  transmute(subject_id = as.character(ParticipantId), age = as.numeric(CDIAge), k = as.numeric(VOCAB)) %>%
  filter(!is.na(k), !is.na(age))
eng_fmw <- ci %>% filter(dataset_name == "fmw_2013") %>% transmute(subject_id = as.character(lab_subject_id), ii)
el <- el %>% left_join(eng_fmw, by = "subject_id")
elena_new <- sort(unique(el$subject_id[is.na(el$ii)]))
new_kids <- tibble(subject_id = elena_new, dataset_name = "fmw_2013", study = 3L, ii = I0 + seq_along(elena_new))
I <- I0 + nrow(new_kids)
study_of_child <- c(sd$study_of_child, rep(3L, nrow(new_kids)))
el <- el %>% left_join(new_kids %>% select(subject_id, ii_new = ii), by = "subject_id") %>%
  mutate(ii = coalesce(ii, ii_new)) %>% select(-ii_new)
cat(sprintf("ELENA-WS count: %d admins, %d kids (%d matched existing + %d new English)\n",
            nrow(el), n_distinct(el$subject_id), n_distinct(el$subject_id) - length(elena_new), length(elena_new)))

## count branch: one form = the English item bank
sum_adm <- el %>% transmute(ii, age, k = pmin(k, J0), form = 1L) %>% filter(!is.na(ii), k >= 0)
stan_data <- modifyList(sd, list(
  I = I, S = S0, study_of_child = study_of_child,
  n_sum = nrow(sum_adm), sum_child = sum_adm$ii, sum_log_age = log(sum_adm$age / a0),
  sum_k = as.integer(sum_adm$k), sum_form = as.integer(sum_adm$form),
  n_forms = 1L, n_form_items = J0, form_item = as.integer(seq_len(J0)),
  form_start = 1L, form_len = as.integer(J0)))

child_info <- bind_rows(ci, new_kids %>% transmute(ii, lab_subject_id = subject_id, dataset_name, study, subject_id)) %>% arrange(ii)
bundle <- list(stan_data = stan_data, word_info = b$word_info, child_info = child_info,
               lwl = b$lwl, datasets = b$datasets, input_recs = b$input_recs,
               language = "English + sumscore count branch (ELENA-WS); no Spanish")
out <- here("fits/joint_io_proc_english_count_subset_data.rds")
saveRDS(bundle, out)
cat(sprintf("\n== English+count bundle ==\n  I=%d A=%d J=%d C=%d S=%d N=%d N_lwl=%d V_obs=%d n_sum=%d\n  Saved %s\n",
            I, sd$A, J0, sd$C, S0, sd$N, sd$N_lwl, sd$V_obs, nrow(sum_adm), out))
