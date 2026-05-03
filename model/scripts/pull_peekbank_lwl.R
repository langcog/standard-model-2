## Pull LWL processing measures from peekbank 2022.1 directly,
## keyed by the lab_subject_id that links cleanly to the Stanford
## TotLot CDI files.
##
## Why this exists: peekbank-development's cached d_sub.Rds is from
## peekbank 2026.1 (not currently hosted), and its subject_id values
## don't bridge to peekbankr 2022.1's subject_id. Without an external
## bridge, joining d_sub to 2022.1's lab_subject_id by anything other
## than subject_id silently mixes children of similar ages. This
## pulls fresh from 2022.1 to get a clean (lab_subject_id, age) ->
## RT/accuracy mapping.
##
## RUN LOCALLY ONLY. Connects to peekbank.stanford.edu.
##
## Output: data/peekbank/peekbank_2022_lwl_summary.csv
##   Columns: dataset_name, lab_subject_id, peekbank_subject_id,
##            peekbank_admin_id, age, sex,
##            n_trials, n_trials_rt,
##            short_window_accuracy, long_window_accuracy,
##            rt, log_rt

source("model/R/config.R")
suppressPackageStartupMessages({
  library(peekbankr); library(dplyr); library(tidyr); library(readr)
})

OUT_DIR <- file.path(PROJECT_ROOT, "data/peekbank")
TARGET_DATASETS <- c("adams_marchman_2018", "fmw_2013")

# RT calculator (ported from peekbank-development helper/rt_helper.R)
get_rt <- function(rle_data, SAMPLING_RATE = 40) {
  if (is.null(rle_data$values) || is.null(rle_data$lengths))
    return(tibble(rt = NA_real_, shift_type = NA_character_))
  onset_aoi <- rle_data$values[1]
  if (!(onset_aoi %in% c("target", "distractor")))
    return(tibble(rt = NA_real_, shift_type = "other"))
  first_landing <- rle_data$values[rle_data$values != onset_aoi &
                                     rle_data$values %in% c("target", "distractor")][1]
  if (is.na(first_landing))
    return(tibble(rt = NA_real_, shift_type = "no shift"))
  shift_type <- if (onset_aoi == "distractor" && first_landing == "target") "D-T"
                else if (onset_aoi == "target" && first_landing == "distractor") "T-D"
                else "other"
  first_landing_idx <- which(rle_data$values == first_landing)[1]
  values_before <- rle_data$lengths[1:(first_landing_idx - 1)]
  rt <- (sum(values_before) + 1) * (1000 / SAMPLING_RATE)
  tibble(rt = rt, shift_type = shift_type)
}

cat("Connecting to peekbank...\n")
con <- connect_to_peekbank(db_version = "current")
admins <- collect(get_administrations(connection = con))
subj   <- collect(get_subjects(connection = con))
trials <- collect(get_trials(connection = con))
trial_types <- collect(get_trial_types(connection = con))
aoi    <- collect(get_aoi_timepoints(connection = con, rle = FALSE))
DBI::dbDisconnect(con)
cat(sprintf("aoi rows: %d  trials: %d  admins: %d  subj: %d\n",
            nrow(aoi), nrow(trials), nrow(admins), nrow(subj)))

# Filter to target datasets + English vanilla trials
admins_t <- admins %>% filter(dataset_name %in% TARGET_DATASETS) %>%
  inner_join(subj %>% select(subject_id, lab_subject_id, sex),
             by = "subject_id")
## Note: peekbank 2022.1 trial_types does not have a `vanilla_trial`
## flag (added in a later version), and trials does not carry
## administration_id (admin attaches via aoi_timepoints). We join from
## the aoi side: aoi has both trial_id and administration_id.
en_trial_types <- trial_types %>% filter(full_phrase_language == "eng") %>%
  select(trial_type_id, target_id, distractor_id, target_side)
aoi_t <- aoi %>%
  inner_join(admins_t %>% select(administration_id, dataset_name),
             by = "administration_id") %>%
  inner_join(trials %>% select(trial_id, trial_type_id),
             by = "trial_id") %>%
  inner_join(en_trial_types, by = "trial_type_id")

cat(sprintf("After filtering: %d aoi rows for %d admins\n",
            nrow(aoi_t), n_distinct(aoi_t$administration_id)))

# Code correctness
aoi_t <- aoi_t %>%
  mutate(correct = case_when(aoi == "target" ~ 1L,
                              aoi == "distractor" ~ 0L,
                              TRUE ~ NA_integer_))

# Per-trial accuracy summaries
trial_summary <- aoi_t %>%
  mutate(short_w = t_norm > 200 & t_norm <= 2000,
         long_w  = t_norm > 200 & t_norm <= 4000) %>%
  group_by(administration_id, trial_id) %>%
  summarise(
    short_window_accuracy = mean(correct[short_w], na.rm = TRUE),
    long_window_accuracy  = mean(correct[long_w],  na.rm = TRUE),
    short_window_prop_data = sum(!is.na(correct[short_w])) /
                              max(1, length(correct[short_w])),
    long_window_prop_data  = sum(!is.na(correct[long_w])) /
                              max(1, length(correct[long_w])),
    .groups = "drop"
  ) %>%
  mutate(across(contains("accuracy"), ~ifelse(is.nan(.x), NA, .x)),
         short_window_accuracy = ifelse(short_window_prop_data >= 0.5,
                                         short_window_accuracy, NA),
         long_window_accuracy = ifelse(long_window_prop_data >= 0.5,
                                        long_window_accuracy, NA))

# Per-trial RT via RLE on aoi sequence (only post-onset, only D-T shifts)
cat("Computing per-trial RTs (this is the slow step)...\n")
rle_data <- aoi_t %>%
  filter(t_norm >= 0) %>%
  group_by(administration_id, trial_id) %>%
  arrange(t_norm, .by_group = TRUE) %>%
  summarise(rle_obj = list(rle(aoi)), .groups = "drop")

rt_summary <- rle_data %>%
  rowwise() %>%
  mutate(rt_info = list(get_rt(rle_obj))) %>%
  ungroup() %>%
  unnest(rt_info) %>%
  mutate(rt = ifelse(shift_type == "D-T" & !is.na(rt) & rt >= 367, rt, NA),
         log_rt = log(rt))

# Combine to per-trial summary
per_trial <- trial_summary %>%
  left_join(rt_summary %>% select(administration_id, trial_id, rt, log_rt),
            by = c("administration_id", "trial_id"))

# Aggregate to per-admin
per_admin <- per_trial %>%
  group_by(administration_id) %>%
  summarise(
    n_trials = sum(!is.na(long_window_accuracy)),
    n_trials_rt = sum(!is.na(rt)),
    short_window_accuracy = mean(short_window_accuracy, na.rm = TRUE),
    long_window_accuracy  = mean(long_window_accuracy,  na.rm = TRUE),
    rt = mean(rt, na.rm = TRUE),
    log_rt = mean(log_rt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(rt, log_rt), ~ifelse(is.nan(.x), NA, .x)),
         across(contains("accuracy"), ~ifelse(is.nan(.x), NA, .x)),
         # Apply min-trial filters as in peekbank-development
         rt = ifelse(n_trials_rt < 2, NA, rt),
         log_rt = ifelse(n_trials_rt < 2, NA, log_rt),
         across(contains("accuracy"), ~ifelse(n_trials < 4, NA, .x)))

# Attach admin metadata
out <- admins_t %>%
  select(dataset_name, peekbank_subject_id = subject_id,
         peekbank_admin_id = administration_id,
         age, lab_subject_id, sex) %>%
  inner_join(per_admin %>% rename(peekbank_admin_id = administration_id),
             by = "peekbank_admin_id") %>%
  mutate(lab_subject_id = as.character(lab_subject_id))

cat(sprintf("\nFinal per-admin rows: %d (kids: %d)\n",
            nrow(out), n_distinct(out$lab_subject_id)))
print(out %>% group_by(dataset_name) %>%
      summarise(n_admins = n(), n_kids = n_distinct(lab_subject_id),
                age_range = sprintf("%d-%d", min(age), max(age)),
                pct_with_rt = round(100 * mean(!is.na(rt)), 1)))

write_csv(out, file.path(OUT_DIR, "peekbank_2022_lwl_summary.csv"))
cat(sprintf("Wrote %s\n", file.path(OUT_DIR, "peekbank_2022_lwl_summary.csv")))
