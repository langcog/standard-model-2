## Cross-link Peekbank LWL data with the Stanford TotLot item-level CDIs.
##
## RUN LOCALLY ONLY. Reads the freshly-pulled LWL summary that already
## carries lab_subject_id (peekbank_2022_lwl_summary.csv); see
## model/scripts/pull_peekbank_lwl.R for how that file is built.
##
## Inputs:
##   data/peekbank/peekbank_2022_lwl_summary.csv  (per LWL admin
##     with lab_subject_id, age, RT, accuracy)
##   data/peekbank/stanford_cdi_items_long.csv     (TL2/TL3 items)
##
## Outputs:
##   data/peekbank/peekbank_stanford_linked.csv     (per LWL
##     admin: processing measures + nearest-age CDI summary)
##   data/peekbank/peekbank_stanford_admin_match.csv (diagnostic)

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

OUT_DIR <- file.path(PROJECT_ROOT, "data/peekbank")

# ---- 1. LWL admins with lab_subject_id + processing measures ---- #
lwl <- read_csv(file.path(OUT_DIR, "peekbank_2022_lwl_summary.csv"),
                show_col_types = FALSE, progress = FALSE) %>%
  filter(!is.na(lab_subject_id), nzchar(lab_subject_id)) %>%
  mutate(lab_subject_id = as.character(lab_subject_id),
         lwl_age = age) %>%
  rename(lwl_short_acc = short_window_accuracy,
         lwl_long_acc  = long_window_accuracy,
         lwl_rt        = rt,
         lwl_log_rt    = log_rt) %>%
  select(dataset_name, peekbank_subject_id, peekbank_admin_id,
         lab_subject_id, sex, lwl_age,
         n_trials, n_trials_rt,
         lwl_short_acc, lwl_long_acc, lwl_rt, lwl_log_rt)

cat(sprintf("LWL admins: %d (kids: %d)\n",
            nrow(lwl), n_distinct(lwl$lab_subject_id)))

# ---- 2. CDI admins (TL3 maps to adams_marchman_2018 via lab IDs) ---- #
cdi_admins <- read_csv(file.path(OUT_DIR, "stanford_cdi_items_long.csv"),
                       show_col_types = FALSE, progress = FALSE) %>%
  mutate(lab_subject_id = as.character(lab_subject_id)) %>%
  group_by(lab_subject_id, age, form) %>%
  summarise(n_items = n(), n_produces = sum(produces), .groups = "drop") %>%
  rename(cdi_age = age)
cat(sprintf("CDI admins: %d (kids: %d)\n",
            nrow(cdi_admins), n_distinct(cdi_admins$lab_subject_id)))

# ---- 3. Inner join LWL <-> CDI by lab_subject_id, then nearest age ---- #
linked <- lwl %>%
  inner_join(cdi_admins, by = "lab_subject_id",
             relationship = "many-to-many") %>%
  mutate(cdi_age_diff = abs(lwl_age - cdi_age)) %>%
  group_by(dataset_name, peekbank_admin_id) %>%
  filter(cdi_age_diff == min(cdi_age_diff)) %>%
  slice_head(n = 1) %>%
  ungroup()

cat(sprintf("\nLinked admins: %d (kids: %d)\n",
            nrow(linked), n_distinct(linked$lab_subject_id)))
cat("By dataset:\n")
print(linked %>% group_by(dataset_name) %>%
      summarise(n_admins = n(), n_kids = n_distinct(lab_subject_id),
                lwl_age_range = sprintf("%d-%d",
                                         min(lwl_age, na.rm=TRUE),
                                         max(lwl_age, na.rm=TRUE))))
cat("\ncdi_age_diff distribution:\n"); print(summary(linked$cdi_age_diff))

# ---- 4. Output: one row per LWL admin ---- #
linked_out <- linked %>%
  rename(cdi_form = form,
         cdi_n_items = n_items, cdi_n_produces = n_produces) %>%
  select(dataset_name, lab_subject_id, sex, lwl_age, cdi_age, cdi_age_diff,
         n_trials, n_trials_rt, lwl_short_acc, lwl_long_acc,
         lwl_rt, lwl_log_rt,
         cdi_form, cdi_n_items, cdi_n_produces,
         peekbank_admin_id)

write_csv(linked_out, file.path(OUT_DIR, "peekbank_stanford_linked.csv"))
cat(sprintf("\nWrote peekbank_stanford_linked.csv: %d admins (%d kids)\n",
            nrow(linked_out), n_distinct(linked_out$lab_subject_id)))

# Diagnostic: every LWL admin's nearest CDI within the same kid
linked %>%
  select(dataset_name, lab_subject_id, lwl_age, cdi_age, cdi_age_diff,
         n_produces) %>%
  arrange(lab_subject_id, lwl_age) %>%
  write_csv(file.path(OUT_DIR, "peekbank_stanford_admin_match.csv"))
cat("Wrote peekbank_stanford_admin_match.csv (diagnostic)\n")
