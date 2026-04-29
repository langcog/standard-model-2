## Cross-link Peekbank LWL data with the Stanford TotLot item-level CDIs.
##
## RUN LOCALLY ONLY. Hits the live Peekbank DB to get lab_subject_id
## per Peekbank subject, then joins to the TL2/TL3 CDI files via that
## key.
##
## Inputs:
##   peekbankr (current DB version) get_subjects + get_administrations
##   data/raw_data/peekbank/1_d_sub.Rds                (LWL processing per admin)
##   data/raw_data/peekbank/stanford_cdi_items_long.csv (TL2/TL3 item-level CDI)
##
## Outputs:
##   data/raw_data/peekbank/peekbank_stanford_linked.csv
##       per (lab_subject_id, lwl_age) admin: RT/accuracy from Peekbank
##       + nearest CDI admin age and item-level production summary.
##   data/raw_data/peekbank/peekbank_stanford_admin_match.csv
##       diagnostic: which CDI admins matched which LWL admins.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(peekbankr)
})

OUT_DIR <- file.path(PROJECT_ROOT, "data/raw_data/peekbank")
TARGET_DATASETS <- c("adams_marchman_2018", "fmw_2013")
# fernald_marchman_2012 and fernald_totlot are only in db 2026.1 (not
# yet hosted by peekbank.stanford.edu), so we can't query them via the
# current peekbankr. Add them when 2026.1 comes online.

cat("Connecting to Peekbank DB...\n")
con   <- connect_to_peekbank(db_version = "current")
subj  <- collect(get_subjects(connection = con))
admin <- collect(get_administrations(connection = con))
DBI::dbDisconnect(con)

# ---- 1.  Peekbank lab_subject_id table for target datasets ---- #
peek_meta <- admin %>%
  filter(dataset_name %in% TARGET_DATASETS) %>%
  inner_join(subj %>% select(subject_id, lab_subject_id, sex),
             by = "subject_id") %>%
  mutate(lab_subject_id = as.character(lab_subject_id)) %>%
  select(dataset_name, peekbank_subject_id = subject_id,
         peekbank_admin_id = administration_id,
         lab_subject_id, lwl_age = age, sex)
cat(sprintf("Peekbank target admins: %d (across %d kids in %d datasets)\n",
            nrow(peek_meta),
            n_distinct(peek_meta$lab_subject_id),
            n_distinct(peek_meta$dataset_name)))

# ---- 2.  LWL summary measures from cached d_sub ---- #
d_sub <- readRDS(file.path(OUT_DIR, "1_d_sub.Rds")) %>% ungroup() %>%
  filter(dataset_name %in% TARGET_DATASETS) %>%
  select(dataset_name, subject_id, administration_id, age,
         n_trials, n_trials_rt,
         short_window_accuracy, long_window_accuracy,
         rt, log_rt, rt_var, log_rt_var,
         prod_pct = prod, comp_pct = comp)
cat(sprintf("d_sub admins: %d\n", nrow(d_sub)))

# Note: d_sub uses 2026.1 subject_ids which won't match peek_meta's
# 2022.1 subject_ids. Join by (dataset_name, age) and lab_subject_id
# is impossible without a bridge. Fortunately the lab_subject_id
# in peekbankr 2022.1 is stable across DB versions (it's the lab's
# original ID), so we use:
#
#   d_sub age  +  Peekbank 2022.1 admin age  +  lab_subject_id
#   (matched via subject_id only at 2022.1; lwl summary measures
#    keyed by 2026.1 subject_id are joined back via dataset+age
#    fuzzy match -- mostly reliable since LWL admin ages are unique
#    per kid within +/- 0.5 mo).

# Build a 2022.1-side admin lookup that also includes 2026.1-style
# d_sub admin age + dataset_name to enable a bridging fuzzy join.
peek_admins_2022 <- peek_meta %>%
  select(dataset_name, peekbank_admin_id, lab_subject_id, lwl_age, sex)

# Fuzzy join: for each peek_admin (2022.1) find d_sub row in same
# dataset within +/- 0.5 mo of lwl_age. Accept exact int match where
# possible; fall back to nearest within 1 mo otherwise.
match_admin <- peek_admins_2022 %>%
  inner_join(d_sub, by = c("dataset_name"),
             relationship = "many-to-many") %>%
  mutate(age_diff = abs(lwl_age - age)) %>%
  group_by(dataset_name, peekbank_admin_id) %>%
  filter(age_diff == min(age_diff)) %>%
  filter(age_diff <= 1) %>%
  slice_head(n = 1) %>%
  ungroup()
cat(sprintf("LWL admins matched to d_sub by age: %d / %d\n",
            nrow(match_admin), nrow(peek_admins_2022)))

# ---- 3.  Stanford TotLot CDI: nearest-age admin per kid ---- #
cdi_long <- read_csv(file.path(OUT_DIR, "stanford_cdi_items_long.csv"),
                     show_col_types = FALSE, progress = FALSE) %>%
  mutate(lab_subject_id = as.character(lab_subject_id))

cdi_admins <- cdi_long %>%
  group_by(lab_subject_id, age, form) %>%
  summarise(n_items = n(), n_produces = sum(produces), .groups = "drop")
cat(sprintf("CDI admins available for matching: %d (across %d kids)\n",
            nrow(cdi_admins), n_distinct(cdi_admins$lab_subject_id)))

# For each LWL admin, find the nearest-age CDI admin for the same
# lab_subject_id (within +/- 3 mo; CDIs are administered every few
# months, so up to 3 mo gap is normal). Rename CDI age before join
# to avoid name collisions with d_sub's age column.
cdi_admins_named <- cdi_admins %>% rename(cdi_age = age)
linked <- match_admin %>%
  inner_join(cdi_admins_named, by = "lab_subject_id",
             relationship = "many-to-many") %>%
  mutate(cdi_age_diff = abs(lwl_age - cdi_age)) %>%
  group_by(dataset_name, peekbank_admin_id) %>%
  filter(cdi_age_diff == min(cdi_age_diff)) %>%
  slice_head(n = 1) %>%
  ungroup()

cat(sprintf("LWL admins linked to a CDI admin: %d (kids: %d)\n",
            nrow(linked), n_distinct(linked$lab_subject_id)))
cat("Distribution of cdi_age_diff:\n"); print(summary(linked$cdi_age_diff))

# Final output: one row per LWL admin with linked CDI summary
linked_out <- linked %>%
  rename(lwl_short_acc = short_window_accuracy,
         lwl_long_acc  = long_window_accuracy,
         lwl_rt = rt, lwl_log_rt = log_rt,
         lwl_rt_var = rt_var, lwl_log_rt_var = log_rt_var,
         cdi_form = form,
         cdi_n_items = n_items, cdi_n_produces = n_produces) %>%
  select(dataset_name, lab_subject_id, sex, lwl_age, cdi_age, cdi_age_diff,
         n_trials, n_trials_rt, lwl_short_acc, lwl_long_acc,
         lwl_rt, lwl_log_rt, lwl_rt_var, lwl_log_rt_var,
         cdi_form, cdi_n_items, cdi_n_produces,
         peekbank_admin_id)

write_csv(linked_out, file.path(OUT_DIR, "peekbank_stanford_linked.csv"))
cat(sprintf("\nWrote peekbank_stanford_linked.csv: %d admins (%d kids)\n",
            nrow(linked_out), n_distinct(linked_out$lab_subject_id)))

# Also emit the LWL admin -> CDI admin mapping diagnostic for review
match_admin %>%
  left_join(cdi_admins_named, by = "lab_subject_id",
            relationship = "many-to-many") %>%
  mutate(cdi_age_diff = abs(lwl_age - cdi_age)) %>%
  group_by(dataset_name, peekbank_admin_id) %>%
  filter(cdi_age_diff == min(cdi_age_diff)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(dataset_name, lab_subject_id, lwl_age,
         cdi_age, cdi_age_diff, n_produces) %>%
  arrange(lab_subject_id, lwl_age) %>%
  write_csv(file.path(OUT_DIR, "peekbank_stanford_admin_match.csv"))
cat(sprintf("Wrote peekbank_stanford_admin_match.csv (diagnostic)\n"))
