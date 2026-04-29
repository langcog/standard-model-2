## Prepare Peekbank longitudinal subjects for the processing-differences
## analysis (per-child RT and accuracy trajectories alongside CDI).
##
## RUN LOCALLY ONLY. Reads pre-built per-admin summaries lifted from
## the peekbank-development repo (see data/raw_data/peekbank/README.md);
## no peekbankr / wordbankr calls happen here.
##
## Usage:   Rscript model/scripts/prepare_peekbank.R [min_admins]
## Default: min_admins = 2 (kept if a subject has >= this many admins
##          on top of the CDI-matching requirement)
##
## Inputs:
##   data/raw_data/peekbank/1_d_sub.Rds       per-(subj, admin) summary
##   data/raw_data/peekbank/0_cdi_subjects.Rds CDI-to-admin fuzzy join
##                                             (already left-joined into
##                                             d_sub; kept around for
##                                             provenance)
##
## Output:  model/fits/peekbank_subset_data.rds
##
## Bundle schema differs from the longitudinal IRT bundles: there's no
## item-level CDI here, only CDI totals (peekbankr returns rawscores).
## Stored as a frame for downstream correlational / SEM analyses.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
MIN_ADMINS <- as.integer(if (length(args) >= 1) args[1] else 2)

PB_DIR <- file.path(PROJECT_ROOT, "data/raw_data/peekbank")

message(sprintf("Preparing Peekbank bundle (min_admins=%d)", MIN_ADMINS))

d_sub <- readRDS(file.path(PB_DIR, "1_d_sub.Rds")) %>% ungroup()
message(sprintf("  d_sub: %d admins across %d subjects, %d datasets",
                nrow(d_sub), n_distinct(d_sub$subject_id),
                n_distinct(d_sub$dataset_name)))

# ---- Longitudinal + CDI filter ----
# Subjects must (a) have >= MIN_ADMINS admins AND (b) have at least one
# admin with a non-NA CDI score (prod or comp).
admin_counts <- d_sub %>%
  count(dataset_name, subject_id, name = "n_admins")

cdi_subjects <- d_sub %>%
  filter(!is.na(prod) | !is.na(comp)) %>%
  distinct(dataset_name, subject_id)

keep <- admin_counts %>%
  filter(n_admins >= MIN_ADMINS) %>%
  inner_join(cdi_subjects, by = c("dataset_name", "subject_id"))

message(sprintf("  longitudinal-with-CDI: %d subjects across %d datasets",
                nrow(keep), n_distinct(keep$dataset_name)))

d_long <- d_sub %>%
  inner_join(keep %>% select(dataset_name, subject_id),
             by = c("dataset_name", "subject_id"))

# ---- Per-subject summary ----
# A subject with stable processing has: low within-subject variance in
# log_rt across admins after detrending age, and stable accuracy. The
# downstream analysis cares about per-child means and per-child slopes.
subject_summary <- d_long %>%
  group_by(dataset_name, subject_id) %>%
  arrange(age, .by_group = TRUE) %>%
  summarise(
    n_admins      = n(),
    age_first     = first(age),
    age_last      = last(age),
    age_span      = last(age) - first(age),
    mean_log_rt   = mean(log_rt, na.rm = TRUE),
    mean_short_acc = mean(short_window_accuracy, na.rm = TRUE),
    mean_long_acc  = mean(long_window_accuracy,  na.rm = TRUE),
    n_admins_with_cdi = sum(!is.na(prod) | !is.na(comp)),
    mean_prod_pct = mean(prod, na.rm = TRUE),
    mean_comp_pct = mean(comp, na.rm = TRUE),
    .groups = "drop"
  )

message(sprintf("  median age span: %.1f mo (range %d-%d)",
                median(subject_summary$age_span),
                min(subject_summary$age_first),
                max(subject_summary$age_last)))

# ---- Bundle ----
bundle <- list(
  d_admin       = d_long,           # per-admin frame (RT, acc, CDI)
  subject_summary = subject_summary,
  description   = paste("Peekbank: longitudinal subjects (>=", MIN_ADMINS,
                        " admins) with at least one CDI measurement"),
  source        = "peekbank-development cached_intermediates 1_d_sub.Rds (2026.1)"
)

out_path <- file.path(PATHS$fits_dir, "peekbank_subset_data.rds")
saveRDS(bundle, out_path)
cat(sprintf("\nSaved %s (%d subjects, %d admins)\n",
            out_path, nrow(subject_summary), nrow(d_long)))
