## Compute per-recording and per-child input-rate estimates from the
## local datasets used in this project (Sperry/Hart-Risley/Weisleder-Fernald
## pooled CSV, BabyView, SEEDLingS).
##
## Outputs
##   input_estimation/local_per_recording.csv  — one row per recording/dyad
##   input_estimation/local_per_child.csv      — one row per child (geometric mean across recordings)
##   input_estimation/local_summary.csv        — one row per (dataset, measure_type) group
##
## Units. Everything is reported in two parallel columns:
##   tokens_per_hour  (raw tokens of adult speech per hour of waking observation)
##   log_r            (natural log of tokens_per_hour)
##
## tokens_per_month is derived in build_validation_set.R via the model's
## H = 365 hr/month convention (12 waking hr/day * 30.44 days/mo).
##
## Usage:  Rscript input_estimation/compute_local_rates.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

ROOT <- Sys.getenv("STANDARD_MODEL_ROOT", unset = "/Users/mcfrank/Projects/standard_model_2")
OUT  <- file.path(ROOT, "input_estimation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Sperry / Hart & Risley / Weisleder & Fernald pooled CSV
# ------------------------------------------------------------
# This is the existing pooled CSV that already drives the model's
# external prior on log r_i. We re-process it here to extract per-sample
# group means and SDs, broken out by:
#   - CDS from mother only      (mother_child_tokens_hr)
#   - CDS from any adult        (adult_child_tokens_hr)
#   - all adult speech in room  (all_tokens_hr; CDS + ODS)

sperry_raw <- read_csv(file.path(ROOT,
  "data/sperry/hourly_tokens_Sperry_HartRisley.csv"),
  show_col_types = FALSE, progress = FALSE)

# Long format so we can summarize by measure_type
sperry_long <- sperry_raw %>%
  pivot_longer(c(mother_child_tokens_hr, adult_child_tokens_hr, all_tokens_hr),
               names_to = "raw_measure", values_to = "tokens_per_hour") %>%
  filter(!is.na(tokens_per_hour), tokens_per_hour > 0) %>%
  mutate(
    measure_type = recode(raw_measure,
      "mother_child_tokens_hr" = "CDS-mother",
      "adult_child_tokens_hr"  = "CDS-any-adult",
      "all_tokens_hr"          = "all-adult"),
    log_r = log(tokens_per_hour),
    source = case_when(
      dataset == "Sperry"                  ~ "Sperry et al. 2019",
      dataset == "Hart & Risley"           ~ "Hart & Risley 1995",
      dataset == "Weisleder & Fernald"     ~ "Weisleder & Fernald 2013",
      dataset == "Soderstrom & Wittebolle" ~ "Soderstrom & Wittebolle 2013",
      TRUE ~ dataset),
    citation_path = case_when(
      dataset == "Sperry"                  ~ "papers/input_estimation/sperry_etal_2019.pdf",
      dataset == "Hart & Risley"           ~ NA_character_,
      dataset == "Weisleder & Fernald"     ~ "papers/input_estimation/weisleder_fernald_2013.pdf",
      dataset == "Soderstrom & Wittebolle" ~ NA_character_,
      TRUE ~ NA_character_),
    language = case_when(
      dataset == "Sperry"                  ~ "English (US)",
      dataset == "Hart & Risley"           ~ "English (US, KS)",
      dataset == "Weisleder & Fernald"     ~ "Spanish (US, low-SES)",
      dataset == "Soderstrom & Wittebolle" ~ "English (CA, Manitoba)",
      TRUE ~ "English (US)"),
    method = "manual transcript (visit recording)") %>%
  rename(sample_label = sample) %>%
  select(source, sample_label, language, measure_type, method,
         tokens_per_hour, log_r, citation_path)

# ------------------------------------------------------------
# 2. SEEDLingS — per-recording LENA AWC
# ------------------------------------------------------------
# 44 children x ~13 monthly LENA recordings, ages 6-17 months. AWC =
# Adult Word Count, an automated count of all adult words near the
# child's microphone (CDS + ODS lumped). We follow the same filter as
# prepare_seedlings.R: month >= 6 & <= 17, !awc_outlier.
lena <- read_csv(file.path(ROOT, "data/seedlings/lena_data.csv"),
                 show_col_types = FALSE, progress = FALSE)

seedlings_rec <- lena %>%
  filter(!awc_outlier, awc_perhr > 0, month >= 6, month <= 17) %>%
  transmute(
    source = "this paper (SEEDLingS; Bergelson 2018; Egan-Dailey & Bergelson 2025)",
    sample_label = "44 longitudinal Rochester subjects",
    language = "English (US, Rochester NY)",
    measure_type = "all-adult",
    method = "LENA AWC (automated)",
    subject_id = subj,
    age_mo = month,
    tokens_per_hour = awc_perhr,
    log_r = log(awc_perhr),
    citation_path = "papers/input_estimation/bergelson_etal_2018.pdf")

# ------------------------------------------------------------
# 3. BabyView — per-video FEM+MAL transcript adult-token rate
# ------------------------------------------------------------
# 20 longitudinal English-dominant subjects x 5,688 head-cam videos,
# ages 8-30 months. Adult tokens = Whisper-transcribed FEM + MAL alpha
# tokens with confidence >= 0.3, normalized by video duration. Counts
# all adult speech in the recording (CDS + ODS).
bv <- readRDS(file.path(ROOT, "model/fits/babyview_subset_data.rds"))

babyview_rec <- bv$videos %>%
  transmute(
    source = "this paper (BabyView; Long et al. 2025)",
    sample_label = "20 longitudinal English-dominant Stanford subjects",
    language = "English (US, CA)",
    measure_type = "all-adult",
    method = "FEM+MAL transcript (Whisper, score>=0.3)",
    subject_id = as.character(subject_id),
    age_mo = age_mo,
    tokens_per_hour = exp(log_r_obs),
    log_r = log_r_obs,
    citation_path = "papers/input_estimation/bergelson_etal_2018.pdf")  # placeholder -- see notes

# Long et al. 2025 isn't in the input_estimation/ folder (it's the
# BabyView dataset paper, not really an input-estimation paper). Cite
# it via the data README. Set a clearer pointer:
babyview_rec$citation_path <- "data/babyview/README.md"

# ------------------------------------------------------------
# 4. Concatenate per-recording table
# ------------------------------------------------------------
# Sperry CSV rows are already per-dyad (one row per family/recording).
sperry_per_rec <- sperry_long %>%
  mutate(subject_id = NA_character_,
         age_mo = NA_real_) %>%
  select(source, sample_label, language, measure_type, method,
         subject_id, age_mo, tokens_per_hour, log_r, citation_path)

per_recording <- bind_rows(sperry_per_rec, seedlings_rec, babyview_rec) %>%
  mutate(across(c(tokens_per_hour, log_r), ~ round(.x, 4)))

write_csv(per_recording, file.path(OUT, "local_per_recording.csv"))
message(sprintf("[per_recording] %d rows -> %s",
                nrow(per_recording), "input_estimation/local_per_recording.csv"))

# ------------------------------------------------------------
# 5. Per-child geometric mean (where multiple recordings per child)
# ------------------------------------------------------------
# For SEEDLingS and BabyView, average log_r across recordings within a
# child to get one per-child estimate. Sperry/HR/WF rows are already
# per-dyad with no within-dyad replication, so they pass through.
per_child_local <- per_recording %>%
  filter(!is.na(subject_id)) %>%
  group_by(source, sample_label, language, measure_type, method, subject_id, citation_path) %>%
  summarise(n_recordings = n(),
            log_r = mean(log_r),
            tokens_per_hour = exp(log_r),
            .groups = "drop")

per_child_sperry <- sperry_per_rec %>%
  mutate(subject_id = paste(source, sample_label, row_number(), sep = "_"),
         n_recordings = 1) %>%
  select(source, sample_label, language, measure_type, method,
         subject_id, n_recordings, log_r, tokens_per_hour, citation_path)

per_child <- bind_rows(per_child_sperry, per_child_local) %>%
  mutate(across(c(tokens_per_hour, log_r), ~ round(.x, 4)))

write_csv(per_child, file.path(OUT, "local_per_child.csv"))
message(sprintf("[per_child]     %d rows -> %s",
                nrow(per_child), "input_estimation/local_per_child.csv"))

# ------------------------------------------------------------
# 6. Per-(source, sample, measure) summary
# ------------------------------------------------------------
# For each subgroup, report:
#   - n_kids, age_min/max
#   - geometric-mean tokens/hr  ( = exp(mean(log_r)) )
#   - SD of log_r across kids (this is the σ_r quantity)
local_summary <- per_child %>%
  group_by(source, sample_label, language, measure_type, method, citation_path) %>%
  summarise(
    n_kids = n_distinct(subject_id),
    log_r_mean = mean(log_r),
    log_r_sd   = if (n() >= 2) sd(log_r) else NA_real_,
    tokens_per_hour_mean = exp(log_r_mean),
    tokens_per_hour_median = median(exp(log_r)),
    age_range_mo = if (any(!is.na(per_child$age_mo))) {
      "see local_per_recording.csv"   # placeholder; resolved below
    } else NA_character_,
    .groups = "drop") %>%
  mutate(across(c(log_r_mean, log_r_sd,
                  tokens_per_hour_mean, tokens_per_hour_median),
                ~ round(.x, 4)))

# Resolve age ranges from the per-recording table
age_ranges <- per_recording %>%
  group_by(source, sample_label, measure_type) %>%
  summarise(age_min = suppressWarnings(min(age_mo, na.rm = TRUE)),
            age_max = suppressWarnings(max(age_mo, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(age_range_mo = ifelse(is.finite(age_min),
                               sprintf("%g-%g", round(age_min,1), round(age_max,1)),
                               NA_character_)) %>%
  select(source, sample_label, measure_type, age_range_mo)

local_summary <- local_summary %>%
  select(-age_range_mo) %>%
  left_join(age_ranges, by = c("source", "sample_label", "measure_type"))

# Pooled aggregate that mirrors what `load_input_rate_prior()` uses
# as the model's external prior on log r_i: all rows where
# adult_child_tokens_hr is non-NA. (HR and WF have NA in that column,
# so this is effectively Sperry only -- documenting it explicitly.)
pooled <- per_child %>%
  filter(measure_type == "CDS-any-adult") %>%
  summarise(
    source = "POOLED (model's external prior on log r_i)",
    sample_label = sprintf("all %d adult_child_tokens_hr non-NA rows", n()),
    language = "English (US, mixed sites)",
    measure_type = "CDS-any-adult",
    method = "manual transcript (visit recording)",
    n_kids = n(),
    log_r_mean = round(mean(log_r), 4),
    log_r_sd   = round(sd(log_r), 4),
    tokens_per_hour_mean = round(exp(mean(log_r)), 4),
    tokens_per_hour_median = round(median(exp(log_r)), 4),
    age_range_mo = NA_character_,
    citation_path = "data/sperry/hourly_tokens_Sperry_HartRisley.csv")

local_summary <- bind_rows(local_summary, pooled) %>%
  arrange(source, measure_type, sample_label)

write_csv(local_summary, file.path(OUT, "local_summary.csv"))
message(sprintf("[summary]       %d rows -> %s",
                nrow(local_summary), "input_estimation/local_summary.csv"))

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------
cat("\n=== Local-data input-rate summary (tokens/hour, log scale) ===\n")
local_summary %>%
  mutate(label = sprintf("%s | %s | %s", source, sample_label, measure_type)) %>%
  select(label, n_kids, log_r_mean, log_r_sd,
         tokens_per_hour_mean, tokens_per_hour_median,
         age_range_mo) %>%
  print(n = Inf, width = Inf)

cat("\nNotes:\n",
"  - log_r_sd is the SD of log(per-child mean tokens/hr) ACROSS children.\n",
"    This is the σ_r quantity that the model treats as the population\n",
"    SD of log r_i (see model_explainer.pdf §4.3).\n",
"  - For the local datasets (BabyView, SEEDLingS), per-child means\n",
"    are geometric across the recordings/videos belonging to that child.\n",
"  - SEEDLingS uses LENA AWC (auto, all adult words).\n",
"  - BabyView uses FEM+MAL transcript tokens (auto, all adult words).\n",
"  - Sperry/HR/WF measure_type='CDS-any-adult' = adult_child_tokens_hr.\n",
sep = "")
