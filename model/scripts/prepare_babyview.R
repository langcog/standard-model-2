## Prepare BabyView (English-dominant longitudinal subjects) for the
## input-observed model. Builds a Stan-ready bundle that includes the
## standard CDI side AND per-video log-token-rate observations.
##
## RUN LOCALLY ONLY. The 875 MB merged_transcripts_parsed.csv is too
## big to ship to Sherlock. The prepared bundle
## `model/fits/babyview_subset_data.rds` is small (~500 KB) and is
## committed to the repo, so Sherlock just reads it directly.
##
## Usage:   Rscript model/scripts/prepare_babyview.R [n_items]
## Default: n_items = 200 (stratified by class x difficulty quartile)
##
## Inputs:
##   data/raw_data/babyview/video_metadata_processed.csv
##   data/raw_data/babyview/merged_transcripts_parsed.csv
##   data/raw_data/babyview/cdi_data_oct_2025/babyview-english-{ws,wg}_items.csv
##   model/fits/long_items.rds  (for English CHILDES p_j; reused)
##
## Output:  model/fits/babyview_subset_data.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
n_items <- as.integer(if (length(args) >= 1) args[1] else 200)

# ---- Constants ----
ENG_PCT_THRESHOLD <- 80    # videos: keep if percent_english >= this
SCORE_THRESHOLD   <- 0.3   # transcripts: drop tokens below this confidence
MIN_ADMINS        <- 2     # subjects: keep only longitudinal
N_DIFF_BINS       <- 4
SEED              <- 20260428

BV_DIR <- file.path(PROJECT_ROOT, "data/raw_data/babyview")
CDI_DIR <- file.path(BV_DIR, "cdi_data_oct_2025")

message(sprintf("Preparing BabyView bundle (n_items=%d)", n_items))

# ---- 1. Video metadata + transcript: per-video log token rate ----
v_meta <- read_csv(file.path(BV_DIR, "video_metadata_processed.csv"),
                   show_col_types = FALSE, progress = FALSE)
v_meta <- v_meta %>%
  filter(percent_english >= ENG_PCT_THRESHOLD,
         !is.na(duration_sec), duration_sec > 60) %>%
  mutate(subject_num = as.integer(gsub("^S0*", "", subject_id)))

message(sprintf("  videos: %d English-dominant, from %d subjects",
                nrow(v_meta), n_distinct(v_meta$subject_num)))

# Per-video adult-token counts. Only FEM + MAL speakers, alpha tokens,
# transcription confidence >= SCORE_THRESHOLD.
message("  reading transcript (slow)...")
tr <- read_csv(file.path(BV_DIR, "merged_transcripts_parsed.csv"),
               show_col_types = FALSE, progress = FALSE,
               col_select = c(video_id, speaker, score, spacy_is_alpha))
tr <- tr %>%
  filter(speaker %in% c("FEM", "MAL"),
         spacy_is_alpha,
         !is.na(score), score >= SCORE_THRESHOLD)

video_tokens <- tr %>%
  count(video_id, name = "n_adult_tokens")

# Combine: video-level log-rate
videos <- v_meta %>%
  inner_join(video_tokens, by = "video_id") %>%
  mutate(log_r_obs = log(pmax(n_adult_tokens, 1) /
                         (duration_sec / 3600)))   # tokens/hour, then log

message(sprintf("  videos with adult tokens: %d (median %.0f tokens, log_r_obs median %.2f)",
                nrow(videos), median(videos$n_adult_tokens),
                median(videos$log_r_obs)))

# ---- 2. CDI: WS + WG completed admins for English-dominant subjects ----
ws <- read_csv(file.path(CDI_DIR, "babyview-english-ws_items.csv"),
               show_col_types = FALSE, progress = FALSE)
wg <- read_csv(file.path(CDI_DIR, "babyview-english-wg_items.csv"),
               show_col_types = FALSE, progress = FALSE)

# WG and WS have different schemas but share a lot of item columns.
# Identify item columns by name pattern: lowercase letter (or letter+suffix)
# and not a known meta column.
META_COLS <- c("opt_out","study_name","subject_id","local_lab_id","repeat_num",
               "administration_id","link","completed","completedBackgroundInfo",
               "due_date","last_modified","created_date","completed_date",
               "event_id","age","sex","country","zip_code","birth_order",
               "birth_weight_lb","birth_weight_kg","multi_birth_boolean",
               "multi_birth","sibling_boolean","sibling_count","sibling_data",
               "born_on_due_date","early_or_late","due_date_diff",
               "form_filler","form_filler_other","primary_caregiver",
               "primary_caregiver_other","mother_yob","mother_education",
               "secondary_caregiver","secondary_caregiver_other",
               "father_yob","father_education","annual_income",
               "child_hispanic_latino","child_ethnicity","caregiver_info",
               "caregiver_other","other_languages_boolean","other_languages",
               "language_from","language_days_per_week","language_hours_per_day",
               "ear_infections_boolean","ear_infections","hearing_loss_boolean",
               "hearing_loss","vision_problems_boolean","vision_problems",
               "illnesses_boolean","illnesses","services_boolean","services",
               "worried_boolean","worried","learning_disability_boolean",
               "learning_disability","children_comforted","show_respect",
               "close_bonds","parents_help_learn","play_learning",
               "explore_experiment","do_as_told","read_at_home","teach_alphbet",
               "rhyming_games","read_for_pleasure","child_asks_for_reading",
               "child_self_reads","child_asks_words_say","place_of_residence",
               "primary_caregiver_occupation","primary_caregiver_occupation_description",
               "secondary_caregiver_occupation","secondary_caregiver_occupation_description",
               "kindergarten_since_when","kindergarten_hpd","kindergarten_dpw")

# Pivot to long form: one row per (admin, item)
to_long <- function(df, form_label) {
  meta_present <- intersect(META_COLS, colnames(df))
  item_cols <- setdiff(colnames(df), meta_present)
  # Drop obviously-non-item trailing columns
  item_cols <- item_cols[!grepl("^Total|^Word|^Combining|^Complexity|benchmark|^How |^Combination |Percentile|Combining",
                                item_cols)]
  df %>%
    filter(completed, !is.na(age), age != 999, age >= 8, age <= 36) %>%
    select(subject_id, age, all_of(item_cols)) %>%
    # Coerce every item column to character so pivot_longer doesn't choke
    # on mixed-type columns ("produces" strings vs. numeric markers).
    mutate(across(-c(subject_id, age), as.character)) %>%
    pivot_longer(-c(subject_id, age),
                 names_to = "item", values_to = "raw") %>%
    mutate(produces = as.integer(!is.na(raw) & raw == "produces"),
           form = form_label) %>%
    select(subject_id, age, form, item, produces)
}

cdi_long <- bind_rows(to_long(ws, "WS"), to_long(wg, "WG")) %>%
  filter(!is.na(produces))

message(sprintf("  CDI long-format rows: %d (admins: %d, items: %d)",
                nrow(cdi_long),
                n_distinct(paste(cdi_long$subject_id, cdi_long$age, cdi_long$form)),
                n_distinct(cdi_long$item)))

# ---- 3. Subset to subjects with both video AND >=2 CDI admins ----
cdi_subjs <- cdi_long %>%
  distinct(subject_id, age) %>%
  count(subject_id, name = "n_admins") %>%
  filter(n_admins >= MIN_ADMINS) %>%
  pull(subject_id)

video_subjs <- videos %>% distinct(subject_num) %>% pull(subject_num)
keep_subjs <- intersect(cdi_subjs, video_subjs)
message(sprintf("  subjects with video AND >=%d admins: %d", MIN_ADMINS,
                length(keep_subjs)))

cdi_long <- cdi_long %>% filter(subject_id %in% keep_subjs)
videos   <- videos %>% filter(subject_num %in% keep_subjs)

# ---- 4. Attach CHILDES p_j to items ----
# Re-use the Wordbank long_items processed English prob table.
long_ws <- readRDS(file.path(PATHS$fits_dir, "long_items.rds")) %>%
  filter(language == "English (American)")
prob_tbl <- long_ws %>% distinct(item, prob) %>% filter(!is.na(prob), prob > 0)

# BabyView item names use underscore-separated lowercase ('peanut_butter');
# Wordbank items use mixed punctuation. Build a normalization fn.
normalize_item <- function(x) {
  x %>% tolower() %>%
    gsub("[\\.\\(\\)]", "_", ., perl = TRUE) %>%
    gsub("\\s+|/", "_", .) %>%
    gsub("__+", "_", .) %>%
    gsub("[^a-z0-9_]", "", .) %>%
    sub("_+$", "", .)
}
prob_tbl <- prob_tbl %>% mutate(item_norm = normalize_item(item))
cdi_long <- cdi_long %>% mutate(item_norm = normalize_item(item))

# Use the most-frequent prob per normalized name
prob_lookup <- prob_tbl %>% group_by(item_norm) %>%
  summarise(prob = max(prob), .groups = "drop")

cdi_long <- cdi_long %>% left_join(prob_lookup, by = "item_norm")
matched <- cdi_long %>% distinct(item, prob) %>% summarise(
  total = n(), with_prob = sum(!is.na(prob)))
message(sprintf("  CDI items with CHILDES match: %d / %d",
                matched$with_prob, matched$total))

# Drop items without prob; lexical_class isn't in the CDI item CSVs, so
# we'll need to fetch it from the Wordbank long_ws_items table.
class_lookup <- long_ws %>% distinct(item, lexical_category) %>%
  mutate(item_norm = normalize_item(item)) %>%
  group_by(item_norm) %>% slice(1) %>% ungroup() %>%
  select(item_norm, lexical_category)
cdi_long <- cdi_long %>% left_join(class_lookup, by = "item_norm") %>%
  filter(!is.na(prob), !is.na(lexical_category))

message(sprintf("  after filtering to items with prob + class: %d obs, %d items",
                nrow(cdi_long), n_distinct(cdi_long$item)))

# ---- 5. Stratified item subsample (class x difficulty quartile) ----
set.seed(SEED)
item_summary <- cdi_long %>%
  distinct(item, lexical_category, prob) %>%
  mutate(log_p = log(prob)) %>%
  group_by(lexical_category) %>%
  mutate(diff_bin = ntile(log_p, N_DIFF_BINS)) %>%
  ungroup()
n_classes <- n_distinct(item_summary$lexical_category)
per_cell <- max(1, floor(n_items / (n_classes * N_DIFF_BINS)))
chosen_items <- item_summary %>%
  group_by(lexical_category, diff_bin) %>%
  slice_sample(n = per_cell) %>%
  ungroup() %>%
  pull(item)
if (length(chosen_items) > n_items) chosen_items <- sample(chosen_items, n_items)
if (length(chosen_items) < n_items) {
  extras <- setdiff(item_summary$item, chosen_items)
  need <- min(n_items - length(chosen_items), length(extras))
  if (need > 0) chosen_items <- c(chosen_items, sample(extras, need))
}
cdi_long <- cdi_long %>% filter(item %in% chosen_items)

message(sprintf("  items kept: %d across %d classes x %d difficulty bins",
                n_distinct(cdi_long$item), n_classes, N_DIFF_BINS))

# ---- 6. Build Stan-ready bundle ----
cdi_long <- cdi_long %>%
  mutate(admin_key = paste(subject_id, age, form, sep = "_"),
         aa = as.integer(factor(admin_key)),
         ii = as.integer(factor(subject_id)),
         jj = as.integer(factor(item)),
         cc = as.integer(factor(lexical_category)))

admin_info <- cdi_long %>% distinct(aa, ii, age, admin_key) %>% arrange(aa)
word_info  <- cdi_long %>% group_by(jj) %>%
  summarise(item = first(item), prob = first(prob), cc = first(cc),
            .groups = "drop") %>% arrange(jj)
class_levels <- levels(factor(cdi_long$lexical_category))

I <- max(cdi_long$ii); A <- max(cdi_long$aa)
J <- max(cdi_long$jj); C <- max(cdi_long$cc)

# Map videos: subject_num -> ii
subj_to_ii <- cdi_long %>% distinct(subject_id, ii) %>%
  rename(subject_num = subject_id)
videos <- videos %>% inner_join(subj_to_ii, by = "subject_num") %>%
  rename(child_ii = ii)
V <- nrow(videos)

message(sprintf("  bundle: I=%d, A=%d, J=%d, C=%d, N=%d, V=%d",
                I, A, J, C, nrow(cdi_long), V))

# Sperry-anchored prior on mu_r
prior_r <- load_input_rate_prior()

stan_data <- c(
  list(
    N = nrow(cdi_long),
    A = A, I = I, J = J, C = C,
    aa = cdi_long$aa, jj = cdi_long$jj,
    admin_to_child = admin_info$ii,
    cc = word_info$cc,
    y  = cdi_long$produces,
    admin_age = admin_info$age,
    log_p = log(word_info$prob),
    log_H = MODEL_CONSTANTS$log_H,
    a0    = MODEL_CONSTANTS$a0,

    # Input-observed side
    V = V,
    video_to_child = videos$child_ii,
    log_r_obs = videos$log_r_obs,
    log_r_obs_weight = rep(0, V),    # unused for now

    # Population-input priors: anchored to Sperry, with weak SD
    mu_r_prior_mean = prior_r$mu_r,
    mu_r_prior_sd   = 1,
    sigma_r_prior_sd = 1,

    # Reactivity prior: videos ~1.4x real life, 95% CI ~[0, +0.8]
    beta_react_prior_mean = 0.4,
    beta_react_prior_sd   = 0.4,

    # Within-child noise (per-video deviation around child's true rate)
    sigma_within_prior_sd = 1
  ),
  DEFAULT_PRIORS
)

bundle <- list(
  stan_data = stan_data,
  admin_info = admin_info,
  word_info = word_info,
  class_levels = class_levels,
  videos = videos,
  df = cdi_long,
  language = "English (BabyView)"
)

saveRDS(bundle, file.path(PATHS$fits_dir, "babyview_subset_data.rds"))
cat("\nSaved model/fits/babyview_subset_data.rds\n")
