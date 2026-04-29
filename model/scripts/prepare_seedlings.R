## Prepare SEEDLingS (Bergelson lab) for the input-observed model.
##
## RUN LOCALLY ONLY. The prepared bundle
## `model/fits/seedlings_subset_data.rds` is small (~50 KB) and gets
## committed to the repo so Sherlock can read it without re-running
## any of the prep upstream.
##
## Usage:   Rscript model/scripts/prepare_seedlings.R [n_items]
## Default: n_items = 200 (stratified by class x difficulty quartile)
##
## Public inputs (already in repo):
##   data/raw_data/seedlings/lena_data.csv      Per-recording LENA AWC
##                                              (44 kids x 13 monthly recordings)
##   data/raw_data/seedlings/seedlings_data.csv Per-child summary (used for QA)
##
## NOT YET PUBLIC, must be obtained from Bergelson lab:
##   data/raw_data/seedlings/cdi_items_long.csv Item-level CDI WG responses at
##                                              1;0 and 1;6 for the 44 SEEDLingS
##                                              subjects. See README in the same
##                                              folder for the expected schema.
##
## Reused inputs:
##   model/fits/long_items.rds  (English CHILDES p_j and lexical_category)
##
## Output:  model/fits/seedlings_subset_data.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
n_items <- as.integer(if (length(args) >= 1) args[1] else 200)

# ---- Constants ----
N_DIFF_BINS <- 4
SEED        <- 20260429
SD_DIR      <- file.path(PROJECT_ROOT, "data/raw_data/seedlings")
CDI_LONG    <- file.path(SD_DIR, "cdi_items_long.csv")

message(sprintf("Preparing Seedlings bundle (n_items=%d)", n_items))

# ---- 1. Per-recording LENA AWC -> log_r_obs ----
lena <- read_csv(file.path(SD_DIR, "lena_data.csv"),
                 show_col_types = FALSE, progress = FALSE) %>%
  filter(!awc_outlier, awc_perhr > 0,
         # Infancy window only; drop the 4;6 followup recording.
         month >= 6, month <= 17)
# subj is "01".."44", left-padded
recordings <- lena %>%
  transmute(subject_id = subj,
            month = month,
            log_r_obs = log(awc_perhr))
message(sprintf("  recordings: %d (children: %d, months %d-%d)",
                nrow(recordings),
                n_distinct(recordings$subject_id),
                min(recordings$month), max(recordings$month)))

# ---- 2. CDI item-level data (gated) ----
if (!file.exists(CDI_LONG)) {
  stop(sprintf(paste0(
    "Item-level CDI data not yet available at %s. ",
    "We have only CDI totals from the public Egan-Dailey github. ",
    "Item-level CDI for SEEDLingS lives in the private cdi_spreadsheet ",
    "repo at BergelsonLab. Once you have a wordlevel-format export, ",
    "drop it at the path above (see data/raw_data/seedlings/README.md ",
    "for expected schema) and rerun."), CDI_LONG))
}

cdi_long <- read_csv(CDI_LONG, show_col_types = FALSE, progress = FALSE)
required_cols <- c("subject_id", "age", "form", "item", "produces")
missing <- setdiff(required_cols, colnames(cdi_long))
if (length(missing)) {
  stop(sprintf("cdi_items_long.csv missing columns: %s",
               paste(missing, collapse = ", ")))
}
cdi_long <- cdi_long %>%
  filter(form == "WG", !is.na(produces)) %>%
  mutate(produces = as.integer(produces))

message(sprintf("  CDI: %d rows (admins: %d, items: %d)",
                nrow(cdi_long),
                n_distinct(paste(cdi_long$subject_id, cdi_long$age, cdi_long$form)),
                n_distinct(cdi_long$item)))

# ---- 3. Subset to subjects with both LENA and >=1 CDI admin ----
keep <- intersect(unique(recordings$subject_id),
                  unique(cdi_long$subject_id))
message(sprintf("  subjects with LENA AND CDI: %d", length(keep)))
cdi_long   <- cdi_long  %>% filter(subject_id %in% keep)
recordings <- recordings %>% filter(subject_id %in% keep)

# ---- 4. Attach English CHILDES p_j + lexical class ----
long_ws <- readRDS(file.path(PATHS$fits_dir, "long_items.rds")) %>%
  filter(language == "English (American)")

normalize_item <- function(x) {
  x %>% tolower() %>%
    gsub("[\\.\\(\\)]", "_", ., perl = TRUE) %>%
    gsub("\\s+|/", "_", .) %>%
    gsub("__+", "_", .) %>%
    gsub("[^a-z0-9_]", "", .) %>%
    sub("_+$", "", .)
}

prob_lookup <- long_ws %>%
  distinct(item, prob) %>% filter(!is.na(prob), prob > 0) %>%
  mutate(item_norm = normalize_item(item)) %>%
  group_by(item_norm) %>% summarise(prob = max(prob), .groups = "drop")

class_lookup <- long_ws %>%
  distinct(item, lexical_category) %>%
  mutate(item_norm = normalize_item(item)) %>%
  group_by(item_norm) %>% slice(1) %>% ungroup() %>%
  select(item_norm, lexical_category)

cdi_long <- cdi_long %>%
  mutate(item_norm = normalize_item(item)) %>%
  left_join(prob_lookup, by = "item_norm") %>%
  left_join(class_lookup, by = "item_norm") %>%
  filter(!is.na(prob), !is.na(lexical_category))

message(sprintf("  after CHILDES match: %d obs, %d items",
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
per_cell  <- max(1, floor(n_items / (n_classes * N_DIFF_BINS)))
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

# Map recordings -> ii via subject_id
subj_to_ii <- cdi_long %>% distinct(subject_id, ii)
recordings <- recordings %>% inner_join(subj_to_ii, by = "subject_id") %>%
  rename(child_ii = ii)
V <- nrow(recordings)

message(sprintf("  bundle: I=%d, A=%d, J=%d, C=%d, N=%d, V=%d",
                I, A, J, C, nrow(cdi_long), V))

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

    # Input-observed side: LENA AWC per recording
    V = V,
    video_to_child = recordings$child_ii,
    log_r_obs = recordings$log_r_obs,
    log_r_obs_weight = rep(0, V),

    # Population-input priors: anchored to Sperry, with weak SD
    mu_r_prior_mean = prior_r$mu_r,
    mu_r_prior_sd   = 1,
    sigma_r_prior_sd = 1,

    # Reactivity prior weaker than BabyView: at-home audio recording
    # has less observer effect than a body-cam video.
    beta_react_prior_mean = 0.1,
    beta_react_prior_sd   = 0.4,

    sigma_within_prior_sd = 1
  ),
  DEFAULT_PRIORS
)

bundle <- list(
  stan_data = stan_data,
  admin_info = admin_info,
  word_info = word_info,
  class_levels = class_levels,
  recordings = recordings,
  df = cdi_long,
  language = "English (Seedlings)"
)

saveRDS(bundle, file.path(PATHS$fits_dir, "seedlings_subset_data.rds"))
cat("\nSaved model/fits/seedlings_subset_data.rds\n")
