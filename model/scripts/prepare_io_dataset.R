## Prepare an IO (input-observed) bundle for the Marchman/Fernald-lab
## cohorts that share the TL3TLO_LENA.csv input file + the parsed
## stanford_cdi_items_long.csv CDI file.
##
## RUN LOCALLY ONLY (needs no wordbankr; reads only repo files).
##
## Usage:  Rscript model/scripts/prepare_io_dataset.R <dataset>
##   <dataset> in {am2018, fmw2013}
##
## Design choices for this re-prep (input-uptake revisited):
##   * WG + WS combined (production) — full age range, no WS-only window.
##   * ALL items (no stratified subsample) — these samples are small, so
##     keep every word the kids were scored on.
##   * beta_react pinned at 0 (mean=0, sd=0.001) — removed for simplicity.
##   * delta_j anchored to the big EN longitudinal fit via an informative
##     per-item prior (delta_j_prior_mean / _sd), so the small IO samples
##     spend their data on the input/efficiency/slope params. Items with
##     no EN anchor (~21 WG-specific words) get a weak prior at the EN
##     population mean.
##
## Output: fits/io_<dataset>_subset_data.rds
##         (kept name "subset" for loader compatibility, though it's now
##          the full item set)

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) tolower(args[1]) else "am2018"
spec <- switch(DATASET,
  am2018  = list(study = "totlot3", lena_study = "TL3", lena_ages = c(16, 18),
                  language = "English (AM2018)"),
  fmw2013 = list(study = "tlo",     lena_study = "TLO", lena_ages = c(18),
                  language = "English (FMW2013)"),
  stop("dataset must be am2018 or fmw2013")
)
cat(sprintf("=== Preparing IO bundle: %s (study=%s, LENA=%s) ===\n",
            DATASET, spec$study, spec$lena_study))

PK_DIR  <- file.path(PROJECT_ROOT, "data/peekbank")
OUT_RDS <- file.path(PATHS$fits_dir, sprintf("io_%s_subset_data.rds", DATASET))

## ---- 1. CDI (WG+WS) for this study --------------------------------
cdi <- read_csv(file.path(PK_DIR, "stanford_cdi_items_long.csv"),
                show_col_types = FALSE) |>
  filter(study == spec$study) |>
  mutate(subject_id = as.character(lab_subject_id))
cat(sprintf("  CDI: %d rows, %d kids, forms %s, ages %d-%d\n",
            nrow(cdi), n_distinct(cdi$subject_id),
            paste(sort(unique(cdi$form)), collapse="+"),
            min(cdi$age), max(cdi$age)))

## ---- 2. LENA input -> per-recording log_r_obs ---------------------
## TL3TLO_LENA.csv is wide: SubjectID1, Study, AWCHr16M, AWCHr18M, ...
lena_wide <- read_csv(file.path(PK_DIR, "TL3TLO_LENA.csv"),
                      show_col_types = FALSE) |>
  filter(Study == spec$lena_study) |>
  mutate(subject_id = as.character(SubjectID1))
# pivot the AWC-per-hour columns we want (16M and/or 18M) to long
awc_cols <- paste0("AWCHr", spec$lena_ages, "M")
recordings <- lena_wide |>
  select(subject_id, all_of(awc_cols)) |>
  pivot_longer(all_of(awc_cols), names_to = "which", values_to = "awc_perhr") |>
  filter(!is.na(awc_perhr), awc_perhr > 0) |>
  mutate(log_r_obs = log(awc_perhr))
cat(sprintf("  LENA: %d recordings, %d kids (AWC/hr median=%.0f)\n",
            nrow(recordings), n_distinct(recordings$subject_id),
            exp(median(recordings$log_r_obs))))

## ---- 3. Restrict to kids with BOTH CDI and LENA -------------------
keep <- intersect(unique(cdi$subject_id), unique(recordings$subject_id))
cdi        <- cdi        |> filter(subject_id %in% keep)
recordings <- recordings |> filter(subject_id %in% keep)
cat(sprintf("  kids with CDI AND LENA: %d\n", length(keep)))

## ---- 4. Attach CHILDES p_j + lexical class ------------------------
long_ws <- readRDS(file.path(PATHS$fits_dir, "long_items.rds")) |>
  filter(language == "English (American)")
normalize_item <- function(x) {
  x |> tolower() |>
    gsub("[\\.\\(\\)]", "_", x = _, perl = TRUE) |>
    gsub("\\s+|/", "_", x = _) |>
    gsub("__+", "_", x = _) |>
    gsub("[^a-z0-9_]", "", x = _) |>
    sub("_+$", "", x = _)
}
freq <- readRDS(file.path(PATHS$fits_dir, "english_word_freq.rds"))
prob_lookup <- freq |>
  mutate(item_norm = normalize_item(w)) |>
  group_by(item_norm) |>
  summarise(count = sum(count), .groups = "drop") |>
  mutate(prob = count / sum(count)) |>
  select(item_norm, prob)
class_lookup <- long_ws |>
  distinct(item, lexical_category) |>
  mutate(item_norm = normalize_item(item)) |>
  group_by(item_norm) |> slice(1) |> ungroup() |>
  select(item_norm, lexical_category)

## Keep ALL items. Frequency (log_p) and class (cc) are vestigial in the
## anchored no_freq model — log_p enters only as beta_c*log_p with
## beta_c pinned at 0, and delta_j is anchored (not class-derived). So
## items that fail the CHILDES freq/class normalization (multi-word like
## "teddy bear", sense-split like "chicken (animal)") get placeholder
## covariates rather than being dropped. They still receive their EN
## delta_j anchor, which joins on the exact item string.
cdi <- cdi |>
  mutate(item_norm = normalize_item(item)) |>
  left_join(prob_lookup, by = "item_norm") |>
  left_join(class_lookup, by = "item_norm")
n_no_prob  <- sum(is.na(cdi$prob[!duplicated(cdi$item)]))
n_no_class <- sum(is.na(cdi$lexical_category[!duplicated(cdi$item)]))
med_prob <- median(cdi$prob, na.rm = TRUE)
cdi <- cdi |>
  mutate(prob = ifelse(is.na(prob), med_prob, prob),
         lexical_category = ifelse(is.na(lexical_category), "other",
                                    lexical_category))
cat(sprintf("  CHILDES match: %d obs, %d items kept (%d placeholder freq, %d placeholder class)\n",
            nrow(cdi), n_distinct(cdi$item), n_no_prob, n_no_class))

## ---- 5. delta_j anchor from the big EN longitudinal fit -----------
## EN long fit delta_j (posterior median) keyed by item_definition.
en_bundle <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))
en_psi    <- read_csv(file.path(PATHS$fits_dir, "summaries",
                                 "long_no_freq_slopes_psi.csv"),
                      show_col_types = FALSE)
en_delta  <- en_bundle$word_info |>
  select(jj, item) |>
  inner_join(en_psi |> select(jj, delta_j_median), by = "jj") |>
  select(item, delta_anchor = delta_j_median)
en_mu_delta <- median(en_delta$delta_anchor)   # population fallback
cat(sprintf("  EN anchor: %d items, mu_delta=%.2f\n",
            nrow(en_delta), en_mu_delta))

## ---- 6. Build Stan-ready bundle (all items) -----------------------
cdi <- cdi |>
  mutate(admin_key = paste(subject_id, age, form, sep = "_"),
         aa = as.integer(factor(admin_key)),
         ii = as.integer(factor(subject_id)),
         jj = as.integer(factor(item)),
         cc = as.integer(factor(lexical_category)))

admin_info <- cdi |> distinct(aa, ii, age, admin_key) |> arrange(aa)
word_info  <- cdi |> group_by(jj) |>
  summarise(item = first(item), prob = first(prob), cc = first(cc),
            .groups = "drop") |> arrange(jj)
class_levels <- levels(factor(cdi$lexical_category))
I <- max(cdi$ii); A <- max(cdi$aa); J <- max(cdi$jj); C <- max(cdi$cc)

# per-item delta anchor aligned to jj order; weak prior for unanchored
word_info <- word_info |>
  left_join(en_delta, by = "item") |>
  mutate(anchored = !is.na(delta_anchor),
         delta_j_prior_mean = ifelse(anchored, delta_anchor, en_mu_delta),
         delta_j_prior_sd   = ifelse(anchored, 0.10, 5.0))
cat(sprintf("  delta_j anchored: %d / %d items (%d weak-prior fallback)\n",
            sum(word_info$anchored), J, sum(!word_info$anchored)))

subj_to_ii <- cdi |> distinct(subject_id, ii)
recordings <- recordings |> inner_join(subj_to_ii, by = "subject_id") |>
  rename(child_ii = ii)
V <- nrow(recordings)

cat(sprintf("  bundle: I=%d, A=%d, J=%d, C=%d, N=%d, V=%d\n",
            I, A, J, C, nrow(cdi), V))

prior_r <- load_input_rate_prior()

stan_data <- c(
  list(
    N = nrow(cdi), A = A, I = I, J = J, C = C,
    aa = cdi$aa, jj = cdi$jj,
    admin_to_child = admin_info$ii,
    cc = word_info$cc,
    y  = cdi$produces,
    admin_age = admin_info$age,
    log_p = log(word_info$prob),
    log_H = MODEL_CONSTANTS$log_H,
    a0    = round(median(admin_info$age)),

    # comprehension / std channels OFF for these datasets
    N_comp = 0L, aa_comp = integer(0), jj_comp = integer(0), y_comp = integer(0),
    N_std = 0L, std_to_child = integer(0), std_score = numeric(0),

    # input channel: LENA AWC/hr per recording
    V = V,
    video_to_child = recordings$child_ii,
    log_r_obs = recordings$log_r_obs,
    log_r_obs_weight = rep(0, V),

    mu_r_prior_mean = prior_r$mu_r,
    mu_r_prior_sd   = 1,
    sigma_r_prior_sd = 1,

    # beta_react REMOVED (pinned at 0): passive LENA, no reactivity term.
    beta_react_prior_mean = 0,
    beta_react_prior_sd   = 0.001,
    sigma_within_prior_sd = 1,

    # delta_j anchor (informative per-item prior from EN longitudinal fit)
    delta_j_prior_mean = word_info$delta_j_prior_mean,
    delta_j_prior_sd   = word_info$delta_j_prior_sd
  ),
  DEFAULT_PRIORS
)

bundle <- list(
  stan_data = stan_data,
  admin_info = admin_info,
  word_info = word_info,
  class_levels = class_levels,
  recordings = recordings,
  df = cdi,
  language = spec$language
)
saveRDS(bundle, OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
