## Build ONE pooled bundle over all four input-observed datasets for the
## hierarchical IO model (input-uptake revisited).
##
## Pools the four cleaned per-dataset bundles (BabyView, SEEDLingS,
## AM2018, FMW2013) into a single Stan data object with:
##   * a study index (1..4) per child + admin + recording
##   * a measurement-type index per study: 1 = head-cam (BabyView),
##     2 = LENA (the other three). sigma_within is grouped by this so
##     SEEDLingS's dense replication pins LENA noise and rescues sigma_r
##     for the sparse AM2018/FMW2013 (the whole point of pooling).
##   * a GLOBAL item index keyed by item string, with delta_j anchored
##     to the EN longitudinal fit (shared across studies).
##   * delta_j_prior_mean/_sd anchor vectors.
##
## Prereqs (rebuild these first with all items / >=1 admin):
##   prepare_babyview.R 9999 1
##   prepare_seedlings.R 9999        (>=1 admin: see its MIN_ADMINS)
##   prepare_io_dataset.R am2018 / fmw2013
##
## Output: fits/io_pooled_subset_data.rds

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

OUT_RDS <- file.path(PATHS$fits_dir, "io_pooled_subset_data.rds")

# Normalize item strings so the SAME word maps to one global item index
# across studies. BabyView writes "baa_baa"/"uh_oh"/"ice_cream"; the LENA
# studies write "baa baa"/"uh oh". Keying jj on the raw string splits
# these into separate items (752), splitting their delta_j anchor and
# failing to pool a kid's production across the variants.
norm_it <- function(x) {
  x <- tolower(x); x <- gsub("[\\.\\(\\)]", "_", x, perl = TRUE)
  x <- gsub("\\s+|/", "_", x); x <- gsub("__+", "_", x)
  x <- gsub("[^a-z0-9_]", "", x); sub("_+$", "", x)
}

# study -> (bundle file, measurement type). meas: 1=head-cam, 2=LENA.
STUDIES <- tibble::tribble(
  ~study,        ~s, ~meas, ~bundle,
  "BabyView",     1L, 1L,   "babyview_subset_data.rds",
  "SEEDLingS",    2L, 2L,   "seedlings_subset_data.rds",
  "AM2018",       3L, 2L,   "io_am2018_subset_data.rds",
  "FMW2013",      4L, 2L,   "io_fmw2013_subset_data.rds"
)

## ---- 1. Pull CDI long + per-recording input from each bundle -------
cdi_all <- list(); rec_all <- list(); subj_all <- list()
for (i in seq_len(nrow(STUDIES))) {
  st <- STUDIES[i, ]
  b  <- readRDS(file.path(PATHS$fits_dir, st$bundle))
  df <- b$df
  # Key children on the bundle's INTERNAL integer index (df$ii ==
  # recordings$child_ii by construction), not subject_id — BabyView's
  # df uses numeric ids while its videos use "S00..." strings, so a
  # subject_id join silently drops all of BabyView.
  df$ckey <- paste(st$study, df$ii, sep = "::")
  # Retain the original subject_id per ckey where the source bundle carries
  # one (SEEDLingS does -> "01".."46" with gaps). Needed for robust joins to
  # external per-child data (e.g. the SEEDLings LWL-RT channel), which must NOT
  # assume ckey "::N" == subject_id "0N" (io ii is a dense factor rank).
  if ("subject_id" %in% names(df))
    subj_all[[st$study]] <- df |> distinct(ckey, subject_id) |>
      mutate(subject_id = as.character(subject_id))
  cdi_all[[st$study]] <- df |>
    transmute(study = st$study, s = st$s, ckey,
              age, form, item, produces,
              prob = prob, lexical_category = lexical_category)

  rec <- if (!is.null(b$recordings)) b$recordings else b$videos
  rec_all[[st$study]] <- rec |>
    transmute(study = st$study, s = st$s, meas = st$meas,
              ckey = paste(st$study, child_ii, sep = "::"),
              log_r_obs)
}
cdi <- bind_rows(cdi_all)
rec <- bind_rows(rec_all)

## keep only children that have BOTH CDI and >=1 input recording
ckeys <- intersect(unique(cdi$ckey), unique(rec$ckey))
cdi <- cdi |> filter(ckey %in% ckeys)
rec <- rec |> filter(ckey %in% ckeys)

cat(sprintf("Pooled: %d children across %d studies\n",
            length(ckeys), nrow(STUDIES)))
print(cdi |> distinct(study, ckey) |> count(study, name = "n_kids"))

## ---- 2. Global indices --------------------------------------------
cdi <- cdi |>
  mutate(admin_key = paste(ckey, age, form, sep = "_"),
         item_norm = norm_it(item),
         ii = as.integer(factor(ckey)),
         aa = as.integer(factor(admin_key)),
         jj = as.integer(factor(item_norm)),
         cc = as.integer(factor(lexical_category)))

child_info <- cdi |> distinct(ii, ckey, s) |> arrange(ii) |>
  left_join(bind_rows(subj_all), by = "ckey")     # subject_id; NA where source lacks it
admin_info <- cdi |> distinct(aa, ii, age, admin_key) |> arrange(aa)
word_info  <- cdi |> group_by(jj) |>
  summarise(item = first(item), item_norm = first(item_norm),
            prob = first(prob), cc = first(cc),
            .groups = "drop") |> arrange(jj)

I <- max(cdi$ii); A <- max(cdi$aa); J <- max(cdi$jj); C <- max(cdi$cc)
S <- nrow(STUDIES)

# recordings -> global child index
rec <- rec |> inner_join(child_info |> select(ii, ckey), by = "ckey")
V <- nrow(rec)

# per-study measurement type, and per-child study
meas_of_study  <- STUDIES$meas
study_of_child <- child_info$s

cat(sprintf("  I=%d, A=%d, J=%d, C=%d, V=%d, S=%d, meas types=%d\n",
            I, A, J, C, V, S, n_distinct(meas_of_study)))

## ---- 3. delta_j anchor from the EN longitudinal fit ---------------
en_b   <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))
en_psi <- read_csv(file.path(PATHS$fits_dir, "summaries",
                              "long_no_freq_slopes_psi.csv"),
                   show_col_types = FALSE)
en_delta <- en_b$word_info |> select(jj, item) |>
  inner_join(en_psi |> select(jj, delta_j_median), by = "jj") |>
  mutate(item_norm = norm_it(item)) |>
  group_by(item_norm) |> slice(1) |> ungroup() |>
  select(item_norm, delta_anchor = delta_j_median)
en_mu <- median(en_delta$delta_anchor)

word_info <- word_info |>
  left_join(en_delta, by = "item_norm") |>
  mutate(anchored = !is.na(delta_anchor),
         delta_j_prior_mean = ifelse(anchored, delta_anchor, en_mu),
         delta_j_prior_sd   = ifelse(anchored, 0.10, 5.0))
cat(sprintf("  delta_j anchored: %d / %d (%d weak-prior)\n",
            sum(word_info$anchored), J, sum(!word_info$anchored)))

## ---- 4. Assemble Stan data ----------------------------------------
prior_r <- tryCatch(
  { source("model/R/helpers.R"); load_input_rate_prior() },
  error = function(e) list(mu_r = 7.34))

a0 <- round(median(admin_info$age))

stan_data <- c(
  list(
    N = nrow(cdi), A = A, I = I, J = J, C = C, S = S,
    n_meas = n_distinct(meas_of_study),
    aa = cdi$aa, jj = cdi$jj,
    admin_to_child = admin_info$ii,
    study_of_child = study_of_child,
    meas_of_study  = meas_of_study,
    cc = word_info$cc,
    y  = cdi$produces,
    admin_age = admin_info$age,
    log_p = log(word_info$prob),
    log_H = MODEL_CONSTANTS$log_H,
    a0    = a0,

    # input channel
    V = V,
    recording_to_child = rec$ii,
    study_of_recording = rec$s,
    log_r_obs = rec$log_r_obs,

    # priors
    mu_r_prior_mean = prior_r$mu_r, mu_r_prior_sd = 3,
    sigma_r_prior_sd = 1, sigma_within_prior_sd = 1,
    sigma_study_delta_prior_sd = 2,

    # delta_j anchor
    delta_j_prior_mean = word_info$delta_j_prior_mean,
    delta_j_prior_sd   = word_info$delta_j_prior_sd
  ),
  DEFAULT_PRIORS
)

bundle <- list(
  stan_data = stan_data,
  child_info = child_info,
  admin_info = admin_info,
  word_info = word_info,
  studies = STUDIES,
  df = cdi,
  recordings = rec,
  language = "pooled-IO (BabyView+SEEDLingS+AM2018+FMW2013)"
)
saveRDS(bundle, OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
