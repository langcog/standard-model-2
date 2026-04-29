## Build the Stanford-linked bundle for the joint vocab + LWL processing
## fit (`log_irt_long_proc.stan`).
##
## RUN LOCALLY ONLY. Reuses Wordbank's English CHILDES p_j table from
## `model/fits/long_items.rds`; no wordbankr/peekbankr calls at runtime.
##
## Inputs:
##   data/raw_data/peekbank/peekbank_stanford_linked.csv  (LWL admins +
##     subject linkage; built by link_peekbank_stanford.R)
##   data/raw_data/peekbank/stanford_cdi_items_long.csv   (item-level CDI
##     in canonical Wordbank item names; built by parse_stanford_cdi.R)
##   model/fits/long_items.rds                             (CHILDES p_j)
##
## Output: model/fits/stanford_linked_subset_data.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
n_items <- as.integer(if (length(args) >= 1) args[1] else 200)

PB_DIR <- file.path(PROJECT_ROOT, "data/raw_data/peekbank")
SEED   <- 20260429
N_DIFF_BINS <- 4

# ---- 1. Linked LWL admins ---- #
lwl <- read_csv(file.path(PB_DIR, "peekbank_stanford_linked.csv"),
                show_col_types = FALSE, progress = FALSE) %>%
  mutate(lab_subject_id = as.character(lab_subject_id)) %>%
  filter(!is.na(lwl_log_rt))    # need RT for the LWL channel

cat(sprintf("Linked LWL admins (with RT): %d (kids: %d)\n",
            nrow(lwl), n_distinct(lwl$lab_subject_id)))

# ---- 2. Item-level CDI for the linked kids ---- #
cdi <- read_csv(file.path(PB_DIR, "stanford_cdi_items_long.csv"),
                show_col_types = FALSE, progress = FALSE) %>%
  mutate(lab_subject_id = as.character(lab_subject_id)) %>%
  filter(lab_subject_id %in% unique(lwl$lab_subject_id),
         form == "WS")          # All TL3 admins are WS; restrict for purity

cat(sprintf("CDI rows for linked kids: %d  (admins: %d, items: %d)\n",
            nrow(cdi),
            n_distinct(paste(cdi$lab_subject_id, cdi$age)),
            n_distinct(cdi$item)))

# ---- 3. Attach CHILDES p_j ---- #
long_ws <- readRDS(file.path(PATHS$fits_dir, "long_items.rds")) %>%
  filter(language == "English (American)")
prob_lookup <- long_ws %>%
  distinct(item, prob) %>% filter(!is.na(prob), prob > 0)
class_lookup <- long_ws %>%
  distinct(item, lexical_category) %>% group_by(item) %>%
  slice(1) %>% ungroup()

cdi <- cdi %>%
  left_join(prob_lookup, by = "item") %>%
  left_join(class_lookup, by = "item") %>%
  filter(!is.na(prob), !is.na(lexical_category))
cat(sprintf("After CHILDES match: %d rows, %d items\n",
            nrow(cdi), n_distinct(cdi$item)))

# ---- 4. Stratified item subsample ---- #
set.seed(SEED)
item_summary <- cdi %>%
  distinct(item, lexical_category, prob) %>%
  mutate(log_p = log(prob)) %>%
  group_by(lexical_category) %>%
  mutate(diff_bin = ntile(log_p, N_DIFF_BINS)) %>%
  ungroup()
n_classes <- n_distinct(item_summary$lexical_category)
per_cell  <- max(1, floor(n_items / (n_classes * N_DIFF_BINS)))
chosen <- item_summary %>%
  group_by(lexical_category, diff_bin) %>%
  slice_sample(n = per_cell) %>%
  ungroup() %>% pull(item)
if (length(chosen) > n_items) chosen <- sample(chosen, n_items)
if (length(chosen) < n_items) {
  extras <- setdiff(item_summary$item, chosen)
  need <- min(n_items - length(chosen), length(extras))
  if (need > 0) chosen <- c(chosen, sample(extras, need))
}
cdi <- cdi %>% filter(item %in% chosen)
cat(sprintf("Items kept: %d across %d classes x %d difficulty bins\n",
            n_distinct(cdi$item), n_classes, N_DIFF_BINS))

# ---- 5. Build CDI-side stan_data ---- #
cdi <- cdi %>%
  mutate(admin_key = paste(lab_subject_id, age, sep = "_"),
         aa = as.integer(factor(admin_key)),
         ii = as.integer(factor(lab_subject_id)),
         jj = as.integer(factor(item)),
         cc = as.integer(factor(lexical_category)))

admin_info <- cdi %>% distinct(aa, ii, lab_subject_id, age, admin_key) %>% arrange(aa)
word_info  <- cdi %>% group_by(jj) %>%
  summarise(item = first(item), prob = first(prob), cc = first(cc),
            .groups = "drop") %>% arrange(jj)
class_levels <- levels(factor(cdi$lexical_category))

I <- max(cdi$ii); A <- max(cdi$aa)
J <- max(cdi$jj); C <- max(cdi$cc)

# Map LWL admins -> child index ii using lab_subject_id
subj_to_ii <- cdi %>% distinct(lab_subject_id, ii)
lwl <- lwl %>% inner_join(subj_to_ii, by = "lab_subject_id") %>%
  rename(child_ii = ii)
N_lwl <- nrow(lwl)

cat(sprintf("Bundle: I=%d, A=%d, J=%d, C=%d, N=%d, N_lwl=%d\n",
            I, A, J, C, nrow(cdi), N_lwl))

prior_r <- load_input_rate_prior()

stan_data <- c(
  list(
    # CDI side
    N = nrow(cdi),
    A = A, I = I, J = J, C = C,
    aa = cdi$aa, jj = cdi$jj,
    admin_to_child = admin_info$ii,
    cc = word_info$cc,
    y  = cdi$produces,
    admin_age = admin_info$age,
    log_p = log(word_info$prob),
    log_H = MODEL_CONSTANTS$log_H,
    a0    = MODEL_CONSTANTS$a0,
    mu_r = prior_r$mu_r,
    sigma_r = prior_r$sigma_r,

    # LWL side
    N_lwl = N_lwl,
    lwl_to_child = lwl$child_ii,
    lwl_log_age  = log(lwl$lwl_age / MODEL_CONSTANTS$a0),
    lwl_log_rt   = lwl$lwl_log_rt,

    # LWL priors. RT ~700 ms is typical at 22 mo, log(700) ~ 6.55.
    # mu_rtslope: log(rt) declines ~0.3 per log(age/a0); allow wide.
    # gamma_rt: log_alpha SD ~ 1; allow gamma in [-2, 2] effectively.
    mu_rt_prior_mean        = 6.5,
    mu_rt_prior_sd          = 1,
    mu_rtslope_prior_sd     = 1,
    gamma_rt_prior_sd       = 1,
    sigma_rtslope_prior_sd  = 1,
    sigma_lwl_prior_sd      = 1
  ),
  DEFAULT_PRIORS
)

bundle <- list(
  stan_data    = stan_data,
  admin_info   = admin_info,
  word_info    = word_info,
  class_levels = class_levels,
  lwl          = lwl,
  df           = cdi,
  language     = "English (Stanford-linked: Adams 2018)"
)

saveRDS(bundle, file.path(PATHS$fits_dir, "stanford_linked_subset_data.rds"))
cat("\nSaved model/fits/stanford_linked_subset_data.rds\n")
