## Prepare the longitudinal Wordbank data (WG + WS combined) for Stan
## fitting.
##
## Usage:
##   Rscript model/scripts/prepare_longitudinal_data.R [language] [n_children] [n_items]
## Defaults: English, 200 children, stratified subsample of 200 items.
##
## Reads:   fits/long_items.rds (from pull_longitudinal.R)
## Writes:  fits/long_subset_data.rds (Stan-ready bundle)
##
## Each (child, age, form) combination is its own admin so a child who
## took both WG and WS at the same age contributes two admin rows.

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
language   <- if (length(args) >= 1) args[1] else "English (American)"
n_children <- as.integer(if (length(args) >= 2) args[2] else 200)
n_items    <- as.integer(if (length(args) >= 3) args[3] else 200)
N_AGE_BINS  <- 4    # children stratified across age (median admin age)
N_DIFF_BINS <- 4    # items stratified across log_p quantiles within class
MIN_ADMINS  <- 2    # require >= this many admins per kept child
SEED <- 20260428

message(sprintf("Preparing longitudinal subset: language=%s, children=%d, items=%d",
                language, n_children, n_items))

long <- readRDS(file.path(PATHS$fits_dir, "long_items.rds"))
freq <- readRDS(file.path(PATHS$fits_dir, "english_word_freq.rds"))

# Filter to language; drop pre-existing prob column (legacy English
# preprocessed source or NA) so the fresh CHILDES freqs win the join.
d <- long %>%
  filter(language == !!language, !is.na(produces)) %>%
  select(-any_of("prob"))

# Wordbank item_definition often has parenthetical qualifiers (e.g.
# "go (verb)", "no (response)"); strip them for the lowercased join,
# parallel to prepare_longitudinal_norwegian.R.
normalize <- function(x) {
  x %>% tolower() %>%
    gsub("\\s*\\(.*\\)\\s*", "", .) %>%
    gsub("\\s+", " ", .) %>% trimws()
}

d <- d %>% mutate(item_norm = normalize(item))
freq_lookup <- freq %>%
  mutate(w_norm = normalize(w)) %>%
  group_by(w_norm) %>%
  summarise(count = sum(count), .groups = "drop")

total_freq <- sum(freq_lookup$count)
freq_lookup <- freq_lookup %>%
  mutate(prob = count / total_freq)

d <- d %>% left_join(freq_lookup, by = c("item_norm" = "w_norm"))

# Diagnostic: how many items matched?
matched <- d %>% distinct(item, item_norm, prob)
cat(sprintf("  items with CHILDES match: %d / %d\n",
            sum(!is.na(matched$prob)), nrow(matched)))

# For unmatched items, use the minimum freq (small but nonzero) so we
# keep the item rather than dropping it. Same convention as Norwegian.
if (sum(is.na(d$prob)) > 0) {
  min_prob <- min(d$prob, na.rm = TRUE)
  d <- d %>% mutate(prob = dplyr::coalesce(prob, min_prob))
}

d <- d %>% filter(prob > 0)

message(sprintf("  input: %d rows, %d children, %d items",
                nrow(d), length(unique(d$child_id)),
                length(unique(d$item))))

## Stratified sampling
##   Children: bin by median admin age into N_AGE_BINS, sample evenly.
##             Require >= MIN_ADMINS admins per child.
##   Items:    bin within (lexical_category x log_p quartile),
##             sample evenly so easy / hard words across all classes
##             are represented.
set.seed(SEED)

# --- Children ---
child_summary <- d %>%
  distinct(child_id, age) %>%
  group_by(child_id) %>%
  summarise(n_admin = n(),
            median_age = median(age),
            span = max(age) - min(age),
            .groups = "drop") %>%
  filter(n_admin >= MIN_ADMINS)

# Bin by median age (quantile bins, so each bin has roughly equal count)
child_summary <- child_summary %>%
  mutate(age_bin = ntile(median_age, N_AGE_BINS))

per_bin_kids <- max(1, floor(n_children / N_AGE_BINS))
chosen_children <- child_summary %>%
  group_by(age_bin) %>%
  # within an age bin, prefer kids with more admins (more longitudinal info)
  arrange(desc(n_admin), desc(span)) %>%
  slice_head(n = per_bin_kids) %>%
  ungroup() %>%
  pull(child_id)
# top-up if some bins were short (rare)
if (length(chosen_children) < n_children) {
  extras <- setdiff(child_summary$child_id, chosen_children)
  need <- min(n_children - length(chosen_children), length(extras))
  if (need > 0) chosen_children <- c(chosen_children, sample(extras, need))
}

message(sprintf("  children: %d kept (median admins=%.1f, median span=%.1f mo)",
                length(chosen_children),
                median(child_summary$n_admin[child_summary$child_id %in% chosen_children]),
                median(child_summary$span[child_summary$child_id %in% chosen_children])))

# --- Items ---
item_summary <- d %>%
  distinct(item, lexical_category, prob) %>%
  filter(!is.na(prob), prob > 0) %>%
  mutate(log_p = log(prob)) %>%
  group_by(lexical_category) %>%
  mutate(diff_bin = ntile(log_p, N_DIFF_BINS)) %>%
  ungroup()

n_classes <- length(unique(item_summary$lexical_category))
per_cell  <- max(1, floor(n_items / (n_classes * N_DIFF_BINS)))
chosen_items <- item_summary %>%
  group_by(lexical_category, diff_bin) %>%
  slice_sample(n = per_cell) %>%
  ungroup() %>%
  pull(item)
if (length(chosen_items) > n_items)
  chosen_items <- sample(chosen_items, n_items)
# top-up: if some cells had < per_cell items, sample to fill
if (length(chosen_items) < n_items) {
  extras <- setdiff(item_summary$item, chosen_items)
  need <- min(n_items - length(chosen_items), length(extras))
  if (need > 0) chosen_items <- c(chosen_items, sample(extras, need))
}

message(sprintf("  items: %d kept across %d classes x %d difficulty bins",
                length(chosen_items), n_classes, N_DIFF_BINS))

d <- d %>% filter(child_id %in% chosen_children, item %in% chosen_items)
message(sprintf("  subset: %d rows, %d children, %d items",
                nrow(d), length(unique(d$child_id)),
                length(unique(d$item))))

# Build admin keys and indices
d <- d %>%
  mutate(admin_key = paste(child_id, age, form, sep = "_"),
         aa = as.integer(factor(admin_key)),
         ii = as.integer(factor(child_id)),
         jj = as.integer(factor(item)),
         cc = as.integer(factor(lexical_category)))

admin_info <- d %>% distinct(aa, ii, age, admin_key) %>% arrange(aa)
# One row per item: when the same word appears across WG and WS with
# slightly different prob or lexical_category, keep the first occurrence
# so length(cc) == J.
word_info  <- d %>% group_by(jj) %>%
  summarise(item = first(item), prob = first(prob), cc = first(cc),
            .groups = "drop") %>%
  arrange(jj)
class_levels <- levels(factor(d$lexical_category))

I <- max(d$ii); A <- max(d$aa); J <- max(d$jj); C <- max(d$cc)
message(sprintf("  I=%d, A=%d, J=%d, C=%d, N=%d", I, A, J, C, nrow(d)))

prior_r <- load_input_rate_prior()

# a_0 is the reference age (where log_age = 0 and theta = xi). Set it to
# the dataset's median admin age so xi is the per-child posterior at the
# natural pivot of the data, not at an arbitrary global constant.
# Using a global a_0 = 20 across datasets caused the per-kid (xi, zeta)
# posterior ridge to tilt by an amount proportional to (median_age - 20),
# producing a spurious cross-language sign flip in rho_xi_zeta.
a0_dataset <- round(median(admin_info$age))
message(sprintf("  a0 (dataset median admin age) = %d", a0_dataset))

stan_data <- c(
  list(
    N = nrow(d),
    A = A, I = I, J = J, C = C,
    aa = d$aa, jj = d$jj,
    admin_to_child = admin_info$ii,
    cc = word_info$cc,
    y  = d$produces,
    admin_age = admin_info$age,
    log_p = log(word_info$prob),
    log_H = MODEL_CONSTANTS$log_H,
    a0    = a0_dataset,
    mu_r = prior_r$mu_r,
    sigma_r = prior_r$sigma_r
  ),
  DEFAULT_PRIORS
)

bundle <- list(
  stan_data = stan_data,
  admin_info = admin_info,
  word_info = word_info,
  class_levels = class_levels,
  df = d,
  language = language
)

saveRDS(bundle, file.path(PATHS$fits_dir, "long_subset_data.rds"))
cat(sprintf("\nSaved fits/long_subset_data.rds\n"))
