## Prepare Norwegian longitudinal WS data for Stan fitting.
##
## Usage:   Rscript model/scripts/prepare_longitudinal_norwegian.R [n_children] [n_items]
## Writes:  model/fits/long_subset_data_nor.rds

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
n_children <- as.integer(if (length(args) >= 1) args[1] else 200)
n_items    <- as.integer(if (length(args) >= 2) args[2] else 200)
N_AGE_BINS  <- 4
N_DIFF_BINS <- 4
MIN_ADMINS  <- 2
SEED <- 20260428

message(sprintf("Preparing Norwegian longitudinal subset: children=%d, items=%d",
                n_children, n_items))

long <- readRDS(file.path(PATHS$fits_dir, "long_items.rds"))
freq <- readRDS(file.path(PATHS$fits_dir, "norwegian_word_freq.rds"))

# Filter to Norwegian; drop pre-existing (English-derived) `prob` column
d <- long %>%
  filter(language == "Norwegian",
         !is.na(produces)) %>%
  select(-any_of("prob"))

message(sprintf("  input: %d rows, %d children, %d items",
                nrow(d), length(unique(d$child_id)),
                length(unique(d$item))))

# Join with Norwegian CHILDES frequencies.
# Wordbank item_definition may contain parenthetical qualifiers (e.g.
# "kjole (mat)"); strip them for the join.
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

# For unmatched items, use minimum freq (small but nonzero)
if (sum(is.na(d$prob)) > 0) {
  min_prob <- min(d$prob, na.rm = TRUE)
  d <- d %>% mutate(prob = dplyr::coalesce(prob, min_prob))
}

# Drop items without freq (prob = 0 / NA) before stratification
d <- d %>% filter(prob > 0)

## Stratified sampling: kids x age, items x (class x difficulty quartile).
set.seed(SEED)

# --- Children ---
child_summary <- d %>%
  distinct(child_id, age) %>%
  group_by(child_id) %>%
  summarise(n_admin = n(),
            median_age = median(age),
            span = max(age) - min(age),
            .groups = "drop") %>%
  filter(n_admin >= MIN_ADMINS) %>%
  mutate(age_bin = ntile(median_age, N_AGE_BINS))

per_bin_kids <- max(1, floor(n_children / N_AGE_BINS))
chosen_children <- child_summary %>%
  group_by(age_bin) %>%
  arrange(desc(n_admin), desc(span)) %>%
  slice_head(n = per_bin_kids) %>%
  ungroup() %>%
  pull(child_id)
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
# Force one row per jj: if an item appears with multiple lexical_category
# labels in the raw data, take the first.
word_info  <- d %>% group_by(jj) %>%
  summarise(item = first(item), prob = first(prob), cc = first(cc),
            .groups = "drop") %>%
  arrange(jj)
class_levels <- levels(factor(d$lexical_category))

I <- max(d$ii); A <- max(d$aa); J <- max(d$jj); C <- max(d$cc)
message(sprintf("  I=%d, A=%d, J=%d, C=%d, N=%d", I, A, J, C, nrow(d)))

# Norwegian input-rate prior. We only have English Sperry/HR/WF; reuse
# it as the best available prior for Norwegian as well. (Note: Norwegian
# parenting styles / CDS rates may differ; this is a modeling choice to
# revisit.)
prior_r <- load_input_rate_prior()

# a_0 = dataset median admin age; see prepare_longitudinal_data.R for rationale.
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
  language = "Norwegian"
)

saveRDS(bundle, file.path(PATHS$fits_dir, "long_subset_data_nor.rds"))
cat("\nSaved model/fits/long_subset_data_nor.rds\n")
