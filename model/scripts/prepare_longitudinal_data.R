## Prepare the longitudinal Wordbank WS data for Stan fitting.
##
## Usage:
##   Rscript model/scripts/prepare_longitudinal_data.R [language] [n_children] [n_items]
## Defaults: English, all longitudinal children, stratified subsample of 200 items.
##
## Reads:   model/fits/long_ws_items.rds (from pull_longitudinal.R)
## Writes:  model/fits/long_subset_data.rds (Stan-ready bundle)

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
language   <- if (length(args) >= 1) args[1] else "English (American)"
n_children <- as.integer(if (length(args) >= 2) args[2] else 600)
n_items    <- as.integer(if (length(args) >= 3) args[3] else 200)
SEED <- 20260421

message(sprintf("Preparing longitudinal subset: language=%s, children=%d, items=%d",
                language, n_children, n_items))

long <- readRDS(file.path(PATHS$fits_dir, "long_ws_items.rds"))
d <- long %>%
  filter(language == !!language,
         !is.na(prob), prob > 0,
         !is.na(produces))

message(sprintf("  input: %d rows, %d children, %d items",
                nrow(d), length(unique(d$child_id)),
                length(unique(d$item))))

# Stratify sample
set.seed(SEED)
# Children by age span (prioritize those with more data)
child_summary <- d %>%
  distinct(child_id, age) %>%
  group_by(child_id) %>%
  summarise(n_admin = n(), span = max(age) - min(age), .groups = "drop")

# Prefer kids with more admins and wider spans
chosen_children <- child_summary %>%
  mutate(priority = n_admin + span / 5) %>%
  arrange(desc(priority)) %>%
  head(n_children) %>%
  pull(child_id)

# Stratify items by lexical_category to match cross-sectional sampling
chosen_items <- d %>%
  distinct(item, lexical_category) %>%
  group_by(lexical_category) %>%
  slice_sample(n = max(1, floor(n_items / length(unique(d$lexical_category))))) %>%
  ungroup() %>%
  pull(item)
if (length(chosen_items) > n_items)
  chosen_items <- sample(chosen_items, n_items)

d <- d %>% filter(child_id %in% chosen_children, item %in% chosen_items)
message(sprintf("  subset: %d rows, %d children, %d items",
                nrow(d), length(unique(d$child_id)),
                length(unique(d$item))))

# Build admin keys and indices
d <- d %>%
  mutate(admin_key = paste(child_id, age, sep = "_"),
         aa = as.integer(factor(admin_key)),
         ii = as.integer(factor(child_id)),
         jj = as.integer(factor(item)),
         cc = as.integer(factor(lexical_category)))

admin_info <- d %>% distinct(aa, ii, age, admin_key) %>% arrange(aa)
word_info  <- d %>% distinct(jj, item, prob, cc) %>% arrange(jj)
class_levels <- levels(factor(d$lexical_category))

I <- max(d$ii); A <- max(d$aa); J <- max(d$jj); C <- max(d$cc)
message(sprintf("  I=%d, A=%d, J=%d, C=%d, N=%d", I, A, J, C, nrow(d)))

prior_r <- load_input_rate_prior()

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
    a0    = MODEL_CONSTANTS$a0,
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
cat(sprintf("\nSaved model/fits/long_subset_data.rds\n"))
