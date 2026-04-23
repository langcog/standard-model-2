## Fit log_irt_v2.stan to real Wordbank (English CDI:WS) + CHILDES data.
##
## Uses the existing preprocessed dataset from standard_model/scripts/data/.
## Starts with a subsample for a real-data shakedown; full-data fit
## controlled by SUBSET flag.

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(ggplot2)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

here <- function(...) file.path("/Users/mcfrank/Projects/standard_model_2", ...)

## ---- Config --------------------------------------------------------------
SUBSET      <- TRUE            # TRUE = fit small subsample first
N_CHILDREN  <- 500             # subsample size (ignored if SUBSET=FALSE)
N_WORDS     <- 200             # subsample size (ignored if SUBSET=FALSE)
N_CHAINS    <- 4
N_ITER      <- 1500
N_WARMUP    <- 750
OUT_TAG     <- if (SUBSET) "subset" else "full"

## ---- 1. Load and shape --------------------------------------------------
load(here("standard_model/scripts/data/engWS_preprocessed.Rdata"))
hr <- read.csv(here("standard_model/scripts/data/hourly_tokens_Sperry_HartRisley.csv"))

## External prior on log CDS input rate (tokens/hour)
hr_clean <- hr$adult_child_tokens_hr[!is.na(hr$adult_child_tokens_hr)]
MU_R    <- mean(log(hr_clean))            # 7.34
SIGMA_R <- sd(log(hr_clean))              # 0.53
cat(sprintf("External input prior: log_r ~ N(%.3f, %.3f^2) from n=%d samples\n",
            MU_R, SIGMA_R, length(hr_clean)))

## d_wf is (5,492 children x 680 items) in long format
## Columns: person, produces, age, item, lexical_class, prob
d <- d_wf %>%
  select(person, age, item, lexical_class, prob, produces) %>%
  filter(!is.na(produces), !is.na(prob), prob > 0) %>%
  mutate(produces = as.integer(produces))

## Subsample for shakedown
if (SUBSET) {
  set.seed(20250420)
  # stratified across age bins
  persons <- d %>% distinct(person, age) %>%
    mutate(age_bin = cut(age, breaks = seq(15, 31, 3))) %>%
    group_by(age_bin) %>%
    slice_sample(n = max(1, floor(N_CHILDREN / 5))) %>%
    ungroup() %>%
    pull(person)
  if (length(persons) > N_CHILDREN) persons <- sample(persons, N_CHILDREN)

  # stratified across lexical class
  items <- d %>% distinct(item, lexical_class) %>%
    group_by(lexical_class) %>%
    slice_sample(n = max(1, floor(N_WORDS / 5))) %>%
    ungroup() %>%
    pull(item)
  if (length(items) > N_WORDS) items <- sample(items, N_WORDS)

  d <- d %>% filter(person %in% persons, item %in% items)
}

cat(sprintf("After subsetting: %d children, %d items, %d obs (mean produces=%.3f)\n",
            length(unique(d$person)), length(unique(d$item)),
            nrow(d), mean(d$produces)))

## Build indices
d <- d %>%
  mutate(ii = as.integer(factor(person)),
         jj = as.integer(factor(item)),
         cc = as.integer(factor(lexical_class)))

I <- max(d$ii); J <- max(d$jj); C <- max(d$cc)
class_levels <- levels(factor(d$lexical_class))
cat("Lexical classes (coded):", paste(seq_along(class_levels), class_levels, sep="=", collapse=", "), "\n")

## Child-level age vector and word-level log_p / class vectors
child_age <- d %>% distinct(ii, age) %>% arrange(ii) %>% pull(age)
word_info <- d %>% distinct(jj, item, prob, cc) %>% arrange(jj)
log_p     <- log(word_info$prob)
cc        <- word_info$cc

stan_data <- list(
  N = nrow(d), I = I, J = J, C = C,
  ii = d$ii, jj = d$jj, cc = cc,
  y  = d$produces,
  age = child_age,
  log_p = log_p,
  log_H = log(365),          # 12 waking hrs * 30.44 days/mo
  a0    = 20,                # reference age
  mu_r = MU_R,
  sigma_r = SIGMA_R,
  mu_mu_c  = 8,
  sigma_mu_c = 3
)
saveRDS(list(data = stan_data,
             word_info = word_info,
             class_levels = class_levels,
             child_age = child_age),
        here(sprintf("model/wordbank_data_%s.rds", OUT_TAG)))

## ---- 2. Fit --------------------------------------------------------------
cat(sprintf("\nFitting log_irt_v2 on %s: I=%d, J=%d, N=%d\n",
            OUT_TAG, I, J, nrow(d)))
t0 <- Sys.time()
fit <- stan(
  file    = here("model/log_irt_v2.stan"),
  data    = stan_data,
  chains  = N_CHAINS, iter = N_ITER, warmup = N_WARMUP,
  seed    = 20250420,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
)
cat(sprintf("Total sampling time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---- 3. Quick diagnostics -----------------------------------------------
cat("\n--- Sampler diagnostics ---\n")
print(check_hmc_diagnostics(fit))

scalar_pars <- c("sigma_alpha", "s", "delta", "pi_alpha", "sigma_xi")
class_pars  <- c(paste0("mu_c[", 1:C, "]"),
                 paste0("tau_c[", 1:C, "]"))
print(summary(fit, pars = c(scalar_pars, class_pars))$summary[
         , c("mean", "2.5%", "50%", "97.5%", "n_eff", "Rhat")])

cat("\nLexical classes:\n")
for (c in seq_len(C)) {
  cat(sprintf("  %d = %s\n", c, class_levels[c]))
}

saveRDS(fit, here(sprintf("model/wordbank_fit_%s.rds", OUT_TAG)))
cat(sprintf("\nDone. Saved fit to model/wordbank_fit_%s.rds\n", OUT_TAG))
