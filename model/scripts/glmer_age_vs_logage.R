## Does linear age earn its keep on top of log_age?
##
## Kachergis (2021) had β_a · a_i + β_t · log(p_j · H · a_i) — i.e.,
## BOTH linear age and log-age in the linear predictor. M_best dropped
## the linear-age term and kept only log-age. Question: with both terms
## free, does linear age carry meaningful weight, or is it ~0?
##
## Three glmers on the EN bundle, all binomial:
##   M_lin      : produces ~ 1 + age + (1 + age | child) + (1 | word)
##   M_log      : produces ~ 1 + log_age + (1 + log_age | child) + (1 | word)  [= M_best equivalent]
##   M_both     : produces ~ 1 + age + log_age + (1 + age + log_age | child) + (1 | word)
##
## Comparison: AIC, BIC, marginal coefficients, signs/sizes.
##
## Expected runtime: M_log already done (27 min); M_lin similar; M_both ~50 min.
## Use nAGQ=0 for speed.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(lme4)
})

bundle <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))
a0 <- bundle$stan_data$a0

df <- bundle$df |>
  mutate(
    age_c    = age - 18,                  # center at a0=18 to put linear and log on similar grids
    log_age  = log(age / a0),
    child    = factor(child_id),
    word     = factor(item),
    produces = as.integer(produces)
  ) |>
  filter(!is.na(produces))

ctrl <- glmerControl(optimizer = "bobyqa",
                      optCtrl = list(maxfun = 200000))

cat("=== M_lin: linear age only ===\n")
t0 <- Sys.time()
m_lin <- glmer(produces ~ 1 + age_c + (1 + age_c | child) + (1 | word),
                data = df, family = binomial(), control = ctrl, nAGQ = 0)
cat(sprintf("M_lin fit %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(m_lin, file.path(PATHS$fits_dir, "glmer_mbest_en_linage.rds"))

cat("\n=== M_both: linear + log age ===\n")
t0 <- Sys.time()
m_both <- glmer(produces ~ 1 + age_c + log_age + (1 + age_c + log_age | child) + (1 | word),
                 data = df, family = binomial(), control = ctrl, nAGQ = 0)
cat(sprintf("M_both fit %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(m_both, file.path(PATHS$fits_dir, "glmer_mbest_en_both.rds"))

m_log <- readRDS(file.path(PATHS$fits_dir, "glmer_mbest_en.rds"))

cat("\n========= COMPARISON =========\n")
ic <- tibble::tribble(
  ~model,    ~AIC,        ~BIC,        ~logLik,
  "M_lin",   AIC(m_lin),  BIC(m_lin),  as.numeric(logLik(m_lin)),
  "M_log",   AIC(m_log),  BIC(m_log),  as.numeric(logLik(m_log)),
  "M_both",  AIC(m_both), BIC(m_both), as.numeric(logLik(m_both))
)
print(ic)

cat("\nM_both fixed effects:\n")
print(summary(m_both)$coefficients)
cat("\nM_both random-effect variances:\n")
print(VarCorr(m_both))

cat("\nM_lin fixed effects:\n")
print(summary(m_lin)$coefficients)
cat("\nM_lin random-effect variances:\n")
print(VarCorr(m_lin))
