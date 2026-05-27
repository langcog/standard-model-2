## glmer sanity check: fit the M_best-equivalent mixed-effects logistic
## regression on the EN longitudinal bundle and compare to the Stan
## posteriors.
##
## Stan M_best linear predictor (s = 0, s_i = 0):
##   eta = xi + log_H + (1 + delta + zeta) * log(age/a0) - delta_j
##
## In glmer notation:
##   glmer(produces ~ 1 + log_age + (1 + log_age | child) + (1 | word),
##         family = binomial())
##
## Expected mapping (Stan posterior medians from long_no_freq_slopes):
##   fixed log_age slope         ≈ kappa_pop = 1 + delta = 11.31
##   sd(Intercept | child)       ≈ sigma_xi  = sqrt(sigma_r^2 + sigma_alpha^2) = 1.89
##   sd(log_age   | child)       ≈ sigma_zeta = 3.81
##   cor((Int, log_age) | child) ≈ rho_xi_zeta = -0.09
##   sd(Intercept | word)        ≈ tau_delta (TBD from Stan fit)
##
## The fixed intercept is harder to map cleanly because it absorbs
##   log_H + mu_r - mu_delta + sum-to-zero centering shifts.
##   Stan: mu_r = 7.34, log_H = 5.90, mu_delta ≈ ?
##
## Output:
##   outputs/glmer_mbest_comparison.csv  -- side-by-side table
##   fits/glmer_mbest_en.rds             -- the fitted glmer object
##   console-printed pretty comparison

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(lme4); library(broom.mixed); library(readr); library(knitr)
})

BUNDLE_PATH <- file.path(PATHS$fits_dir, "long_subset_data.rds")
STAN_SUMM   <- file.path(PATHS$fits_dir, "summaries", "long_no_freq_slopes.summary.rds")
OUT_RDS     <- file.path(PATHS$fits_dir, "glmer_mbest_en.rds")
OUT_CSV     <- file.path(PATHS$outputs_dir, "glmer_mbest_comparison.csv")

bundle <- readRDS(BUNDLE_PATH)
sd_b   <- bundle$stan_data
a0     <- sd_b$a0
log_H  <- sd_b$log_H
cat(sprintf("EN bundle: I=%d, J=%d, N=%d, a0=%g, log_H=%.3f\n",
            sd_b$I, sd_b$J, sd_b$N, a0, log_H))

## ---- Build the per-observation data frame ----------------------
## bundle$df already has one row per (admin, item) with age, item,
## child_id, produces. Stan treats ~3% duplicate rows as independent
## observations; we leave them in so the comparison is apples-to-apples.
df <- bundle$df |>
  mutate(
    log_age  = log(age / a0),
    child    = factor(child_id),
    word     = factor(item),
    produces = as.integer(produces)
  ) |>
  filter(!is.na(produces))
cat(sprintf("Data: %d rows, %d kids, %d words\n",
            nrow(df), nlevels(df$child), nlevels(df$word)))

## ---- Fit glmer --------------------------------------------------
## bobyqa is the default; nAGQ = 0 trades a bit of statistical accuracy
## for substantial speed at N=1.1M. We're after a sanity check, not a
## publication-grade fit -- the σ estimates will be very close to nAGQ=1
## but the runtime is order-of-magnitude better.
##
## If you want the slower-but-more-accurate version, set nAGQ = 1.
## Expected runtime: nAGQ=0 ~15 min; nAGQ=1 ~1-2 hr at N=1.1M.

t0 <- Sys.time()
m <- glmer(
  produces ~ 1 + log_age + (1 + log_age | child) + (1 | word),
  data    = df,
  family  = binomial(),
  control = glmerControl(optimizer = "bobyqa",
                          optCtrl   = list(maxfun = 200000)),
  nAGQ    = 0
)
dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nglmer fit took %.1f min\n", dt))

saveRDS(m, OUT_RDS)
cat(sprintf("Saved %s\n", OUT_RDS))

cat("\n=== glmer summary ===\n")
print(summary(m), correlation = FALSE)

## ---- Extract glmer estimates -----------------------------------
fixed <- broom.mixed::tidy(m, effects = "fixed")
ranef <- broom.mixed::tidy(m, effects = "ran_pars")
cat("\n=== Fixed effects ===\n"); print(fixed)
cat("\n=== Random-effect parameters ===\n"); print(ranef)

get_fix <- function(term) fixed$estimate[fixed$term == term]
get_ran <- function(group, term) {
  ranef$estimate[ranef$group == group & ranef$term == term]
}

glmer_results <- list(
  fixed_intercept = get_fix("(Intercept)"),
  fixed_log_age   = get_fix("log_age"),
  sd_int_child    = get_ran("child", "sd__(Intercept)"),
  sd_logage_child = get_ran("child", "sd__log_age"),
  cor_child       = get_ran("child", "cor__(Intercept).log_age"),
  sd_int_word     = get_ran("word",  "sd__(Intercept)")
)

## ---- Pull Stan M_best posteriors -------------------------------
stan <- readRDS(STAN_SUMM)
get_stan <- function(p) {
  r <- stan[stan$variable == p, ]
  list(mean = r$mean[1], q025 = r$q025[1], q975 = r$q975[1])
}
sigma_alpha <- get_stan("sigma_alpha")
sigma_xi    <- get_stan("sigma_xi")
sigma_zeta  <- get_stan("sigma_zeta")
rho_xi_zeta <- get_stan("rho_xi_zeta")
delta_post  <- get_stan("delta")
mu_r        <- sd_b$mu_r
sigma_r     <- sd_b$sigma_r

## ---- Side-by-side comparison table ------------------------------
comp <- tibble::tribble(
  ~quantity,              ~glmer_estimate,           ~stan_mean,           ~stan_95ci,                         ~notes,
  "Fixed log_age slope",  sprintf("%.3f", glmer_results$fixed_log_age),
                          sprintf("%.3f", 1 + delta_post$mean),
                          sprintf("[%.2f, %.2f]", 1 + delta_post$q025, 1 + delta_post$q975),
                          "Should match 1 + delta = kappa_pop",
  "sd(Intercept | child)", sprintf("%.3f", glmer_results$sd_int_child),
                          sprintf("%.3f", sigma_xi$mean),
                          sprintf("[%.2f, %.2f]", sigma_xi$q025, sigma_xi$q975),
                          "Should match sigma_xi = sqrt(sigma_r^2 + sigma_alpha^2)",
  "sd(log_age | child)",  sprintf("%.3f", glmer_results$sd_logage_child),
                          sprintf("%.3f", sigma_zeta$mean),
                          sprintf("[%.2f, %.2f]", sigma_zeta$q025, sigma_zeta$q975),
                          "Should match sigma_zeta",
  "cor((Int, log_age) | child)", sprintf("%+.3f", glmer_results$cor_child),
                          sprintf("%+.3f", rho_xi_zeta$mean),
                          sprintf("[%+.2f, %+.2f]", rho_xi_zeta$q025, rho_xi_zeta$q975),
                          "Should match rho_xi_zeta",
  "sd(Intercept | word)",  sprintf("%.3f", glmer_results$sd_int_word),
                          "(not reported)",
                          "—",
                          "Stan reports tau_delta as a hyperprior; not in scalar summary"
)
cat("\n=== Stan vs glmer comparison ===\n")
print(knitr::kable(comp, format = "pipe", align = "lrrll"))

write_csv(comp, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))

## ---- Implied pi_alpha from glmer -------------------------------
## In glmer-land, sigma_xi is data-identified directly. If we BORROW
## sigma_r = 0.534 from the external pin, we get an implied
## sigma_alpha^2 = sigma_xi^2 - sigma_r^2 and pi_alpha = ratio.
sd_int_child <- glmer_results$sd_int_child
sigma_alpha_implied <- sqrt(pmax(sd_int_child^2 - sigma_r^2, 0))
pi_alpha_implied    <- sigma_alpha_implied^2 / sd_int_child^2

cat(sprintf(
"\n=== Implied pi_alpha from glmer + external sigma_r pin ===
  sd(Intercept | child)   (glmer)      = %.3f
  sigma_r                 (external)   = %.3f
  implied sigma_alpha     (residual)   = %.3f
  implied pi_alpha                     = %.3f
  Stan posterior:
    sigma_alpha = %.3f [%.2f, %.2f]
    pi_alpha    = %.3f
",
  sd_int_child, sigma_r, sigma_alpha_implied, pi_alpha_implied,
  sigma_alpha$mean, sigma_alpha$q025, sigma_alpha$q975,
  sigma_alpha$mean^2 / (sigma_alpha$mean^2 + sigma_r^2)
))
