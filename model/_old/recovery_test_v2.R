## Re-fit v2 (reparameterized) model on the SAME sim used for v1.
## Run after recovery_test.R to compare recovery of sigma_alpha and pi_alpha
## between the two parameterizations.

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(ggplot2)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

here <- function(...) file.path("/Users/mcfrank/Projects/standard_model_2/model", ...)
source(here("simulate.R"))

sim <- readRDS(here("recovery_sim.rds"))
stan_data <- build_stan_data(sim)

cat(sprintf("Fitting v2 on same sim: I=%d, J=%d, N=%d\n",
            sim$constants$I, sim$constants$J, nrow(sim$obs)))

fit2 <- stan(
  file   = here("log_irt_v2.stan"),
  data   = stan_data,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 123,
  control = list(adapt_delta = 0.9)
)

cat("\n--- v2 sampler diagnostics ---\n")
print(check_hmc_diagnostics(fit2))

scalar_pars <- c("sigma_alpha", "s", "delta", "pi_alpha")
class_pars  <- c(paste0("mu_c[", 1:3, "]"),
                 paste0("tau_c[", 1:3, "]"))
s2 <- summary(fit2, pars = c(scalar_pars, class_pars))$summary
print(s2[, c("mean", "2.5%", "50%", "97.5%", "n_eff", "Rhat")])

## Side-by-side comparison with v1
fit1 <- readRDS(here("recovery_fit.rds"))
s1 <- summary(fit1, pars = c(scalar_pars, class_pars))$summary

truth_tbl <- data.frame(
  param = c(scalar_pars, class_pars),
  truth = c(sim$true$sigma_alpha, sim$true$s, sim$true$delta,
            sim$true$sigma_alpha^2 / (sim$true$sigma_alpha^2 + sim$true$sigma_r^2),
            sim$true$mu_c, sim$true$tau_c),
  v1_median = s1[, "50%"],
  v1_lo = s1[, "2.5%"], v1_hi = s1[, "97.5%"],
  v2_median = s2[, "50%"],
  v2_lo = s2[, "2.5%"], v2_hi = s2[, "97.5%"]
)
truth_tbl$v1_ci <- truth_tbl$truth >= truth_tbl$v1_lo & truth_tbl$truth <= truth_tbl$v1_hi
truth_tbl$v2_ci <- truth_tbl$truth >= truth_tbl$v2_lo & truth_tbl$truth <= truth_tbl$v2_hi

cat("\n--- v1 vs v2 recovery (scalars) ---\n")
print(truth_tbl)

saveRDS(fit2, here("recovery_fit_v2.rds"))
cat("\nDone.\n")
