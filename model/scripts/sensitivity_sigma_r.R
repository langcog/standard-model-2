## Sensitivity analysis on the externally-pinned input-rate prior sigma_r.
##
## The baseline posterior claims pi_alpha ~ 0.93 (93% of child-level variance
## is learning efficiency, 7% input). That conclusion depends on trusting
## sigma_r = 0.53 (pooled Sperry/HR/WF). Real within-population variation in
## CDS may be wider (100k-1.5M words/month spans ~2.7 log units).
##
## This script refits a named variant with sigma_r in a sweep and reports
## how pi_alpha, sigma_alpha, sigma_xi move.
##
## Usage:
##   Rscript model/scripts/sensitivity_sigma_r.R <variant_name> [sigma_r_list]
## Example:
##   Rscript model/scripts/sensitivity_sigma_r.R 2pl 0.3,0.53,0.8,1.0

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript ... <variant> [sigma_r_list]")
nm <- args[1]
sr_values <- if (length(args) >= 2) {
  as.numeric(strsplit(args[2], ",")[[1]])
} else {
  c(0.3, 0.53, 0.8, 1.0)
}

bundle <- readRDS(file.path(PATHS$fits_dir, "subset_data.rds"))
overrides <- variant_hyperpriors(nm)

results <- list()
for (sr in sr_values) {
  tag <- sprintf("wordbank_%s_sigmaR_%.2f", nm, sr)
  cat(sprintf("\n===== variant=%s, sigma_r=%.3f =====\n", nm, sr))
  stan_data <- modifyList(bundle$stan_data, overrides)
  stan_data$sigma_r <- sr
  fit <- fit_variant(stan_data, tag = tag,
                     cfg = modifyList(DEFAULT_FIT_CONFIG,
                                      list(chains = 2, iter = 1500, warmup = 750)))
  pars <- c("sigma_alpha", "pi_alpha", "sigma_xi", "s", "delta")
  if (grepl("2pl", nm))    pars <- c(pars, "sigma_lambda")
  if (grepl("slopes", nm)) pars <- c(pars, "sigma_zeta")
  s <- summarize_fit(fit, pars = pars)
  s$sigma_r <- sr
  results[[length(results) + 1]] <- s
  cat("\n--- Scalars ---\n")
  print(s, digits = 3)
}

sens_tbl <- do.call(rbind, results)
saveRDS(sens_tbl,
        file.path(PATHS$fits_dir,
                  sprintf("sensitivity_sigma_r_%s.rds", nm)))

cat("\n===== Summary across sigma_r values =====\n")
pi_tbl <- sens_tbl[sens_tbl$param == "pi_alpha",
                   c("sigma_r", "median", "lo95", "hi95")]
print(pi_tbl, digits = 3)
alp_tbl <- sens_tbl[sens_tbl$param == "sigma_alpha",
                    c("sigma_r", "median", "lo95", "hi95")]
print(alp_tbl, digits = 3)

cat(sprintf("\nSaved: fits/sensitivity_sigma_r_%s.rds\n", nm))
