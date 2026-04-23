## Smoke test: verifies that config + helpers load cleanly and core
## functions are callable. Does NOT fit anything. Safe to run anytime.
##
## Usage: Rscript model/scripts/00_smoke.R

source("model/R/config.R")
source("model/R/helpers.R")

cat("PATHS:\n")
for (k in names(PATHS)) cat(sprintf("  %-10s = %s\n", k, PATHS[[k]]))

cat("\nStan model exists:", file.exists(PATHS$stan_model), "\n")
cat("Wordbank data exists:", file.exists(PATHS$wordbank), "\n")
cat("Input-rate data exists:", file.exists(PATHS$input_rate), "\n")

prior_r <- load_input_rate_prior()
cat(sprintf("\nInput prior: log_r ~ N(%.3f, %.3f^2) from n=%d\n",
            prior_r$mu_r, prior_r$sigma_r, prior_r$n))

for (v in c("baseline", "fix_delta", "fix_s", "both_fixed",
            "2pl", "2pl_fix_delta", "slopes", "2pl_slopes")) {
  overrides <- variant_hyperpriors(v)
  cat(sprintf("\nvariant=%s overrides:\n", v))
  if (length(overrides) == 0) cat("  (none)\n")
  else for (k in names(overrides))
    cat(sprintf("  %-18s = %s\n", k, overrides[[k]]))
}

cat("\nSmoke test OK.\n")
