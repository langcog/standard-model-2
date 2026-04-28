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

## Parse-check every R script in the project to catch bare-else and
## other Rscript-specific parse failures before they hit a SLURM queue.
cat("\nParse-checking all R scripts...\n")
all_R <- c(list.files("model/R", pattern = "\\.R$", full.names = TRUE),
           list.files("model/scripts", pattern = "\\.R$", full.names = TRUE))
ok_count <- 0
for (f in all_R) {
  res <- tryCatch(parse(file = f), error = function(e) e)
  if (inherits(res, "error")) {
    cat(sprintf("  FAIL  %s\n          %s\n", f, conditionMessage(res)))
  } else {
    ok_count <- ok_count + 1
  }
}
cat(sprintf("  parsed %d / %d scripts OK\n", ok_count, length(all_R)))
if (ok_count < length(all_R)) stop("Some scripts failed to parse.")

cat("\nSmoke test OK.\n")
