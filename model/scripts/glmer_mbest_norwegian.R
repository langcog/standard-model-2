## Fit glmer M_best on the NO longitudinal bundle so we have ranef BLUPs
## for ξ_i and ζ_i (efficiency level and change-in-efficiency).
##
## The Stan draws.rds files only saved global scalars; refitting via glmer
## is much faster and per the §28 equivalency check the BLUPs are point-
## for-point matches of Stan posterior medians.
##
## Same formula as glmer_mbest_comparison.R but on NO bundle.
##
## Expected runtime at N=2.2M, nAGQ=0: ~1 hr.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(lme4)
})

BUNDLE_PATH <- file.path(PATHS$fits_dir, "long_subset_data_nor.rds")
OUT_RDS     <- file.path(PATHS$fits_dir, "glmer_mbest_no.rds")

bundle <- readRDS(BUNDLE_PATH)
sd_b   <- bundle$stan_data
a0     <- sd_b$a0
cat(sprintf("NO bundle: I=%d, J=%d, N=%d, a0=%g\n",
            sd_b$I, sd_b$J, sd_b$N, a0))

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
