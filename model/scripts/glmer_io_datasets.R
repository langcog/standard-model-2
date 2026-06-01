## Per-IO-dataset glmer sanity check on kappa heterogeneity.
##
## The pooled hierarchical model gives kappa_study = 3.4 (BabyView),
## 6.6 (SEEDLingS), 8.2 (AM2018), 7.5 (FMW2013). The diagnostic spaghetti
## shows BabyView trajectories that look healthy and broadly similar to
## the other studies, so the question is: does a vanilla glmer recover
## the same per-study kappa pattern? If yes, heterogeneity is real;
## if no, the pooled IO model is doing something to BabyView specifically.
##
## Fits C_log = logistic glmer with random intercept per kid + per item,
## no random slope (fast — population-mean slope only). Slope coefficient
## on log(age/a0) is the population kappa per dataset.
##
## Output: fits/glmer_io_perstudy.rds + console table.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(lme4)
})

A0 <- 15  # pivot age (mo); choice doesn't change the slope, only the intercept

bundles <- list(
  BabyView  = readRDS(file.path(PATHS$fits_dir, "babyview_subset_data.rds")),
  SEEDLingS = readRDS(file.path(PATHS$fits_dir, "seedlings_subset_data.rds")),
  AM2018    = readRDS(file.path(PATHS$fits_dir, "io_am2018_subset_data.rds")),
  FMW2013   = readRDS(file.path(PATHS$fits_dir, "io_fmw2013_subset_data.rds"))
)

prep <- function(b) {
  df <- b$df
  df$la <- log(df$age / A0)
  df$kid  <- factor(df$ii)
  df$item <- factor(df$item)
  df
}

fit_one <- function(name, b) {
  d <- prep(b)
  cat(sprintf("[%s] N=%d  I=%d  J=%d  -- fitting C_log ...",
              name, nrow(d), nlevels(d$kid), nlevels(d$item)))
  t0 <- Sys.time()
  m <- glmer(produces ~ la + (1 | kid) + (1 | item),
             data = d, family = binomial,
             control = glmerControl(optimizer = "bobyqa",
                                    optCtrl = list(maxfun = 2e5)))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf(" done (%.1f min)\n", dt))
  m
}

results <- list()
for (n in names(bundles)) results[[n]] <- fit_one(n, bundles[[n]])

cat("\n=== C_log fixed-effect log(age/a0) slope per IO dataset ===\n")
tab <- data.frame(study = names(results), kappa_glmer = NA_real_, se = NA_real_,
                  sigma_intercept_kid = NA_real_, sigma_intercept_item = NA_real_)
for (i in seq_along(results)) {
  m <- results[[i]]
  fe <- fixef(m); se <- sqrt(diag(vcov(m)))
  vc <- VarCorr(m)
  tab$kappa_glmer[i] <- fe["la"]
  tab$se[i]          <- se["la"]
  tab$sigma_intercept_kid[i]  <- attr(vc$kid,  "stddev")["(Intercept)"]
  tab$sigma_intercept_item[i] <- attr(vc$item, "stddev")["(Intercept)"]
}
print(tab, digits = 3)

cat("\n=== compare to pooled Stan kappa_study ===\n")
cat("  BabyView   pooled  3.4   glmer ", round(tab$kappa_glmer[1], 2), "\n")
cat("  SEEDLingS  pooled  6.6   glmer ", round(tab$kappa_glmer[2], 2), "\n")
cat("  AM2018     pooled  8.2   glmer ", round(tab$kappa_glmer[3], 2), "\n")
cat("  FMW2013    pooled  7.5   glmer ", round(tab$kappa_glmer[4], 2), "\n")
cat("  (longitudinal Wordbank pooled: kappa_pop ~ 11.3)\n")

saveRDS(results, file.path(PATHS$fits_dir, "glmer_io_perstudy.rds"))
cat("\nSaved glmer_io_perstudy.rds\n")
