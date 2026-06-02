## Per-IO-dataset glmer D_log fit (random slope per kid).
##
## Apples-to-apples comparison to Stan's kappa_pop: same random-slope
## structure. The fixed-effect log(age/a0) coefficient under D_log is
## the proper benchmark for whether Stan's kappa_study (3.4-8.2) is
## low because of model bias OR because C_log was overstating the
## population slope (Jensen's-inequality-style bias when forcing
## random-slope data into a fixed-slope model).
##
## Output: fits/glmer_io_perstudy_Dlog.rds + console table.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(lme4)
})

A0 <- 15

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
  cat(sprintf("[%s] N=%d  I=%d  J=%d  -- fitting D_log ...",
              name, nrow(d), nlevels(d$kid), nlevels(d$item)))
  t0 <- Sys.time()
  m <- glmer(produces ~ la + (la | kid) + (1 | item),
             data = d, family = binomial,
             control = glmerControl(optimizer = "bobyqa",
                                    optCtrl = list(maxfun = 5e5)))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf(" done (%.1f min)\n", dt))
  m
}

results <- list()
for (n in names(bundles)) results[[n]] <- fit_one(n, bundles[[n]])

cat("\n=== D_log per IO dataset (random slope per kid) ===\n")
tab <- data.frame(study = names(results),
                  kappa_D     = NA_real_, se_D        = NA_real_,
                  sigma_kid_slope = NA_real_, sigma_kid_int = NA_real_,
                  sigma_item_int = NA_real_)
for (i in seq_along(results)) {
  m <- results[[i]]
  fe <- fixef(m); se <- sqrt(diag(vcov(m)))
  vc <- VarCorr(m)
  tab$kappa_D[i]          <- fe["la"]
  tab$se_D[i]             <- se["la"]
  tab$sigma_kid_slope[i]  <- attr(vc$kid, "stddev")["la"]
  tab$sigma_kid_int[i]    <- attr(vc$kid, "stddev")["(Intercept)"]
  tab$sigma_item_int[i]   <- attr(vc$item, "stddev")["(Intercept)"]
}
print(tab, digits = 3)

cat("\n=== three-way comparison ===\n")
cat("  study      C_log kappa   D_log kappa   Stan kappa_study   Stan sigma_zeta\n")
glmer_C <- c(BabyView = 9.81, SEEDLingS = 11.67, AM2018 = 10.42, FMW2013 = 9.61)
stan_k  <- c(BabyView = 3.40, SEEDLingS = 6.59,  AM2018 = 8.15,  FMW2013 = 7.48)
for (n in tab$study) {
  cat(sprintf("  %-10s %9.2f   %9.2f   %14.2f   %14.2f\n",
              n, glmer_C[n], tab$kappa_D[tab$study==n],
              stan_k[n], tab$sigma_kid_slope[tab$study==n]))
}
cat("  (Stan baseline sigma_zeta = 5.22 across all)\n")

saveRDS(results, file.path(PATHS$fits_dir, "glmer_io_perstudy_Dlog.rds"))
cat("\nSaved fits/glmer_io_perstudy_Dlog.rds\n")
