## Fit ONE glmer ladder model on ONE language's longitudinal data.
##
## Designed for one-fit-per-SLURM-array-task: single core, deterministic
## inputs, writes a single summary CSV + a fit RDS.
##
## Usage:
##   Rscript glmer_ladder/02_fit_one.R <lang_slug> <model_id>
##
## lang_slug: matches the file fits/glmer_ladder/data_<lang_slug>.rds
##            (e.g. "english_american_", "norwegian", "japanese")
## model_id : one of A, B_log, B_lin, C_log, C_lin, D_log, D_lin
##
## Output:
##   fits/glmer_ladder/fit_<lang_slug>_<model_id>.rds      -- the glmer fit
##   fits/glmer_ladder/summary_<lang_slug>_<model_id>.csv  -- one-row summary
##
## Models on cumulative nested ladder:
##   A     : produces ~ offset(log_age) + (1 | item)
##                                  -- pure unit accumulator (kappa=1)
##   B_log : produces ~ 1 + log_age + (1 | item)
##                                  -- free fixed kappa on log_age
##   B_lin : produces ~ 1 + age_c   + (1 | item)
##                                  -- free fixed slope on linear age
##   C_log : + (1 | child)          -- random kid intercept
##   C_lin : + (1 | child)
##   D_log : + (1 + log_age | child)-- random kid intercept + slope (M_best)
##   D_lin : + (1 + age_c   | child)

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(lme4); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
LANG_SLUG <- args[1]
MODEL_ID  <- args[2]

stopifnot(
  !is.na(LANG_SLUG), !is.na(MODEL_ID),
  MODEL_ID %in% c("A", "B_log", "B_lin", "C_log", "C_lin", "D_log", "D_lin")
)

LADDER_DIR <- file.path(PATHS$fits_dir, "glmer_ladder")
DATA_RDS   <- file.path(LADDER_DIR, sprintf("data_%s.rds", LANG_SLUG))
OUT_RDS    <- file.path(LADDER_DIR, sprintf("fit_%s_%s.rds", LANG_SLUG, MODEL_ID))
OUT_CSV    <- file.path(LADDER_DIR, sprintf("summary_%s_%s.csv", LANG_SLUG, MODEL_ID))
OUT_RANEF  <- file.path(LADDER_DIR, sprintf("ranef_%s_%s.csv", LANG_SLUG, MODEL_ID))

if (!file.exists(DATA_RDS)) stop(sprintf("missing data: %s", DATA_RDS))
data <- readRDS(DATA_RDS)
cat(sprintf("[%s / %s] loaded: %d kids, %d items, %d obs\n",
            LANG_SLUG, MODEL_ID,
            data$n_kids, data$n_items, data$n_obs))

## ---- Build modeling frame ----
## Use the median age across all admins as the centering pivot. This
## puts age_c roughly on the same numerical scale as log_age (both
## comparable in magnitude over the age window) and centers fixed
## intercept near 0.
a0 <- median(data$df$age, na.rm = TRUE)
df <- data$df |>
  mutate(
    age_c   = age - a0,
    log_age = log(age / a0),
    child   = factor(child_id),
    item    = factor(item),
    produces = as.integer(produces)
  ) |>
  filter(!is.na(produces))

cat(sprintf("Pivot a0 = %g mo; modeling N = %d rows\n", a0, nrow(df)))

## ---- Build glmer formula based on MODEL_ID ----
formula_for <- function(id) {
  switch(id,
    A     = produces ~ offset(log_age) + (1 | item),
    B_log = produces ~ 1 + log_age + (1 | item),
    B_lin = produces ~ 1 + age_c   + (1 | item),
    C_log = produces ~ 1 + log_age + (1 | child) + (1 | item),
    C_lin = produces ~ 1 + age_c   + (1 | child) + (1 | item),
    D_log = produces ~ 1 + log_age + (1 + log_age | child) + (1 | item),
    D_lin = produces ~ 1 + age_c   + (1 + age_c   | child) + (1 | item)
  )
}
form <- formula_for(MODEL_ID)
cat("Formula: "); print(form)

ctrl <- glmerControl(optimizer = "bobyqa",
                      optCtrl = list(maxfun = 200000))

t0 <- Sys.time()
fit <- glmer(form, data = df, family = binomial(),
              control = ctrl, nAGQ = 0)
dt_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Fit took %.1f min\n", dt_min))

## ---- Save the fit + a one-row summary ----
saveRDS(fit, OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))

vc <- VarCorr(fit)
fe <- fixef(fit)
summary_row <- tibble(
  language = data$language,
  lang_slug = LANG_SLUG,
  model = MODEL_ID,
  formula = paste(deparse(form), collapse = " "),
  a0 = a0,
  n_kids = data$n_kids,
  n_items = data$n_items,
  n_obs = data$n_obs,
  logLik = as.numeric(logLik(fit)),
  AIC = AIC(fit),
  BIC = BIC(fit),
  df = attr(logLik(fit), "df"),
  fixed_intercept = ifelse("(Intercept)" %in% names(fe), fe[["(Intercept)"]], NA_real_),
  fixed_age_slope = if ("log_age" %in% names(fe)) fe[["log_age"]]
                    else if ("age_c" %in% names(fe)) fe[["age_c"]] else NA_real_,
  sd_item = ifelse("item" %in% names(vc), attr(vc$item, "stddev")[1], NA_real_),
  sd_child_intercept = ifelse("child" %in% names(vc), attr(vc$child, "stddev")[1], NA_real_),
  sd_child_slope = if ("child" %in% names(vc) && length(attr(vc$child, "stddev")) >= 2)
                    attr(vc$child, "stddev")[2] else NA_real_,
  cor_child_int_slope = if ("child" %in% names(vc) && !is.null(attr(vc$child, "correlation"))
                            && nrow(attr(vc$child, "correlation")) >= 2)
                         attr(vc$child, "correlation")[1, 2] else NA_real_,
  fit_minutes = dt_min,
  fit_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write_csv(summary_row, OUT_CSV)
cat(sprintf("Wrote %s\n", OUT_CSV))
cat("\nSummary:\n")
print(summary_row)

## ---- Per-kid BLUPs (slim CSV; portable across machines) ----------
## Extract ranef(fit)$child for downstream demographic regressions
## without needing to pull back the full fit RDS. Only meaningful when
## there's a child random effect (C and D models).
re_child <- tryCatch(ranef(fit)$child, error = function(e) NULL)
if (!is.null(re_child)) {
  ranef_df <- tibble(
    child_id = rownames(re_child),
    intercept_blup = re_child[, 1],
    slope_blup     = if (ncol(re_child) >= 2) re_child[, 2] else NA_real_,
    language       = data$language,
    lang_slug      = LANG_SLUG,
    model          = MODEL_ID
  )
  write_csv(ranef_df, OUT_RANEF)
  cat(sprintf("Wrote %s (%d kids)\n", OUT_RANEF, nrow(ranef_df)))
} else {
  cat("No child random effects in this model — no ranef CSV written.\n")
}
