## 07_input_stability.R -- is the constant-input-rate assumption defensible?
##
## The SI converts age to cumulative input by assuming a fixed rate: a child's cumulative
## input at age t is t x 561K tokens. That is an assumption, and a reviewer asked what
## supports it. The strongest available evidence is our own longitudinal input data rather
## than a citation, because the published rate estimates are cross-sectional.
##
## Four datasets with repeated within-child input measurements:
##   babyview   head-cam transcripts, 22 children, ~261 recordings each, 6-37 months
##   SEEDLingS  LENA daylongs, 44 children, ~13 each, 6-54 months
##   AM2018     LENA, 66 children, ~2 each, 14-20 months
##   FMW2013    LENA, 89 children, ~2 each, 17-26 months
##
## The quantity that matters is the WITHIN-child slope of log input rate on age. A
## between-child slope would confound age with who is in the sample at each age, so we fit
## log_input ~ age + (1 | child_id) per dataset and report the age coefficient. Since the
## response is logged, that coefficient is a proportional change per month, and the useful
## summary is what it implies across each dataset's observed span: exp(slope x span) - 1.
##
## AGE WINDOW. Slopes are fitted over 8-36 months, the range across which the constant-rate
## assumption is actually applied. This is not a cosmetic restriction, and it changes one
## result, so it is reported both ways below. SEEDLingS is measured MONTHLY from 6-17
## months and then once more at 54 months, with no observation in between; fitting a
## linear-in-age trend over the full range therefore interpolates across a 37-month gap,
## and the resulting slope (+0.0035/mo, CI excluding zero, +18% over the span) is driven
## entirely by the contrast between the dense first-year block and that isolated
## preschool follow-up. Within its dense region the SEEDLingS slope is -0.0037
## [-0.0156, +0.0081]. Since 54 months is far outside the CDI range, the windowed fit is
## the one that speaks to the assumption; the full-range fit is kept in `stats_full`.
##
## Usage:  Rscript studies/bayes_long/07_input_stability.R
## Output: paper/cache/si_input_stability.rds
##         (points, fitted lines, stats [8-36 mo], stats_full [all ages])

suppressPackageStartupMessages({library(dplyr); library(here); library(lme4)})

LAB <- c(babyview = "BabyView (head-cam)", SEEDLingS = "SEEDLingS (LENA)",
         AM2018 = "AM2018 (LENA)", FMW2013 = "FMW2013 (LENA)")

d <- readr::read_csv(here("data", "harmonized", "input_level.csv"), show_col_types = FALSE) |>
  filter(paper_code %in% names(LAB), is.finite(log_input), age > 0) |>
  transmute(paper_code, child_id = as.character(child_id), age, log_input)

AGE_LO <- 8; AGE_HI <- 36                      # the range the assumption is applied over
win <- function(x, use_window) if (use_window) filter(x, age >= AGE_LO, age <= AGE_HI) else x

fit_one <- function(pc, use_window = TRUE) {
  dd <- win(filter(d, paper_code == pc), use_window)
  ## Random intercept per child: the fixed age term is then the within-child trend,
  ## not a mixture of within- and between-child variation.
  m  <- lmer(log_input ~ age + (1 | child_id), data = dd, REML = TRUE)
  b  <- unname(fixef(m)["age"]); se <- sqrt(diag(vcov(m)))[["age"]]
  span <- diff(range(dd$age))
  data.frame(
    paper_code = pc, lang = unname(LAB[pc]),
    n_obs = nrow(dd), n_kids = n_distinct(dd$child_id),
    age_lo = min(dd$age), age_hi = max(dd$age),
    slope = b, se = se, lo = b - 1.96 * se, hi = b + 1.96 * se,
    ## proportional change in rate implied across the observed age span
    pct_change = 100 * (exp(b * span) - 1),
    pct_lo = 100 * (exp((b - 1.96 * se) * span) - 1),
    pct_hi = 100 * (exp((b + 1.96 * se) * span) - 1),
    row.names = NULL)
}
stats      <- bind_rows(lapply(names(LAB), fit_one, use_window = TRUE)) |>
  mutate(lang = factor(lang, levels = unname(LAB)))
stats_full <- bind_rows(lapply(names(LAB), fit_one, use_window = FALSE)) |>
  mutate(lang = factor(lang, levels = unname(LAB)))

## Fitted within-child line over each dataset's own age range, with a confidence band on
## the FIXED effects. The band is what makes "non-significant" visible rather than merely
## asserted: where a horizontal line fits inside it, the data do not distinguish the
## dataset from a constant rate. `flat` is that constant-rate null -- the fitted value at
## the midpoint age, i.e. the same mean rate with the age term set to zero.
fitted_lines <- bind_rows(lapply(seq_len(nrow(stats)), function(i) {
  s <- stats[i, ]; dd <- win(filter(d, paper_code == s$paper_code), TRUE)
  m <- lmer(log_input ~ age + (1 | child_id), data = dd, REML = TRUE)
  ages <- seq(s$age_lo, s$age_hi, length.out = 60)
  X  <- cbind(1, ages)                       # design matrix for the fixed effects
  fit <- as.vector(X %*% fixef(m))
  se  <- sqrt(rowSums((X %*% as.matrix(vcov(m))) * X))   # diag(X V X') without the full matrix
  data.frame(lang = s$lang, age = ages, log_input = fit,
             lo = fit - 1.96 * se, hi = fit + 1.96 * se,
             flat = as.vector(c(1, mean(range(dd$age))) %*% fixef(m)))
}))

pts <- d |> mutate(lang = factor(unname(LAB[paper_code]), levels = unname(LAB)))

## ---- what a drift of this size would actually do -----------------------
## The assumption enters only through cumulative input. Constant rate gives cumulative
## proportional to t; a rate growing at beta per month gives (exp(beta t) - 1)/beta. The
## practical question is how far that moves a child along the log-token axis of the
## fractional-return figure, and the natural yardstick is the +/-1 SD between-child input
## band already drawn there, which spans a factor of 1.71 (0.536 log units).
CDI_LO <- AGE_LO; CDI_HI <- AGE_HI; BAND <- log(1.71)
## vectorised over b; the b -> 0 limit of (exp(bt)-1)/b is t
cum <- function(b, t) ifelse(abs(b) < 1e-12, t, (exp(b * t) - 1) / b)
stats <- stats |> mutate(
  shift_log = log(cum(slope, CDI_HI) / CDI_HI) - log(cum(slope, CDI_LO) / CDI_LO),
  shift_pct_of_band = 100 * abs(shift_log) / BAND)
cat("\nConsequence for cumulative input over the CDI range (", CDI_LO, "-", CDI_HI, " months):\n", sep = "")
for (i in seq_len(nrow(stats))) with(stats[i, ],
  cat(sprintf("  %-20s horizontal shift %.3f log units = %.0f%% of the +/-1 SD input band\n",
              lang, shift_log, shift_pct_of_band)))

saveRDS(list(points = pts, fitted = fitted_lines, stats = stats, stats_full = stats_full,
             age_window = c(AGE_LO, AGE_HI)),
        here("paper", "cache", "si_input_stability.rds"))

show <- function(x) print(x |> transmute(
  Dataset = lang, n_obs, n_kids,
  ages = sprintf("%.0f-%.0f", age_lo, age_hi),
  `slope/mo` = sprintf("%+.4f", slope),
  `95% CI` = sprintf("[%+.4f, %+.4f]", lo, hi),
  `change over span` = sprintf("%+.1f%% [%+.1f, %+.1f]", pct_change, pct_lo, pct_hi),
  flat = ifelse(lo < 0 & hi > 0, "yes", "NO")) |> as.data.frame(), row.names = FALSE)

cat("\n=== within-child slope, fitted over", AGE_LO, "-", AGE_HI, "months (REPORTED) ===\n")
show(stats)
cat("\n=== same fit over each dataset's full observed range (not reported) ===\n")
show(stats_full)
cat("\nThe two differ only for SEEDLingS, whose full range reaches a lone 54-month\n")
cat("follow-up with no observation between 17 and 54 months; see the header note.\n")
cat("\nwrote", here("paper", "cache", "si_input_stability.rds"), "\n")
