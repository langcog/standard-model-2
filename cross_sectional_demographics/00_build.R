## Build the cross-sectional demographic analysis.
##
## For each language: pull one admin per child (monolingual-TD), score
## item-level production, and fit two Rasch GLMMs with RAW log(age/a0):
##   produces ~ predictor * log(age/a0) + (1|item) + (1|child)
## where predictor is sex (Male vs Female) or maternal ed (z-scored years).
## The predictor MAIN effect = effect on EFFICIENCY (intercept at a0);
## the predictor:log_age interaction = effect on ACCELERATION (slope).
## Units match the longitudinal raw-BLUP regressions (paper fig-demographics),
## so cross-sectional and longitudinal estimates are directly comparable.
##
## Caching (all under cache/, frames/ + fits/ are gitignored, regenerable):
##   cache/frames/<slug>.rds   per-language item-level modeling frame
##   cache/fits/<slug>_<k>.rds  per-(language,predictor) fit summary
##   cache/fits.rds            COMMITTED: xsec fits + meta + longitudinal
##   cache/scatter.rds         COMMITTED: per-child proportion for scatters
##
## Run: Rscript cross_sectional_demographics/00_build.R
## Re-run is cheap (frames + fits cached); delete a cache file to recompute.

suppressPackageStartupMessages({
  library(here); library(wordbankr); library(dplyr); library(tidyr)
  library(lme4); library(metafor)
})

DIR     <- here("cross_sectional_demographics")
FRAMES  <- file.path(DIR, "cache", "frames")
FITDIR  <- file.path(DIR, "cache", "fits")
dir.create(FRAMES, recursive = TRUE, showWarnings = FALSE)
dir.create(FITDIR, recursive = TRUE, showWarnings = FALSE)

ED_YEARS <- c("None"=0, "Primary"=6, "Some Secondary"=9, "Secondary"=12,
              "Some College"=14, "College"=16, "Some Graduate"=18, "Graduate"=20)
CLINICAL  <- c("Edgin", "Byers")
FORM_KEEP <- c("WG","WGProd","WGProdShort","WGShort","WS","WSShort")
slug <- function(s) gsub("[^a-z0-9]+", "_", tolower(s))

## Monolingual-TD languages with sex >= 200. English (British) excluded
## (TEDS twins, short form); French (Quebecois) excluded (bilingual).
LANGS <- c("English (American)","Norwegian","Danish","Portuguese (European)","Turkish",
  "Mandarin (Taiwanese)","Korean","Spanish (Mexican)","Russian","Slovak","Dutch",
  "English (Australian)","Catalan","Italian","Swedish","French (French)","Estonian",
  "Cantonese","German","Mandarin (Beijing)","Spanish (European)","Japanese",
  "Spanish (Argentinian)","Latvian","Croatian","Hebrew","Czech","Arabic (Saudi)",
  "Hungarian","Finnish","Kigiriama")

## ---- per-language item frame (one admin per child, linked by data_id) ----
build_frame <- function(L, n_sub = 1200, seed = 2026) {
  cf <- file.path(FRAMES, paste0(slug(L), ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  ad <- get_administration_data(language = L, include_demographic_info = TRUE) |>
    filter(!grepl("Bilingual", dataset_origin_name, ignore.case = TRUE),
           !(dataset_name %in% CLINICAL),
           form %in% FORM_KEEP, !is.na(age), !is.na(sex)) |>
    mutate(ed_years = unname(ED_YEARS[as.character(caregiver_education)]))
  set.seed(seed)
  ad1 <- ad |> group_by(child_id) |> slice_sample(n = 1) |> ungroup()
  keep <- sample(unique(ad1$child_id), min(n_sub, n_distinct(ad1$child_id)))
  ad1 <- ad1 |> filter(child_id %in% keep) |>
    select(data_id, child_id, age, form, sex, ed_years)
  prod <- bind_rows(lapply(intersect(FORM_KEEP, unique(ad1$form)), function(f) {
    d <- tryCatch(get_instrument_data(language = L, form = f,
                                      administration_info = TRUE, item_info = TRUE),
                  error = function(e) NULL)
    if (is.null(d)) return(NULL)
    out <- d |> filter(data_id %in% ad1$data_id, item_kind == "word") |>  # link by ADMIN
      transmute(data_id, item = item_definition,
                produces = as.integer(!is.na(value) & tolower(value) == "produces"))
    rm(d); gc(); out
  })) |> inner_join(ad1, by = "data_id") |> filter(!is.na(item), !is.na(produces))
  ik <- prod |> count(item) |> filter(n >= 100) |> pull(item)
  a0 <- median(ad1$age)
  fr <- prod |> filter(item %in% ik) |>
    mutate(la = log(age / a0), language = L, a0 = a0) |>
    select(language, a0, child_id, item, age, la, sex, ed_years, produces)
  saveRDS(fr, cf); fr
}

## ---- one predictor fit ----
fit_one <- function(fr, kind, L) {
  ff <- file.path(FITDIR, sprintf("%s_%s.rds", slug(L), kind))
  if (file.exists(ff)) return(readRDS(ff))
  if (kind == "sex") {
    d <- fr |> filter(!is.na(sex)) |> mutate(p = factor(sex, levels = c("Female","Male")))
    if (n_distinct(d$p) < 2 || nrow(d) < 1000) return(NULL)
    pn <- "pMale"
  } else {
    d <- fr |> filter(!is.na(ed_years))
    if (nrow(d) < 1000 || n_distinct(d$ed_years) < 2) return(NULL)
    d <- d |> mutate(p = as.numeric(scale(ed_years))); pn <- "p"
  }
  d$item <- factor(d$item); d$child_id <- factor(d$child_id)
  m <- tryCatch(glmer(produces ~ la * p + (1|item) + (1|child_id), data = d,
                      family = binomial(),
                      control = glmerControl(optimizer = "bobyqa"), nAGQ = 0),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  fe <- fixef(m); se <- sqrt(diag(as.matrix(vcov(m))))
  out <- data.frame(language = L, predictor = kind, n_kids = n_distinct(d$child_id),
                    a0 = fr$a0[1],
                    eff = unname(fe[pn]), eff_se = unname(se[pn]),
                    acc = unname(fe[paste0("la:", pn)]),
                    acc_se = unname(se[paste0("la:", pn)]))
  saveRDS(out, ff); out
}

## ---- run: frames, fits, scatter data ----
fits <- list(); scatter <- list()
for (L in LANGS) {
  cat(sprintf("== %s ==\n", L))
  fr <- tryCatch(build_frame(L), error = function(e) { cat(" pull-err:", conditionMessage(e), "\n"); NULL })
  if (is.null(fr)) next
  for (k in c("sex", "matEd")) {
    r <- fit_one(fr, k, L); if (!is.null(r)) fits[[paste(L, k)]] <- r
  }
  nI <- n_distinct(fr$item)
  scatter[[L]] <- fr |> group_by(child_id) |>
    summarise(age = first(age), sex = first(sex), ed_years = first(ed_years),
              prop = sum(produces) / nI, .groups = "drop") |>
    mutate(language = L, n_items = nI)
}
xsec <- do.call(rbind, fits)
scatter_df <- do.call(rbind, scatter)

## ---- random-effects meta-analysis per predictor x component ----
as_meta <- function(b, s) {
  ok <- is.finite(b) & is.finite(s) & s > 0
  if (sum(ok) < 2) return(data.frame(estimate=NA, ci.lb=NA, ci.ub=NA, k=sum(ok), pval=NA, tau2=NA))
  m <- tryCatch(rma(yi = b[ok], sei = s[ok], method = "REML",
                    control = list(stepadj = 0.5, maxiter = 1000)),
                error = function(e) NULL)
  if (is.null(m)) return(data.frame(estimate=NA, ci.lb=NA, ci.ub=NA, k=sum(ok), pval=NA, tau2=NA))
  data.frame(estimate = as.numeric(m$beta), ci.lb = m$ci.lb, ci.ub = m$ci.ub,
             k = m$k, pval = m$pval, tau2 = m$tau2)
}
meta <- bind_rows(lapply(c("sex","matEd"), function(p) {
  s <- xsec |> filter(predictor == p)
  bind_rows(cbind(predictor = p, component = "efficiency",   as_meta(s$eff, s$eff_se)),
            cbind(predictor = p, component = "acceleration", as_meta(s$acc, s$acc_se)))
}))

## ---- longitudinal raw-BLUP effects (same units) for the combined figure ----
blups <- readRDS(here("paper", "cache", "blups_demographics.rds"))
long <- bind_rows(lapply(unique(blups$lang_slug), function(u) {
  d <- blups |> filter(lang_slug == u)
  rows <- list()
  ds <- d |> filter(!is.na(sex), sex %in% c("Female","Male"))
  if (nrow(ds) >= 30 && n_distinct(ds$sex) == 2) {
    mx <- summary(lm(xi ~ sex, ds))$coef; mz <- summary(lm(zeta ~ sex, ds))$coef
    rows$sex <- data.frame(language = u, predictor = "sex",
      eff = mx["sexMale","Estimate"], eff_se = mx["sexMale","Std. Error"],
      acc = mz["sexMale","Estimate"], acc_se = mz["sexMale","Std. Error"])
  }
  dm <- d |> mutate(ed = unname(ED_YEARS[as.character(caregiver_ed)])) |> filter(!is.na(ed))
  if (nrow(dm) >= 30 && n_distinct(dm$ed) >= 2) {
    mx <- summary(lm(xi ~ scale(ed), dm))$coef; mz <- summary(lm(zeta ~ scale(ed), dm))$coef
    rows$matEd <- data.frame(language = u, predictor = "matEd",
      eff = mx["scale(ed)","Estimate"], eff_se = mx["scale(ed)","Std. Error"],
      acc = mz["scale(ed)","Estimate"], acc_se = mz["scale(ed)","Std. Error"])
  }
  bind_rows(rows)
}))
long_meta <- bind_rows(lapply(c("sex","matEd"), function(p) {
  s <- long |> filter(predictor == p)
  bind_rows(cbind(predictor = p, component = "efficiency",   as_meta(s$eff, s$eff_se)),
            cbind(predictor = p, component = "acceleration", as_meta(s$acc, s$acc_se)))
}))

saveRDS(list(xsec = xsec, meta = meta, long = long, long_meta = long_meta),
        file.path(DIR, "cache", "fits.rds"))
saveRDS(scatter_df, file.path(DIR, "cache", "scatter.rds"))
cat(sprintf("\nWrote cache/fits.rds (%d xsec fits) and cache/scatter.rds (%d kids).\n",
            nrow(xsec), nrow(scatter_df)))
