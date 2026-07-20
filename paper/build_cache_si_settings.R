## build_cache_si_settings.R -- the SI "does it hold across settings?" table.
##
## The main text reports ONE setting: 3+ administrations, separate per-dataset M3 fits.
## That is the conservative choice -- a two-administration child yields a single slope
## with zero residual df, so their kappa_i is confounded with the noise in xi_i, and the
## resulting sigma_b is inflated by measurement error rather than by real heterogeneity.
## This table is the transparency claim behind that choice: the same qualitative result
## (acceleration exists and varies) appears in every setting we can fit.
##
## Rows = dataset x setting:
##   "2+ separate"  all longitudinal children, per-dataset M3   (SFX "")
##   "3+ separate"  >=3 administrations, per-dataset M3         (SFX "_a3")  <- main text
##   "mega (3+)"    one cross-dataset model, partially pooled scales (pool_a3)
##
## Columns: n_kids/n_admins/n_obs, kappa, sigma_a, sigma_b, rho, and the M2->M3 LOO gap.
##
## NB on the gap: raw elpd differences scale with the NUMBER OF OBSERVATIONS, so a 2+
## sample with 3x the children shows a bigger gap almost mechanically. Comparing raw
## gaps across thresholds is therefore meaningless. We report gap_per_obs (elpd
## difference per binary response) alongside the raw value; only that is comparable.
##
## Run:  Rscript paper/build_cache_si_settings.R
suppressPackageStartupMessages({library(dplyr); library(here); library(loo)})

SUMM <- here("fits", "bayes_long", "summaries")
BL   <- here("fits", "bayes_long")
CACHE <- here("paper", "cache")

DATASETS <- c(thal = "English (Thal)", smith = "English (Smith)",
              marchman = "English (Marchman)", norwegian = "Norwegian",
              japanese = "Japanese")
SETTINGS <- c("2+ separate" = "", "3+ separate" = "_a3")

## ---- per-dataset separate fits, at each threshold -----------------------
one <- function(slug, label, sfx, setting) {
  sf <- file.path(SUMM, paste0(slug, sfx, "_m3.summary.rds"))
  bf <- file.path(BL,   paste0("bundle_", slug, sfx, ".rds"))
  if (!file.exists(sf) || !file.exists(bf)) {
    cat("  skip ", slug, " (", setting, ": missing fit/bundle)\n", sep = ""); return(NULL) }
  s <- as.data.frame(readRDS(sf)); g <- function(v) s$median[s$variable == v]
  sd0 <- readRDS(bf)$stan_data

  ## M2 -> M3 gap: pairwise loo_compare, so se_diff is the SE of THIS difference
  ## (the 4-way ladder's se_diff is relative to the best model, a different quantity).
  f2 <- file.path(SUMM, paste0(slug, sfx, "_m2.loo.rds"))
  f3 <- file.path(SUMM, paste0(slug, sfx, "_m3.loo.rds"))
  gap <- se <- n_obs <- NA_real_
  if (file.exists(f2) && file.exists(f3)) {
    l2 <- readRDS(f2); l3 <- readRDS(f3)
    cmp <- loo::loo_compare(list(M2 = l2, M3 = l3))
    gap <- l3$estimates["elpd_loo", "Estimate"] - l2$estimates["elpd_loo", "Estimate"]
    se  <- cmp[2, "se_diff"]                       # SE of the pairwise difference
    n_obs <- nrow(l3$pointwise)
  }
  data.frame(dataset = label, setting = setting,
             n_kids = sd0$I, n_admins = sd0$A, n_obs_fit = sd0$N,
             kappa = 1 + g("delta"), sigma_a = g("sigma_a"),
             sigma_b = g("sigma_b"), rho = g("rho_ab"),
             gap_m2_m3 = gap, gap_se = se, n_obs_loo = n_obs,
             row.names = NULL)
}

cat("SI settings table: per-dataset separate fits\n")
sep <- bind_rows(unname(unlist(recursive = FALSE, lapply(names(SETTINGS), function(nm)
  Map(function(s, l) one(s, l, SETTINGS[[nm]], nm), names(DATASETS), DATASETS)))))

## ---- mega-model row (one fit, all datasets, partially pooled scales) ----
## Skipped silently until the fit lands; rerun this script to fill it in.
mega <- NULL
pf <- file.path(SUMM, "pool_a3.summary.rds")
pb <- file.path(BL, "bundle_pool_a3.rds")
if (file.exists(pf) && file.exists(pb)) {
  s <- as.data.frame(readRDS(pf)); sd0 <- readRDS(pb)$stan_data
  gi <- function(v, i) s$median[s$variable == sprintf("%s[%d]", v, i)]
  ## dataset order in the pooled bundle -- carried by the bundle, not assumed
  dnames <- if (!is.null(sd0$ds_names)) sd0$ds_names else names(DATASETS)
  mega <- bind_rows(lapply(seq_along(dnames), function(d) {
    slug <- dnames[d]
    data.frame(dataset = unname(DATASETS[slug]), setting = "mega (3+)",
               n_kids = sum(sd0$child_ds == d),
               n_admins = sum(sd0$child_ds[sd0$admin_to_child] == d),
               n_obs_fit = NA_real_,
               kappa = gi("kappa_pop", d), sigma_a = gi("sigma_a", d),
               sigma_b = gi("sigma_b", d), rho = gi("rho_ab", d),
               ## the mega has no M2 rung, so no ladder gap
               gap_m2_m3 = NA_real_, gap_se = NA_real_, n_obs_loo = NA_real_,
               row.names = NULL) }))
  cat(sprintf("  mega: pooled sigma_b = %.2f, pooled sigma_a = %.2f\n",
              s$median[s$variable == "sigma_b_pop"], s$median[s$variable == "sigma_a_pop"]))
} else {
  cat("  mega: pool_a3 fit not present yet -- rerun to add the mega rows\n")
}

si_settings <- bind_rows(sep, mega) |>
  mutate(gap_per_obs = gap_m2_m3 / n_obs_loo,          # the only cross-setting-comparable gap
         dataset = factor(dataset, levels = unname(DATASETS)),
         setting = factor(setting, levels = c("2+ separate", "3+ separate", "mega (3+)"))) |>
  arrange(dataset, setting)

saveRDS(si_settings, file.path(CACHE, "si_settings.rds"))
cat(sprintf("\nWrote %s (%d rows, %d settings)\n",
            file.path(CACHE, "si_settings.rds"), nrow(si_settings),
            dplyr::n_distinct(si_settings$setting)))
print(si_settings |>
      transmute(dataset, setting, n_kids, kappa = round(kappa, 1),
                sigma_b = round(sigma_b, 2), gap = round(gap_m2_m3),
                se = round(gap_se), per_obs = signif(gap_per_obs, 3)) |>
      as.data.frame())
