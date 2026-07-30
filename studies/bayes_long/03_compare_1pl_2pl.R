## 03_compare_1pl_2pl.R -- is the acceleration result robust to relaxing equal
## discrimination? Compares M3 (1PL/Rasch) with M3-2PL on the same _a3 data.
##
## WHY RAW kappa CANNOT BE COMPARED DIRECTLY. Standard IRT makes 1PL and 2PL abilities
## comparable by pinning the population scale (theta ~ N(0,1)). We cannot do that: theta
## here is structured, theta_i(t) = xi_i + kappa_i log(t/a0) + log H, and its population
## SD, sqrt(sigma_a^2 + sigma_b^2 L^2 + 2 rho sigma_a sigma_b L), GROWS WITH AGE as the
## fan opens -- there is no single sigma_theta to fix. In the 1PL the scale is pinned by
## the link instead (lambda == 1, so one theta unit = one logit); in the 2PL one theta
## unit = lambda_j logits, item-specific, and the scale is set only by our normalization
## of lambda (geometric mean pinned to 1). Different normalizations give different kappa.
##
## THREE SCALE-FREE COMPARISONS, in increasing order of how much they settle:
##
## 1. delta log kappa across datasets. A pure change of units is a single multiplicative
##    factor c, i.e. an ADDITIVE and CONSTANT shift -log(c) in log kappa. So if the five
##    datasets shift log kappa by the same amount, the difference is convention; if they
##    shift by different amounts, the 2PL is revising the growth story. This is the test
##    that distinguishes the two readings.
##
## 2. Age of acquisition, t_50 -- EXACTLY invariant to the lambda scale. Solving
##    theta(t) = delta_j gives log(t_50/a0) = (delta_j - log_H - mu_xi)/kappa. Under the
##    invariance (lambda -> c lambda, theta - delta -> (theta - delta)/c, absorbed by the
##    free intercept since log_H is fixed data) every term rescales so that t_50 is
##    unchanged. Note lambda_j does not enter: the threshold where p = 0.5 is theta =
##    delta_j whatever the discrimination. So matching AoAs mean the models agree on the
##    substance and the kappa gap is pure units. This is also the quantity Fig 2 uses.
##
## 3. LOO. Does the 2PL actually predict better? (Requires the LOO branch for the 2PL
##    linear predictor, added to 01_fit.R.)
##
## Usage:  Rscript studies/bayes_long/03_compare_1pl_2pl.R
## Output: paper/cache/si_2pl.rds + console report

suppressPackageStartupMessages({library(dplyr)})
SUMM  <- file.path("fits", "bayes_long", "summaries")
CACHE <- file.path("paper", "cache")
DATASETS <- c(thal = "English (Thal)", smith = "English (Smith)",
              marchman = "English (Marchman)", norwegian = "Norwegian",
              japanese = "Japanese")
AGES <- seq(10, 36, by = 1)

one <- function(slug, label) {
  f1 <- file.path(SUMM, sprintf("%s_a3_m3.summary.rds", slug))
  f2 <- file.path(SUMM, sprintf("%s_a3_m32pl.summary.rds", slug))
  p1 <- file.path(SUMM, sprintf("%s_a3_m3_psi.csv", slug))
  p2 <- file.path(SUMM, sprintf("%s_a3_m32pl_psi.csv", slug))
  if (!all(file.exists(c(f1, f2, p1, p2)))) { cat("  pending:", slug, "\n"); return(NULL) }
  s1 <- as.data.frame(readRDS(f1)); s2 <- as.data.frame(readRDS(f2))
  g <- function(s, v, q = "median") { x <- s[[q]][s$variable == v]; if (!length(x)) NA_real_ else x }
  b  <- readRDS(file.path("fits", "bayes_long", sprintf("bundle_%s_a3.rds", slug)))
  a0 <- b$stan_data$a0; lH <- b$stan_data$log_H
  d1 <- read.csv(p1); d2 <- read.csv(p2)

  k1 <- g(s1, "kappa_pop"); k2 <- g(s2, "kappa_pop")
  m1 <- g(s1, "mu_xi");     m2 <- g(s2, "mu_xi")

  ## ---- (2) AoA, the lambda-scale-invariant observable ----
  ## t_50 = a0 * exp((delta_j - log_H - mu_xi)/kappa); lambda_j does not enter.
  aoa <- inner_join(transmute(d1, item, d_1pl = delta_j),
                    transmute(d2, item, d_2pl = delta_j, lambda), by = "item") |>
    mutate(t50_1pl = a0 * exp((d_1pl - lH - m1) / k1),
           t50_2pl = a0 * exp((d_2pl - lH - m2) / k2)) |>
    filter(is.finite(t50_1pl), is.finite(t50_2pl), t50_1pl > 4, t50_1pl < 120,
           t50_2pl > 4, t50_2pl < 120)

  ## ---- model-implied population trajectory (proportion of items produced) ----
  ## Median child: xi = mu_xi, kappa = kappa_pop. 1PL vs 2PL over the SAME item set.
  aoa_items <- aoa$item
  dd1 <- d1$delta_j[match(aoa_items, d1$item)]
  dd2 <- d2$delta_j[match(aoa_items, d2$item)]
  lam <- d2$lambda[match(aoa_items, d2$item)]
  traj <- data.frame(age = AGES) |> rowwise() |>
    mutate(L    = log(age / a0),
           p1pl = mean(plogis(m1 + k1 * L + lH - dd1)),
           p2pl = mean(plogis(lam * (m2 + k2 * L + lH - dd2)))) |>
    ungroup() |> mutate(slug = slug, lang = label)

  ## ---- (3) LOO ----
  l1 <- file.path(SUMM, sprintf("%s_a3_m3.loo.rds", slug))
  l2 <- file.path(SUMM, sprintf("%s_a3_m32pl.loo.rds", slug))
  elpd_diff <- se_diff <- NA_real_; winner <- NA_character_
  if (file.exists(l1) && file.exists(l2)) {
    cmp <- loo::loo_compare(list(`1PL` = readRDS(l1), `2PL` = readRDS(l2)))
    winner <- rownames(cmp)[1]
    elpd_diff <- cmp[2, "elpd_diff"]; se_diff <- cmp[2, "se_diff"]
  }

  list(
    summary = data.frame(
      slug = slug, lang = label,
      kappa_1pl = k1, kappa_2pl = k2,
      log_kappa_1pl = log(k1), log_kappa_2pl = log(k2),
      d_log_kappa = log(k2) - log(k1),          # constant across datasets => units, not substance
      sigma_a_1pl = g(s1,"sigma_a"), sigma_a_2pl = g(s2,"sigma_a"),
      sigma_b_1pl = g(s1,"sigma_b"), sigma_b_2pl = g(s2,"sigma_b"),
      tau_delta_1pl = g(s1,"tau_delta"), tau_delta_2pl = g(s2,"tau_delta"),
      sigma_lambda = g(s2,"sigma_lambda"),
      lambda_p10 = g(s2,"lambda_p10"), lambda_p90 = g(s2,"lambda_p90"),
      cor_lambda_delta = g(s2,"cor_lambda_delta"),
      cor_lambda_delta_q5 = g(s2,"cor_lambda_delta","q5"),
      cor_lambda_delta_q95 = g(s2,"cor_lambda_delta","q95"),
      ## AoA agreement: the substantive check
      n_items = nrow(aoa),
      aoa_cor = cor(aoa$t50_1pl, aoa$t50_2pl),
      aoa_mad = median(abs(aoa$t50_2pl - aoa$t50_1pl)),
      aoa_med_1pl = median(aoa$t50_1pl), aoa_med_2pl = median(aoa$t50_2pl),
      ## trajectory agreement (max abs difference in proportion produced)
      traj_maxdiff = max(abs(traj$p1pl - traj$p2pl)),
      rhat_2pl = max(s2$rhat, na.rm = TRUE),
      loo_winner = winner, loo_elpd_diff = elpd_diff, loo_se_diff = se_diff,
      row.names = NULL),
    aoa = mutate(aoa, slug = slug, lang = label), traj = traj)
}

cat("1PL vs 2PL comparison (M3 on the _a3 data)\n")
R <- Filter(Negate(is.null), Map(one, names(DATASETS), DATASETS))
if (!length(R)) { cat("no dataset has both fits yet\n"); quit(save = "no") }
S <- bind_rows(lapply(R, `[[`, "summary")) |>
  mutate(lang = factor(lang, levels = unname(DATASETS))) |> arrange(lang)
out <- list(summary = S,
            aoa  = bind_rows(lapply(R, `[[`, "aoa")),
            traj = bind_rows(lapply(R, `[[`, "traj")))
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
saveRDS(out, file.path(CACHE, "si_2pl.rds"))

cat("\n=== discrimination actually varies? ===\n")
print(S |> transmute(Dataset = lang,
  sigma_lambda = sprintf("%.3f", sigma_lambda),
  `lambda p10-p90` = sprintf("%.2f-%.2f", lambda_p10, lambda_p90),
  `cor(log lambda, delta)` = sprintf("%+.3f [%+.3f, %+.3f]",
      cor_lambda_delta, cor_lambda_delta_q5, cor_lambda_delta_q95)) |> as.data.frame(),
  row.names = FALSE)

cat("\n=== (1) is the kappa shift a constant change of units? ===\n")
print(S |> transmute(Dataset = lang,
  kappa_1pl = sprintf("%.2f", kappa_1pl), kappa_2pl = sprintf("%.2f", kappa_2pl),
  `d log kappa` = sprintf("%+.3f", d_log_kappa),
  `implied factor` = sprintf("%.2fx", exp(d_log_kappa))) |> as.data.frame(), row.names = FALSE)
if (nrow(S) > 1) {
  cat(sprintf("\n  d log kappa: mean %+.3f, sd %.3f, range %+.3f to %+.3f\n",
              mean(S$d_log_kappa), sd(S$d_log_kappa), min(S$d_log_kappa), max(S$d_log_kappa)))
  cat("  A SMALL sd means one common rescaling -> kappa's magnitude is convention-dependent,\n")
  cat("  and the acceleration claim is unaffected. A LARGE sd means the 2PL revises the\n")
  cat("  growth story per dataset, which WOULD matter substantively.\n")
}

cat("\n=== (2) AoA (lambda-invariant) and trajectory agreement ===\n")
print(S |> transmute(Dataset = lang, n_items,
  `cor(t50)` = sprintf("%.4f", aoa_cor),
  `median |dt50| (mo)` = sprintf("%.2f", aoa_mad),
  `median t50 1PL/2PL` = sprintf("%.1f / %.1f", aoa_med_1pl, aoa_med_2pl),
  `max |dP(t)|` = sprintf("%.4f", traj_maxdiff)) |> as.data.frame(), row.names = FALSE)

cat("\n=== (3) LOO: does the 2PL predict better? ===\n")
print(S |> transmute(Dataset = lang, winner = loo_winner,
  `dELPD (loser)` = ifelse(is.na(loo_elpd_diff), "pending",
                           sprintf("%.0f (%.0f)", loo_elpd_diff, loo_se_diff)),
  `max rhat 2PL` = sprintf("%.3f", rhat_2pl)) |> as.data.frame(), row.names = FALSE)
cat("\nwrote", file.path(CACHE, "si_2pl.rds"), "\n")
