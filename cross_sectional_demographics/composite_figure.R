## Shared builder for the demographic composite figure (cross-sectional +
## longitudinal). Sourced by both the notebook and the manuscript's
## fig-demographics chunk so the paper and the analysis stay in sync.
##
## Input: `fits` = readRDS("cross_sectional_demographics/cache/fits.rds")
##   $xsec      per-(language,predictor) cross-sectional eff/acc (+ se)
##   $long      longitudinal raw-BLUP effects per by-study unit
## Sex and maternal ed use INDEPENDENT language sets (each: n >= the
## predictor's threshold), so sex keeps its full breadth even where maternal
## ed is missing. Both predictors report BOTH efficiency and acceleration.
## Returns a patchwork object: sex forest / maternal-ed forest / meta summary.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(forcats)
  library(patchwork); library(metafor)
})

make_demographics_composite <- function(fits, sex_min_n = 300, mated_min_n = 300,
                                         drop_matEd = c("French (French)")) {
  COL_XS <- "#1e88e5"; COL_LO <- "#c41e37"; comp_lev <- c("Efficiency", "Acceleration")

  sex_langs   <- fits$xsec |> filter(predictor == "sex",   n_kids >= sex_min_n) |> pull(language)
  mated_langs <- fits$xsec |> filter(predictor == "matEd", n_kids >= mated_min_n,
                                     !language %in% drop_matEd) |> pull(language)

  ## longitudinal by-study units -> cross-sectional language names (pool English)
  long_map <- c(thal="English (American)", smith="English (American)",
                marchman="English (American)", norwegian="Norwegian",
                japanese="Japanese")
  long_lang <- fits$long |> mutate(language = unname(long_map[language])) |>
    filter(!is.na(language)) |>
    group_by(language, predictor) |>
    summarise(eff = weighted.mean(eff, 1/eff_se^2), eff_se = sqrt(1/sum(1/eff_se^2)),
              acc = weighted.mean(acc, 1/acc_se^2), acc_se = sqrt(1/sum(1/acc_se^2)),
              .groups = "drop")

  to_long <- function(df) df |>
    pivot_longer(c(eff, acc), names_to = "component", values_to = "estimate") |>
    mutate(se = if_else(component == "eff", eff_se, acc_se),
           lo = estimate - 1.96*se, hi = estimate + 1.96*se,
           component = factor(recode(component, eff = "Efficiency", acc = "Acceleration"),
                              levels = comp_lev)) |>
    select(language, predictor, component, estimate, lo, hi)

  ## shared x-limits per component (so sex & maternal-ed panels align and the
  ## efficiency crossover is visible). Acceleration is capped so a couple of
  ## small-n languages with very wide CIs don't blow out the axis.
  all_ll <- bind_rows(
    fits$xsec |> filter(predictor == "sex",   language %in% sex_langs)   |> to_long(),
    fits$xsec |> filter(predictor == "matEd", language %in% mated_langs) |> to_long(),
    long_lang |> filter(predictor == "sex",   language %in% sex_langs)   |> to_long(),
    long_lang |> filter(predictor == "matEd", language %in% mated_langs) |> to_long())
  eff_lim <- all_ll |> filter(component == "Efficiency") |>
    summarise(lo = min(lo, na.rm = TRUE), hi = max(hi, na.rm = TRUE))
  eff_lim <- c(eff_lim$lo, eff_lim$hi) + c(-0.05, 0.05)
  acc_hw  <- max(abs(all_ll$estimate[all_ll$component == "Acceleration"]), na.rm = TRUE) + 0.7
  acc_lim <- c(-acc_hw, acc_hw)
  cap_acc <- function(d) d |>
    mutate(lo = if_else(component == "Acceleration", pmax(lo, acc_lim[1]), lo),
           hi = if_else(component == "Acceleration", pmin(hi, acc_lim[2]), hi))

  ## random-effects meta per (predictor, component) over the DISPLAYED languages
  meta_for <- function(pred, langs) {
    s <- fits$xsec |> filter(predictor == pred, language %in% langs)
    one <- function(b, se) {
      m <- tryCatch(rma(yi = b, sei = se, method = "REML",
                        control = list(stepadj = 0.5, maxiter = 1000)), error = function(e) NULL)
      if (is.null(m)) c(NA, NA, NA) else c(as.numeric(m$beta), m$ci.lb, m$ci.ub)
    }
    e <- one(s$eff, s$eff_se); a <- one(s$acc, s$acc_se)
    tibble(predictor = pred,
           component = factor(comp_lev, levels = comp_lev),
           est = c(e[1], a[1]), ci.lb = c(e[2], a[2]), ci.ub = c(e[3], a[3]))
  }
  meta_sex <- meta_for("sex", sex_langs); meta_med <- meta_for("matEd", mated_langs)

  ## one by-language forest (efficiency | acceleration), x-sec + longitudinal
  forest <- function(pred, langs, meta, title) {
    xs <- fits$xsec |> filter(predictor == pred, language %in% langs) |> to_long()
    lo <- long_lang |> filter(predictor == pred, language %in% langs) |> to_long()
    ord <- xs |> filter(component == "Efficiency") |> arrange(estimate) |> pull(language)
    xs$language <- factor(xs$language, levels = ord); lo$language <- factor(lo$language, levels = ord)
    xs <- cap_acc(xs); lo <- cap_acc(lo)
    lim_df <- tibble(component = factor(rep(comp_lev, each = 2), levels = comp_lev),
                     x = c(eff_lim, acc_lim), language = ord[1])
    ggplot(xs, aes(estimate, language)) +
      geom_rect(data = meta, aes(xmin = ci.lb, xmax = ci.ub, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = COL_XS, alpha = 0.08) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
      geom_blank(data = lim_df, aes(x = x, y = language)) +
      geom_pointrange(aes(xmin = lo, xmax = hi), colour = COL_XS, shape = 16, size = 0.22) +
      geom_point(data = lo, aes(estimate, language), colour = COL_LO, shape = 18, size = 2.4) +
      facet_wrap(~ component, scales = "free_x") +
      labs(x = NULL, y = NULL, title = title) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 10),
            axis.text.y = element_text(size = 7), strip.text = element_text(face = "bold", size = 8.5),
            panel.grid.minor = element_blank())
  }
  pSex <- forest("sex", sex_langs, meta_sex, "A. Sex (female vs male)")
  pMed <- forest("matEd", mated_langs, meta_med, "B. Maternal education (per SD)") +
    labs(x = "effect on latent ability (logits, per SD predictor)")

  ## condensed meta panel (cross-sectional vs longitudinal)
  long_meta <- fits$long_meta |> filter(predictor %in% c("sex","matEd")) |>
    transmute(predictor, component = factor(recode(component, efficiency="Efficiency",
              acceleration="Acceleration"), levels = comp_lev),
              est = estimate, ci.lb, ci.ub, method = "Longitudinal")
  comb <- bind_rows(bind_rows(meta_sex, meta_med) |> mutate(method = "Cross-sectional"), long_meta) |>
    mutate(predictor = factor(recode(predictor, sex = "Sex", matEd = "Maternal ed."),
                              levels = c("Sex", "Maternal ed.")),
           method = factor(method, levels = c("Cross-sectional", "Longitudinal")),
           ci.lb = if_else(component == "Acceleration", pmax(ci.lb, acc_lim[1]), ci.lb),
           ci.ub = if_else(component == "Acceleration", pmin(ci.ub, acc_lim[2]), ci.ub))
  # same grid as A/B (efficiency | acceleration columns; sex / maternal-ed rows)
  # and the same per-component x-limits, so all three panels align.
  limC <- tidyr::crossing(
    predictor = factor(c("Sex","Maternal ed."), levels = c("Sex","Maternal ed.")),
    tibble(component = factor(rep(comp_lev, each = 2), levels = comp_lev),
           x = c(eff_lim, acc_lim))) |>
    mutate(method = factor("Cross-sectional", levels = c("Cross-sectional","Longitudinal")))
  pMeta <- ggplot(comb, aes(est, fct_rev(method), colour = method)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    geom_blank(data = limC, aes(x = x, y = method)) +
    geom_pointrange(aes(xmin = ci.lb, xmax = ci.ub), size = 0.5) +
    facet_grid(predictor ~ component, scales = "free_x") +
    scale_colour_manual(values = c("Cross-sectional" = COL_XS, "Longitudinal" = COL_LO), name = NULL) +
    labs(x = "meta-analytic estimate (logits, per SD)", y = NULL, title = "C. Meta-analytic summary") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10), legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 8.5), panel.grid.minor = element_blank())

  (pSex / pMed / pMeta) +
    plot_layout(heights = c(length(sex_langs) + 3, length(mated_langs) + 3, 15)) &
    theme(plot.margin = margin(3, 8, 3, 4))
}
