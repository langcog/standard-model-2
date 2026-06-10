## Shared builder for the demographic composite figure (cross-sectional +
## longitudinal). Sourced by both the notebook and the manuscript's
## fig-demographics chunk so the paper and the analysis stay in sync.
##
## Input: `fits` = readRDS("cross_sectional_demographics/cache/fits.rds")
##   $xsec      per-(language,predictor) cross-sectional eff/acc (+ se)
##   $meta      cross-sectional random-effects meta (per predictor x component)
##   $long      longitudinal raw-BLUP effects per by-study unit
##   $long_meta longitudinal meta
## Returns a patchwork object (Panel A by-language over Panel B meta).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(patchwork)
})

make_demographics_composite <- function(fits, mated_min_n = 300,
                                         drop_langs = c("French (French)")) {
  COL_XS <- "#1e88e5"; COL_LO <- "#c41e37"
  pred_lab <- c(sex = "Sex (male vs female)", matEd = "Maternal ed. (per SD)")
  comp_lev <- c("Efficiency", "Acceleration")

  ## consistent language set: matEd languages with n >= threshold, minus drops
  keep_langs <- fits$xsec |>
    filter(predictor == "matEd", n_kids >= mated_min_n,
           !language %in% drop_langs) |>
    pull(language)

  ## longitudinal by-study units -> cross-sectional language names (pool English)
  long_map <- c(thal="English (American)", smith="English (American)",
                marchman="English (American)", norwegian="Norwegian",
                japanese="Japanese")
  long_lang <- fits$long |>
    mutate(language = unname(long_map[language])) |>
    filter(!is.na(language), language %in% keep_langs) |>
    group_by(language, predictor) |>
    summarise(eff = weighted.mean(eff, 1/eff_se^2), eff_se = sqrt(1/sum(1/eff_se^2)),
              acc = weighted.mean(acc, 1/acc_se^2), acc_se = sqrt(1/sum(1/acc_se^2)),
              .groups = "drop")

  to_long <- function(df) {
    df |>
      pivot_longer(c(eff, acc), names_to = "component", values_to = "estimate") |>
      mutate(se = if_else(component == "eff", eff_se, acc_se),
             lo = estimate - 1.96*se, hi = estimate + 1.96*se,
             component = recode(component, eff = "Efficiency", acc = "Acceleration"),
             component = factor(component, levels = comp_lev),
             predictor = factor(pred_lab[predictor], levels = pred_lab)) |>
      select(language, predictor, component, estimate, lo, hi)
  }
  xs <- fits$xsec |> filter(predictor %in% c("sex","matEd"), language %in% keep_langs) |> to_long()
  lo <- to_long(long_lang)

  ## language order: by cross-sectional SEX efficiency (most female-advantaged top)
  ord <- xs |> filter(predictor == pred_lab["sex"], component == "Efficiency") |>
    arrange(estimate) |> pull(language)
  xs$language <- factor(xs$language, levels = ord)
  lo$language <- factor(lo$language, levels = ord)

  meta <- fits$meta |>
    filter(predictor %in% c("sex","matEd")) |>
    mutate(component = recode(component, efficiency = "Efficiency", acceleration = "Acceleration"),
           component = factor(component, levels = comp_lev),
           predictor = factor(pred_lab[predictor], levels = pred_lab))

  ## ---- Panel A: by-language, efficiency + acceleration, both methods ----
  pA <- ggplot(xs, aes(estimate, language)) +
    geom_rect(data = meta, aes(xmin = ci.lb, xmax = ci.ub, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = COL_XS, alpha = 0.08) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = "Cross-sectional"),
                    shape = 16, size = 0.25, fatten = 2.2) +
    geom_point(data = lo, aes(estimate, language, colour = "Longitudinal"),
               shape = 18, size = 2.6) +
    facet_grid(predictor ~ component, scales = "free_x") +
    scale_colour_manual(values = c("Cross-sectional" = COL_XS, "Longitudinal" = COL_LO),
                        name = NULL) +
    labs(x = "effect on latent ability (logits, per SD predictor)", y = NULL,
         title = "A. Per-language demographic effects") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 7.5),
          strip.text = element_text(face = "bold", size = 8.5),
          legend.position = "none",
          panel.grid.minor = element_blank())

  ## ---- Panel B: condensed meta (cross-sectional vs longitudinal) ----
  comb <- bind_rows(meta |> mutate(method = "Cross-sectional"),
                    fits$long_meta |> filter(predictor %in% c("sex","matEd")) |>
                      mutate(component = recode(component, efficiency = "Efficiency",
                                                acceleration = "Acceleration"),
                             component = factor(component, levels = comp_lev),
                             predictor = factor(pred_lab[predictor], levels = pred_lab),
                             method = "Longitudinal"))
  pB <- ggplot(comb, aes(estimate, fct_rev(component), colour = method)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    geom_pointrange(aes(xmin = ci.lb, xmax = ci.ub),
                    position = position_dodge(width = 0.5), size = 0.45, fatten = 2.5) +
    facet_wrap(~ predictor, scales = "free_x") +
    scale_colour_manual(values = c("Cross-sectional" = COL_XS, "Longitudinal" = COL_LO),
                        name = NULL) +
    labs(x = "meta-analytic estimate (logits, per SD predictor)", y = NULL,
         title = "B. Meta-analytic summary") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 11),
          strip.text = element_text(face = "bold", size = 8.5),
          legend.position = "bottom", panel.grid.minor = element_blank())

  (pA / pB) + plot_layout(heights = c(2.5, 1)) &
    theme(plot.margin = margin(4, 8, 4, 4))
}
