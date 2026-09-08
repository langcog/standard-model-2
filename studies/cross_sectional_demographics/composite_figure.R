## Shared builder for the demographic composite figure (cross-sectional +
## longitudinal). Sourced by both the notebook and the paper's figure script
## (paper/figures/fig_demographics.R) so the paper and the analysis stay in sync.
##
## Input: `fits` = readRDS("cross_sectional_demographics/cache/fits.rds")
##   $xsec      per-(language,predictor) cross-sectional eff/acc (+ se)
##   $long      longitudinal raw-BLUP effects per by-study unit
## Sex and maternal ed use INDEPENDENT language sets (each: n >= the
## predictor's threshold), so sex keeps its full breadth even where maternal
## ed is missing. Both predictors report BOTH efficiency and acceleration.
## make_demographics_composite() returns a patchwork: sex forest (tall, left) |
## maternal-ed forest over the meta summary (right); forests_only = TRUE drops the
## meta panel (SI figure). make_demographics_main() is the one-grid main-text
## version (longitudinal datasets + summaries only). The forests share per-component x-limits; the meta
## panel has its own scale, common to both components. Colour = trait (efficiency / acceleration), shape =
## method (filled circle cross-sectional, open diamond longitudinal); the
## paper passes its palette / theme / tinted strips, the notebook uses the
## defaults. Panel letters are added by the caller (patchwork tags).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(forcats)
  library(patchwork); library(metafor)
})

## Cross-sectional meta-analytic main effects shown in Panel C (and the
## per-panel meta bands of the figure). Uses the SAME language filtering as
## make_demographics_composite(), so values match the plot exactly. Returns
## predictor x component with estimate, 95% CI, and p-value. The manuscript's
## fig-demographics inline beta/p values read from this.
panelC_meta <- function(fits, sex_min_n = 300, mated_min_n = 300,
                        drop_matEd = c("French (French)")) {
  sex_langs   <- fits$xsec |> filter(predictor == "sex",   n_kids >= sex_min_n) |> pull(language)
  mated_langs <- fits$xsec |> filter(predictor == "matEd", n_kids >= mated_min_n,
                                     !language %in% drop_matEd) |> pull(language)
  one <- function(b, s) {
    m <- tryCatch(rma(yi = b, sei = s, method = "REML",
                      control = list(stepadj = 0.5, maxiter = 1000)), error = function(e) NULL)
    if (is.null(m)) c(NA, NA, NA, NA) else c(as.numeric(m$beta), m$ci.lb, m$ci.ub, m$pval)
  }
  rows <- function(pred, langs) {
    s <- fits$xsec |> filter(predictor == pred, language %in% langs)
    e <- one(s$eff, s$eff_se); a <- one(s$acc, s$acc_se)
    data.frame(predictor = pred, component = c("efficiency", "acceleration"),
               estimate = c(e[1], a[1]), conf.low = c(e[2], a[2]),
               conf.high = c(e[3], a[3]), p = c(e[4], a[4]), stringsAsFactors = FALSE)
  }
  rbind(rows("sex", sex_langs), rows("matEd", mated_langs))
}

## ---- shared helpers (used by every builder below) ----
LONG_LAB <- c(thal = "English (Thal)", marchman = "English (Marchman)",
              norwegian = "Norwegian", japanese = "Japanese")

## languages displayed per predictor (n >= threshold; maternal-ed drops)
demo_langs <- function(fits, sex_min_n = 300, mated_min_n = 300, drop_matEd = c("French (French)")) list(
  sex   = fits$xsec |> filter(predictor == "sex",   n_kids >= sex_min_n) |> pull(language),
  matEd = fits$xsec |> filter(predictor == "matEd", n_kids >= mated_min_n,
                              !language %in% drop_matEd) |> pull(language))

## random-effects meta per component over a set of languages (cross-sectional)
xs_meta <- function(fits, pred, langs, comp_lev) {
  s <- fits$xsec |> filter(predictor == pred, language %in% langs)
  one <- function(b, se) {
    m <- tryCatch(rma(yi = b, sei = se, method = "REML",
                      control = list(stepadj = 0.5, maxiter = 1000)), error = function(e) NULL)
    if (is.null(m)) c(NA, NA, NA) else c(as.numeric(m$beta), m$ci.lb, m$ci.ub)
  }
  e <- one(s$eff, s$eff_se); a <- one(s$acc, s$acc_se)
  tibble(predictor = pred, component = factor(comp_lev, levels = comp_lev),
         est = c(e[1], a[1]), ci.lb = c(e[2], a[2]), ci.ub = c(e[3], a[3]), k = nrow(s))
}

make_demographics_composite <- function(fits, sex_min_n = 300, mated_min_n = 300,
                                         drop_matEd = c("French (French)"),
                                         pal = c(Efficiency = "#0072B2", Acceleration = "#D55E00"),
                                         base_theme = theme_minimal(base_size = 9),
                                         strips = NULL, forests_only = FALSE) {
  comp_lev <- names(pal)
  SHP <- c("Cross-sectional" = 16, "Longitudinal" = 23)   # filled circle / open diamond
  facet_comp <- function() {                               # efficiency | acceleration columns
    if (is.null(strips)) facet_wrap(~ component, scales = "free_x")
    else ggh4x::facet_wrap2(~ component, scales = "free_x", strip = strips)
  }

  L <- demo_langs(fits, sex_min_n, mated_min_n, drop_matEd)
  sex_langs <- L$sex; mated_langs <- L$matEd

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
  meta_sex <- xs_meta(fits, "sex", sex_langs, comp_lev)
  meta_med <- xs_meta(fits, "matEd", mated_langs, comp_lev)

  ## one by-language forest (efficiency | acceleration), x-sec + longitudinal;
  ## languages ordered by their cross-sectional efficiency effect.
  forest <- function(pred, langs, meta, title, xlab) {
    xs <- fits$xsec |> filter(predictor == pred, language %in% langs) |> to_long() |>
      mutate(method = "Cross-sectional")
    lo <- long_lang |> filter(predictor == pred, language %in% langs) |> to_long() |>
      mutate(method = "Longitudinal")
    ord <- xs |> filter(component == "Efficiency") |> arrange(estimate) |> pull(language)
    lv  <- function(d) mutate(d, language = factor(language, levels = ord),
                              method = factor(method, levels = names(SHP)))
    xs <- cap_acc(lv(xs)); lo <- cap_acc(lv(lo))
    lim_df <- tibble(component = factor(rep(comp_lev, each = 2), levels = comp_lev),
                     x = c(eff_lim, acc_lim), language = factor(ord[1], levels = ord))
    ggplot(xs, aes(estimate, language)) +
      geom_rect(data = meta, aes(xmin = ci.lb, xmax = ci.ub, ymin = -Inf, ymax = Inf, fill = component),
                inherit.aes = FALSE, alpha = 0.15) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
      geom_errorbar(aes(xmin = lo, xmax = hi, colour = component), width = 0, linewidth = 0.35) +
      geom_point(aes(colour = component, shape = method), size = 1.1) +
      geom_point(data = lo, aes(estimate, language, colour = component, shape = method),
                 fill = "white", size = 1.9, stroke = 0.7) +
      # x-limits layer LAST: a single-value layer trained first makes ggplot
      # re-sort the discrete axis (this is what made the languages alphabetical)
      geom_blank(data = lim_df, aes(x = x, y = language)) +
      scale_colour_manual(values = pal, guide = "none") +
      scale_fill_manual(values = pal, guide = "none") +
      scale_shape_manual(values = SHP, name = NULL, drop = FALSE) +
      guides(shape = guide_legend(override.aes = list(colour = "grey20", fill = "white",
                                                      size = c(1.4, 1.9)))) +
      facet_comp() +
      labs(x = xlab, y = NULL, title = title) +
      base_theme + theme(panel.grid.minor = element_blank())
  }
  pSex <- forest("sex", sex_langs, meta_sex, "Sex", "effect of female vs male (logits)") +
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
  pMed <- forest("matEd", mated_langs, meta_med, "Maternal education",
                 "effect per SD of maternal education (logits)") +
    theme(legend.position = "none")

  ## condensed meta panel (cross-sectional vs longitudinal)
  long_meta <- fits$long_meta |> filter(predictor %in% c("sex","matEd")) |>
    transmute(predictor, component = factor(recode(component, efficiency="Efficiency",
              acceleration="Acceleration"), levels = comp_lev),
              est = estimate, ci.lb, ci.ub, method = "Longitudinal")
  comb <- bind_rows(bind_rows(meta_sex, meta_med) |> mutate(method = "Cross-sectional"), long_meta) |>
    mutate(predictor = factor(recode(predictor, sex = "Sex", matEd = "Maternal ed."),
                              levels = c("Sex", "Maternal ed.")),
           method = factor(method, levels = names(SHP)),
           row = factor(method, levels = rev(names(SHP))))   # cross-sectional on top
  # same grid as the forests (efficiency | acceleration columns; sex / maternal-ed
  # rows) but its OWN x-scale, shared by the two components: the forests'
  # acceleration axis is set by the widest per-language CIs (+/- 2.4) and would
  # squash the meta estimates at zero; a common scale also makes efficiency and
  # acceleration directly comparable here.
  facet_meta <- if (is.null(strips)) facet_grid(predictor ~ component) else
    ggh4x::facet_grid2(predictor ~ component, strip = strips)
  pMeta <- ggplot(comb, aes(est, row, colour = component, shape = method)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), width = 0, linewidth = 0.45) +
    geom_point(fill = "white", size = 2, stroke = 0.7) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_shape_manual(values = SHP, guide = "none") +
    facet_meta +
    labs(x = "meta-analytic estimate (logits)", y = NULL, title = "Meta-analytic summary") +
    base_theme + theme(panel.grid.minor = element_blank())

  if (forests_only)   # SI version: the two by-language forests, no meta panel
    return((pSex | ((pMed / plot_spacer()) +
                      plot_layout(heights = c(length(mated_langs) + 3,
                                              length(sex_langs) - length(mated_langs))))) +
             plot_layout(widths = c(1, 1)))
  (pSex | ((pMed / pMeta) + plot_layout(heights = c(length(mated_langs) + 4, 10)))) +
    plot_layout(widths = c(1, 1))
}

## Main-text version (2026-09-04 prototype): ONE grid -- columns = trait,
## row facets = predictor -- showing the longitudinal datasets by name (raw-BLUP
## regressions; filled circles with 95% CI lines) above two summary rows drawn
## in forest-plot convention as open diamonds whose width IS the 95% CI: the
## longitudinal meta-analysis over datasets and the cross-sectional
## meta-analysis over the displayed Wordbank languages. No bands, no legend
## (the row labels name every marker). The cross-sectional by-language forests
## are the SI figure (make_demographics_composite(..., forests_only = TRUE)).
## Shared x-scale; sized for a single column.
make_demographics_main <- function(fits, sex_min_n = 300, mated_min_n = 300,
                                   drop_matEd = c("French (French)"),
                                   pal = c(Efficiency = "#0072B2", Acceleration = "#D55E00"),
                                   base_theme = theme_minimal(base_size = 9),
                                   strips = NULL) {
  comp_lev <- names(pal)
  L <- demo_langs(fits, sex_min_n, mated_min_n, drop_matEd)
  PRED <- c(sex = "Sex\n(female \u2212 male)", matEd = "Maternal ed.\n(per SD)")
  KIND <- c("Longitudinal dataset", "Longitudinal summary", "Cross-sectional summary")

  ## longitudinal datasets (95% CI from the BLUP regressions)
  dat <- fits$long |> filter(predictor %in% names(PRED), language %in% names(LONG_LAB)) |>
    transmute(predictor, row = unname(LONG_LAB[language]), kind = KIND[1],
              Efficiency_est = eff, Efficiency_lo = eff - 1.96 * eff_se, Efficiency_hi = eff + 1.96 * eff_se,
              Acceleration_est = acc, Acceleration_lo = acc - 1.96 * acc_se, Acceleration_hi = acc + 1.96 * acc_se) |>
    pivot_longer(-c(predictor, row, kind), names_to = c("component", ".value"), names_sep = "_")
  ## longitudinal meta (k datasets) and cross-sectional meta (k languages)
  lmeta <- fits$long_meta |> filter(predictor %in% names(PRED)) |>
    transmute(predictor, row = sprintf("Longitudinal (%d datasets)", k), kind = KIND[2],
              component = recode(component, efficiency = "Efficiency", acceleration = "Acceleration"),
              est = estimate, lo = ci.lb, hi = ci.ub)
  xmeta <- bind_rows(xs_meta(fits, "sex", L$sex, comp_lev), xs_meta(fits, "matEd", L$matEd, comp_lev)) |>
    transmute(predictor, row = sprintf("Cross-sectional (%d languages)", k), kind = KIND[3],
              component = as.character(component), est, lo = ci.lb, hi = ci.ub)
  ## row order, bottom -> top within each facet: the two summaries, then datasets
  lev <- c(unique(xmeta$row), unique(lmeta$row), rev(unname(LONG_LAB)))
  d <- bind_rows(dat, lmeta, xmeta) |>
    mutate(component = factor(component, levels = comp_lev),
           predictor = factor(unname(PRED[predictor]), levels = unname(PRED)),
           kind = factor(kind, levels = KIND), row = factor(row, levels = lev))
  ## summary diamonds: width = 95% CI, at numeric rows 1 (cross-sectional) and 2
  ## (longitudinal) of each facet -- the two bottom rows by construction of `lev`
  H <- 0.36
  dia <- d |> filter(kind != KIND[1]) |>
    mutate(y = if_else(kind == KIND[3], 1, 2), id = paste(predictor, component, row)) |>
    rowwise() |>
    reframe(predictor, component, id, x = c(lo, est, hi, est), y = c(y, y + H, y, y - H))
  facet <- if (is.null(strips)) facet_grid(predictor ~ component, scales = "free_y", space = "free_y") else
    ggh4x::facet_grid2(predictor ~ component, scales = "free_y", space = "free_y", strip = strips)
  ggplot(d |> filter(kind == KIND[1]), aes(est, row, colour = component)) +
    geom_blank(data = d, aes(est, row)) +   # FIRST: trains the axis on all rows in `lev` order
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = 2.5, colour = "grey75", linewidth = 0.3) +   # summaries | datasets
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4) +
    geom_point(shape = 16, size = 1.6) +
    geom_polygon(data = dia, aes(x, y, group = id, colour = component), inherit.aes = FALSE,
                 fill = "white", linewidth = 0.5) +
    scale_colour_manual(values = pal, guide = "none") +
    facet +
    labs(x = "effect, 95% CI\n(logits; acceleration in \u03ba units)", y = NULL) +
    base_theme + theme(panel.grid.minor = element_blank(),
                       strip.text.y = element_text(lineheight = 0.9),
                       axis.title.x = element_text(lineheight = 0.95))
}
