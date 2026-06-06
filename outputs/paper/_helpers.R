## _helpers.R — utility functions for chunks in standard_model.qmd
##
## Pattern: each `make_*` function returns a ggplot (or patchwork) object.
## Chunks call these; the actual data prep happens here or in the
## cache builders under outputs/paper/cache/.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(latex2exp); library(here); library(knitr); library(scales)
})

# ---- consistent palettes -----------------------------------------------
LANG_PAL <- c("english_american" = "#E69F00",
              "norwegian"        = "#56B4E9",
              "french_quebecois" = "#009E73",
              "japanese"         = "#D55E00",
              "finnish"          = "#CC79A7",
              "english_british"  = "#999999")

LANG_LABELS <- c("english_american" = "English (American)",
                 "norwegian"        = "Norwegian",
                 "french_quebecois" = "French (Quebecois)",
                 "japanese"         = "Japanese",
                 "finnish"          = "Finnish",
                 "english_british"  = "English (British)")

STUDY_PAL <- c("BabyView"  = "#E69F00", "SEEDLingS" = "#56B4E9",
               "AM2018"    = "#009E73", "FMW2013"   = "#D55E00")

# Primary four languages for the main-text figures
PAPER_LANGS <- c("english_american", "norwegian", "french_quebecois", "japanese")

# ---- repo-root helper --------------------------------------------------
# All paths inside this file are resolved relative to repo root so the
# qmd can be rendered from anywhere.
repo_path <- function(...) here::here(...)

# ---- caching --------------------------------------------------------
# Tiny pattern: load_or_build(path, builder) caches expensive computation.
load_or_build <- function(path, builder) {
  if (file.exists(path)) return(readRDS(path))
  obj <- builder()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path)
  obj
}

# ---- Figure 1: schematic --------------------------------------------
# Five-panel schematic walking through the model ladder. Each panel
# shows three illustrative quantile trajectories (25 / 50 / 75) plus
# its governing equation as a latex annotation.
make_fig1_schematic <- function() {
  age <- seq(8, 30, length.out = 80)
  a0  <- 19
  L   <- log(age / a0)

  # parameters chosen to be illustrative, not estimated
  variants <- list(
    list(name = "(A) pure accumulator", xi = 0,  kappa = 1,
         eq = r"($\theta = \log(r) + \log(t)$)"),
    list(name = "(B) + efficiency variance", xi = c(-1.5, 0, 1.5), kappa = 1,
         eq = r"($\theta = \log(\alpha_i) + \log(r) + \log(t)$)"),
    list(name = "(C) + acceleration", xi = c(-1.5, 0, 1.5), kappa = 2.5,
         eq = r"($\theta = \log(\alpha_i) + \log(r) + \kappa \cdot \log(t)$)"),
    list(name = "(D) + per-kid acceleration", xi = c(-1.5, 0, 1.5),
         kappa = c(1.5, 2.5, 3.5),
         eq = r"($\theta = \log(\alpha_i) + \log(r) + \kappa_i \cdot \log(t)$)"),
    list(name = "(E) + per-kid input variance", xi = c(-1.5, 0, 1.5),
         kappa = c(1.5, 2.5, 3.5), input = c(-0.6, 0, 0.6),
         eq = r"($\theta = \log(\alpha_i) + \log(r_i) + \kappa_i \cdot \log(t)$)")
  )

  build_panel <- function(v) {
    xi  <- if (length(v$xi) == 1) rep(v$xi, 3) else v$xi
    kap <- if (length(v$kappa) == 1) rep(v$kappa, 3) else v$kappa
    inp <- if (is.null(v$input)) rep(0, 3) else v$input
    qd  <- bind_rows(lapply(seq_along(xi), function(i) {
      data.frame(quantile = c("25%", "50%", "75%")[i],
                 age = age,
                 theta = xi[i] + inp[i] + kap[i] * L)
    }))
    qd$quantile <- factor(qd$quantile, levels = c("25%", "50%", "75%"))
    ggplot(qd, aes(age, theta, color = quantile, linetype = quantile)) +
      geom_line(linewidth = 0.7) +
      annotate("label", x = 8.5, y = 7.5, hjust = 0, vjust = 1, size = 2.6,
               label.size = 0, fill = scales::alpha("white", 0.7),
               label = latex2exp::TeX(v$eq, output = "character"),
               parse = TRUE) +
      scale_color_manual(values = c("25%" = "#888888", "50%" = "black",
                                    "75%" = "#888888"),
                         name = "child quantile") +
      scale_linetype_manual(values = c("25%" = "dashed", "50%" = "solid",
                                       "75%" = "dotted"),
                            name = "child quantile") +
      coord_cartesian(ylim = c(-3, 8)) +
      labs(x = "age (months)", y = expression(theta[i](t)),
           title = v$name) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold"),
            legend.position = "bottom",
            plot.margin = margin(2, 4, 2, 4))
  }

  panels <- lapply(variants, build_panel)
  # 5-panel row; collect legends to one (will appear at bottom)
  wrap_plots(panels, ncol = 5) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# ---- Figure 2: glmer ladder predictions per language ----------------
make_fig2_glmer_ladder <- function() {
  pred_path <- repo_path("outputs", "paper", "cache", "fig2_glmer_ladder.rds")
  d <- load_or_build(pred_path, function() {
    raw <- read.csv(repo_path("outputs", "glmer_ladder",
                              "predictions_quantiles.csv"))
    raw |>
      filter(lang_slug %in% PAPER_LANGS,
             model %in% c("A", "B_log", "C_log", "D_log")) |>
      mutate(model_label = recode(model,
               "A"     = "(A) accumulator",
               "B_log" = "(B) + pop. acceleration",
               "C_log" = "(C) + per-kid efficiency",
               "D_log" = "(D) + per-kid acceleration"))
  })
  d$lang_slug <- factor(d$lang_slug, levels = PAPER_LANGS,
                        labels = LANG_LABELS[PAPER_LANGS])
  ggplot(d, aes(age, vocab, group = interaction(model, qprob),
                color = factor(qprob))) +
    geom_line(linewidth = 0.55) +
    facet_grid(lang_slug ~ model_label, scales = "free_y") +
    scale_color_manual(values = c("0.1" = "#cbd5e8", "0.25" = "#a6bddb",
                                  "0.5" = "#2b8cbe", "0.75" = "#a6bddb",
                                  "0.9" = "#cbd5e8"),
                       name = "quantile") +
    labs(x = "age (months)", y = "vocabulary (items produced)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8))
}

# ---- Table 1: dataset characteristics --------------------------------
# Pulled from bundles + glmer survey. Citation column intentionally
# placeholder — the author fills citation markers.
make_table1_datasets <- function() {
  cache_path <- repo_path("outputs", "paper", "cache", "table1_datasets.rds")
  load_or_build(cache_path, function() {
    rows <- list(
      # Wordbank longitudinal
      data.frame(panel = "Wordbank longitudinal", dataset = "English (American)",
                 citation = "[@CITE]", n_kids = NA_integer_,
                 admins_per_kid = NA_character_, age_range = NA_character_,
                 input_info = "—"),
      data.frame(panel = "Wordbank longitudinal", dataset = "Norwegian",
                 citation = "[@CITE]", n_kids = NA_integer_,
                 admins_per_kid = NA_character_, age_range = NA_character_,
                 input_info = "—"),
      data.frame(panel = "Wordbank longitudinal", dataset = "French (Quebecois)",
                 citation = "[@CITE]", n_kids = NA_integer_,
                 admins_per_kid = NA_character_, age_range = NA_character_,
                 input_info = "—"),
      data.frame(panel = "Wordbank longitudinal", dataset = "Japanese",
                 citation = "[@CITE]", n_kids = NA_integer_,
                 admins_per_kid = NA_character_, age_range = NA_character_,
                 input_info = "—"),
      # IO datasets
      data.frame(panel = "Input-observed (IO)", dataset = "BabyView",
                 citation = "[@CITE]", n_kids = 22L,
                 admins_per_kid = "5 (1–8)", age_range = "6–32 mo",
                 input_info = "Head-cam, n=5739 videos"),
      data.frame(panel = "Input-observed (IO)", dataset = "SEEDLingS",
                 citation = "[@CITE]", n_kids = 44L,
                 admins_per_kid = "12 (6–13)", age_range = "6–17 mo",
                 input_info = "LENA, n=525 day recordings"),
      data.frame(panel = "Input-observed (IO)", dataset = "AM2018",
                 citation = "[@CITE]", n_kids = 66L,
                 admins_per_kid = "3 (1–3)", age_range = "16–24 mo",
                 input_info = "LENA, n=126 recordings"),
      data.frame(panel = "Input-observed (IO)", dataset = "FMW2013",
                 citation = "[@CITE]", n_kids = 51L,
                 admins_per_kid = "3 (1–3)", age_range = "18–30 mo",
                 input_info = "LENA, n=51 recordings")
    )
    # try to fill Wordbank panel n_kids from the glmer ladder survey csv
    survey_n <- tryCatch({
      s <- read.csv(repo_path("outputs", "glmer_ladder",
                              "00_language_survey.csv"))
      # accommodate either schema
      lang_col <- intersect(c("language", "lang_label", "display", "name"),
                            names(s))[1]
      n_col    <- intersect(c("n_kids", "kids", "n"), names(s))[1]
      if (is.na(lang_col) || is.na(n_col)) NULL else
        setNames(s[[n_col]], s[[lang_col]])
    }, error = function(e) NULL)
    out <- bind_rows(rows)
    if (!is.null(survey_n))
      out$n_kids <- ifelse(is.na(out$n_kids),
                            unname(survey_n[out$dataset]),
                            out$n_kids)
    out
  })
}

# ---- Figure 3: demographic predictors of slope / intercept ---------
# Initially uses the previously generated panel PNG; chunked rebuild
# from cached per-kid estimates is the planned upgrade.
fig3_png_path <- function() {
  repo_path("outputs", "figs", "longitudinal", "predictors_alpha_zeta_panel.png")
}

# ---- Figure 5: IO ladder — growth curves + variance partition -----
make_fig5_io_panel <- function() {
  cache_path <- repo_path("outputs", "paper", "cache", "fig5_io_summary.rds")
  s <- load_or_build(cache_path, function() {
    b <- readRDS(repo_path("fits", "io_pooled_subset_data.rds"))
    fit <- readRDS(repo_path("fits", "io_pooled_widedelta.rds"))
    fitg <- readRDS(repo_path("fits", "io_pooled_gamma_widedelta_add.rds"))
    d <- fit$draws(format = "df")
    dg <- fitg$draws(format = "df")
    list(
      obs = b$df |> group_by(study, ii, aa, age) |>
        summarise(prop = mean(produces), .groups = "drop") |>
        mutate(study = factor(study, levels = names(STUDY_PAL))),
      params = list(
        sigma_r     = median(d$sigma_r),
        sigma_alpha = median(d$sigma_alpha),
        sigma_zeta_baseline = median(d$sigma_zeta),
        sigma_zeta_add = median(dg$sigma_zeta),
        gamma_add   = median(dg$gamma),
        kappa_pop   = median(d$delta) + 1
      )
    )
  })

  obs  <- s$obs
  p    <- s$params
  a0   <- 15
  ages <- seq(8, 30, length.out = 60)
  L    <- log(ages / a0)

  pA <- ggplot(obs, aes(age, prop, color = study)) +
    geom_point(alpha = 0.35, size = 0.6) +
    geom_smooth(se = FALSE, method = "loess", span = 0.7,
                linewidth = 0.7, formula = y ~ x) +
    scale_color_manual(values = STUDY_PAL, name = "study") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "age (months)", y = "vocabulary proportion produced",
         title = "(A) IO datasets: observed growth") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 10))

  # Variance partition (additive form, at typical mean child)
  parts <- data.frame(
    component = factor(c("Input: intercept", "Input: slope",
                         "Efficiency: intercept", "Efficiency: slope"),
                       levels = c("Efficiency: slope", "Input: slope",
                                  "Efficiency: intercept", "Input: intercept")),
    age = rep(ages, each = 4),
    L   = rep(L, each = 4)
  )
  parts$var <- with(parts, {
    sr2 <- p$sigma_r^2; sa2 <- p$sigma_alpha^2
    sz2 <- p$sigma_zeta_add^2; g <- p$gamma_add
    ifelse(component == "Input: intercept", sr2,
    ifelse(component == "Input: slope",     L^2 * g^2 * sr2 + 2*L*g*sr2,
    ifelse(component == "Efficiency: intercept", sa2,
                                            L^2 * sz2)))
  })
  parts$var <- pmax(parts$var, 0)  # cross can be slightly negative early; clip

  fill_pal <- c("Input: intercept"      = "#B4D7E8",
                "Input: slope"          = "#1F77B4",
                "Efficiency: intercept" = "#003D7A",
                "Efficiency: slope"     = "#E8A0A8")

  pB <- ggplot(parts, aes(age, var, fill = component)) +
    geom_area(alpha = 0.95) +
    scale_fill_manual(values = fill_pal, name = NULL) +
    labs(x = "age (months)",
         y = expression(Var[i](theta[it])),
         title = "(B) Variance partition") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 8))

  pA + pB + plot_layout(widths = c(1, 1))
}

# ---- Figure 6: LLM acceleration --------------------------------
# Two-panel placeholder using existing PNGs. Author iterates.
fig6_panels <- function() {
  list(
    A = repo_path("outputs", "figs", "longitudinal",
                  "chang_bergen_slope_comparison.png"),
    B = repo_path("outputs", "figs", "schematic", "D1_scaling_disanalogy.png")
  )
}
