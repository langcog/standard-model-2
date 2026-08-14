## _helpers.R — thin shared utilities sourced by standard_model.qmd.
## Most plotting and data manipulation lives inline in the chunks; this
## file just holds palettes, label maps, and the Figure 1 schematic
## (kept here because the panel layout is verbose).

# ---- palettes & labels ------------------------------------------------
# By-study units (English split into Thal/Smith/Marchman) plus the legacy
# whole-language keys (kept for any non-ladder reference).
LANG_PAL <- c("thal"             = "#E69F00",
              "smith"            = "#009E73",
              "marchman"         = "#CC79A7",
              "english_american" = "#E69F00",
              "norwegian"        = "#56B4E9",
              "french_quebecois" = "#009E73",
              "japanese"         = "#D55E00",
              "finnish"          = "#CC79A7",
              "english_british"  = "#999999")

LANG_LABELS <- c("thal"             = "English (Thal)",
                 "smith"            = "English (Smith)",
                 "marchman"         = "English (Marchman)",
                 "english_american" = "English (American)",
                 "norwegian"        = "Norwegian",
                 "french_quebecois" = "French (Quebecois)",
                 "japanese"         = "Japanese",
                 "finnish"          = "Finnish",
                 "english_british"  = "English (British)")

STUDY_PAL <- c("BabyView"  = "#E69F00",
               "SEEDLingS" = "#56B4E9",
               "AM2018"    = "#009E73",
               "FMW2013"   = "#D55E00")

# Wordbank quantile palette (matches glmer_ladder/04b_plot.R) — used in
# Fig 1 schematic + Fig 2 model ladder so the 10/25/50/75/90 quantile
# fans line up visually across the manuscript.
# Okabe-Ito. The previous ramp paired green (50th) with red (90th), which Science
# prohibits outright ("avoid using red and green together") and which measured as the
# worst pair in the figure: under simulated deuteranopia those two sat at dE 13.8 in CIE
# Lab, where anything under ~15 starts to be confusable. This ordering keeps the same
# low-to-high blue-to-warm direction and lifts the worst case to 17.6 with no green.
# Ordered ramps that score higher (RdYlBu 28.4, viridis 18.2) all place a near-white
# colour at one end -- in RdYlBu that lands on the emphasised median.
WORDBANK_PALETTE <- c("0.1"  = "#0072B2",   # blue
                      "0.25" = "#56B4E9",   # sky blue
                      "0.5"  = "#CC79A7",   # reddish purple
                      "0.75" = "#E69F00",   # orange
                      "0.9"  = "#D55E00")   # vermillion

# The five main-text by-study units (English split into its 3 datasets)
PAPER_LANGS <- c("thal", "smith", "marchman", "norwegian", "japanese")

# ---- Figure 1 schematic ----------------------------------------------
# Stand-alone illustrative ggplot — no data needed; pure model display.
make_fig1_schematic <- function() {
  age <- seq(8, 30, length.out = 80)
  a0  <- 19
  L   <- log(age / a0)

  variants <- list(
    list(name = "(A) pure accumulator", xi = 0, kappa = 1,
         eq = r"($\theta = \log(r) + \log(t)$)"),
    list(name = "(B) + efficiency variance", xi = c(-1.5, 0, 1.5), kappa = 1,
         eq = r"($\theta = \log(\alpha_i) + \log(r) + \log(t)$)"),
    list(name = "(C) + acceleration", xi = c(-1.5, 0, 1.5), kappa = 2.5,
         eq = r"($\theta = \log(\alpha_i) + \log(r) + \kappa\,\log(t)$)"),
    list(name = "(D) + per-kid acceleration", xi = c(-1.5, 0, 1.5),
         kappa = c(1.5, 2.5, 3.5),
         eq = r"($\theta = \log(\alpha_i) + \log(r) + \kappa_i\,\log(t)$)"),
    list(name = "(E) + per-kid input variance", xi = c(-1.5, 0, 1.5),
         kappa = c(1.5, 2.5, 3.5), input = c(-0.6, 0, 0.6),
         eq = r"($\theta = \log(\alpha_i) + \log(r_i) + \kappa_i\,\log(t)$)")
  )

  build_panel <- function(v) {
    xi  <- if (length(v$xi)    == 1) rep(v$xi,    3) else v$xi
    kap <- if (length(v$kappa) == 1) rep(v$kappa, 3) else v$kappa
    inp <- if (is.null(v$input)) rep(0, 3) else v$input
    qd  <- bind_rows(lapply(seq_along(xi), function(i) {
      data.frame(quantile = c("25%", "50%", "75%")[i], age = age,
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

  patchwork::wrap_plots(lapply(variants, build_panel), ncol = 5) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}
