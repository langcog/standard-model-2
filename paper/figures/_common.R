## Shared setup for the standalone figure scripts in paper/figures/.
## Every fig_*.R sources this, reads ONLY committed caches (paper/cache/ and
## the committed study caches), and writes out/<name>.pdf (vector, what goes
## to Overleaf) + out/<name>.png (preview). No numbering: figure order is
## decided in the manuscript, not here.
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(ggh4x)
})

CACHE   <- here("paper", "cache")
FIG_OUT <- here("paper", "figures", "out")
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

## PNAS figure widths (inches): single column 8.7 cm, 1.5 column 11.4 cm,
## double column 17.8 cm; max height 22.5 cm. Figures are saved at their
## final width, so the type sizes below are the printed sizes.
PNAS_1COL  <- 8.7  / 2.54
PNAS_15COL <- 11.4 / 2.54
PNAS_2COL  <- 17.8 / 2.54
PNAS_MAXH  <- 22.5 / 2.54

## ---- One palette, one meaning per hue, across every figure ----
## The two child-level traits and the two predictor channels keep their hues
## wherever they appear (graphical model, demographic forests, io-proc
## schematic + coefficients, data-figure rows, imputed-share curves); the
## vocabulary outcome is neutral. Okabe-Ito blue / vermillion / green plus
## Tol's purple: every pair passes the deutan, tritan and normal-vision
## separation checks and 3:1 contrast on white (dataviz validate_palette.js).
## fig_graphical_model.tex carries the same four hex values.
TRAIT_PAL   <- c(Efficiency = "#0072B2", Acceleration = "#D55E00")
CHANNEL_PAL <- c(Input = "#009E73", Processing = "#AA3377")
VOCAB_COL   <- "#404040"
tint <- function(col, f = 0.2) scales::alpha(col, f)      # light fills (strips, bands)

## Facet strips for the trait columns (efficiency | acceleration), tinted
## with the trait hues so the outcome is colour-coded the same way everywhere.
trait_strips <- function(levels = names(TRAIT_PAL)) {
  ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(
    fill = unname(tint(TRAIT_PAL[levels], 0.22)), colour = NA))
}

## ---- Type: the same sizes in every figure ----
## axis text 7 pt, axis titles / strips / panel titles 8 pt, panel tags 10 pt
## bold, in-plot annotations 7 pt (TXT is that size in geom_text's mm units).
TXT <- 7 / .pt
theme_fig <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(panel.border     = element_rect(colour = "grey55", fill = NA, linewidth = 0.3),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = 7, colour = "grey30"),
          axis.title       = element_text(size = 8),
          strip.text       = element_text(size = 8, face = "bold", margin = margin(2, 2, 2, 2)),
          strip.background = element_rect(fill = "grey93", colour = NA),
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 8),
          legend.key.size  = unit(0.7, "lines"),
          plot.title       = element_text(size = 8, face = "bold", hjust = 0),
          plot.tag         = element_text(size = 10, face = "bold"),
          plot.margin      = margin(4, 6, 3, 3))
}

## Letter the panels of a composed patchwork (A, B, ...) at the top-left of
## each panel's own box. Single-panel figures are not tagged.
tag_panels <- function(p) p + plot_annotation(tag_levels = "A")

## Save vector PDF (cairo, for non-ASCII glyphs like theta/kappa) + PNG preview.
save_fig <- function(p, name, width, height, dpi = 300) {
  pdf_path <- file.path(FIG_OUT, paste0(name, ".pdf"))
  ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(FIG_OUT, paste0(name, ".png")), p, width = width, height = height,
         dpi = dpi, bg = "white")
  cat(sprintf("wrote %s (%.2f x %.2f in)\n", pdf_path, width, height))
  invisible(pdf_path)
}
