## Shared setup for the standalone figure scripts in paper/figures/.
## Every fig_*.R sources this, reads ONLY committed caches (paper/cache/ and
## the committed study caches), and writes out/<name>.pdf (vector, what goes
## to Overleaf) + out/<name>.png (preview). No numbering: figure order is
## decided in the manuscript, not here.
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

CACHE   <- here("paper", "cache")
FIG_OUT <- here("paper", "figures", "out")
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

## PNAS figure widths (inches): single column 8.7 cm, 1.5 column 11.4 cm,
## double column 17.8 cm; max height 22.5 cm. Type must stay >= 6 pt at
## final size, so base_size 8 (2-col) / 7-8 (1-col) in ggplot.
PNAS_1COL  <- 8.7  / 2.54
PNAS_15COL <- 11.4 / 2.54
PNAS_2COL  <- 17.8 / 2.54
PNAS_MAXH  <- 22.5 / 2.54

## Palettes (Okabe-Ito; the same hues the manuscript chunks used)
FACTOR_PAL <- c(Input = "#009E73", Processing = "#CC79A7")      # io-proc channels
MODEL_PAL  <- c(`English (D)` = "#E69F00", `Norwegian (D)` = "#56B4E9")
KID_PAL    <- c("lower-ability child" = "#1b9e77", "higher-ability child" = "#7570b3")

theme_fig <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = base_size + 1),
          strip.text = element_text(face = "bold"),
          legend.key.size = unit(0.8, "lines"))
}

## Save vector PDF (cairo, for non-ASCII glyphs like theta/kappa) + PNG preview.
save_fig <- function(p, name, width, height, dpi = 300) {
  pdf_path <- file.path(FIG_OUT, paste0(name, ".pdf"))
  ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(FIG_OUT, paste0(name, ".png")), p, width = width, height = height,
         dpi = dpi, bg = "white")
  cat(sprintf("wrote %s (%.2f x %.2f in)\n", pdf_path, width, height))
  invisible(pdf_path)
}
