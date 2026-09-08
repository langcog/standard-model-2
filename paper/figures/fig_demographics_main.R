## Demographics, main-text prototype (1.5 column, single panel): ONE grid --
## columns = trait (efficiency | acceleration, tinted strips), row facets =
## predictor (sex; maternal education). Within each facet: the longitudinal
## datasets by name (open diamonds, raw-BLUP regressions), then the
## longitudinal meta-analysis (filled diamond) and the cross-sectional
## meta-analysis over the displayed Wordbank languages (filled circle, its 95%
## CI shaded). The cross-sectional by-language forests are
## fig_S_demographics_xsec.R. Prototype: does not replace out/demographics.pdf.
##
## Reads: studies/cross_sectional_demographics/cache/fits.rds (committed).
## Run:   Rscript paper/figures/fig_demographics_main.R
source(here::here("paper", "figures", "_common.R"))
source(here("studies", "cross_sectional_demographics", "composite_figure.R"))

xsec_fits <- readRDS(here("studies", "cross_sectional_demographics", "cache", "fits.rds"))
p <- make_demographics_main(xsec_fits, pal = TRAIT_PAL, base_theme = theme_fig(),
                            strips = trait_strips())
save_fig(p, "demographics_main", width = PNAS_15COL, height = 3.3)
