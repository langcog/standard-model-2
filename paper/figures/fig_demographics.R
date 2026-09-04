## Demographics figure (two column): sex and maternal-education effects on
## efficiency and acceleration, per language, cross-sectional + longitudinal,
## with the meta-analytic summary panel. Layout: (A) sex forest, tall, on the
## left; (B) maternal-education forest over (C) the meta summary on the right.
## Colour = trait (TRAIT_PAL, tinted strips), shape = method.
##
## Reads: studies/cross_sectional_demographics/cache/fits.rds (committed;
##        xsec fits + meta, and the longitudinal BLUP effects rebuilt from
##        paper/cache/blups_demographics.rds by 00_build.R's last section).
## Run:   Rscript paper/figures/fig_demographics.R
source(here::here("paper", "figures", "_common.R"))
source(here("studies", "cross_sectional_demographics", "composite_figure.R"))

xsec_fits <- readRDS(here("studies", "cross_sectional_demographics", "cache", "fits.rds"))
p <- make_demographics_composite(xsec_fits, pal = TRAIT_PAL, base_theme = theme_fig(),
                                 strips = trait_strips())
save_fig(tag_panels(p), "demographics", width = PNAS_2COL, height = 5.6)
