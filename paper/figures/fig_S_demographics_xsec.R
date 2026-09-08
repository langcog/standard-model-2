## SI figure: the cross-sectional by-language forests -- (A) sex, (B) maternal
## education -- per Wordbank language, with the longitudinal estimates as open
## diamonds where a language has them and the cross-sectional meta-analytic 95%
## CI shaded. Moved out of the main demographics figure (fig_demographics_main.R).
##
## Reads: studies/cross_sectional_demographics/cache/fits.rds (committed).
## Run:   Rscript paper/figures/fig_S_demographics_xsec.R
source(here::here("paper", "figures", "_common.R"))
source(here("studies", "cross_sectional_demographics", "composite_figure.R"))

xsec_fits <- readRDS(here("studies", "cross_sectional_demographics", "cache", "fits.rds"))
p <- make_demographics_composite(xsec_fits, pal = TRAIT_PAL, base_theme = theme_fig(),
                                 strips = trait_strips(), forests_only = TRUE)
save_fig(tag_panels(p), "S_demographics_xsec", width = PNAS_2COL, height = 5.2)
