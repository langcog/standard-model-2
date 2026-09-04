## Demographics figure (two column): sex and maternal-education effects on
## efficiency and acceleration, per language, cross-sectional + longitudinal,
## with the meta-analytic summary panel. Content unchanged from the
## manuscript's fig-demographics chunk; only the output path/size differ.
##
## Reads: studies/cross_sectional_demographics/cache/fits.rds (committed;
##        xsec fits + meta, and the longitudinal BLUP effects rebuilt from
##        paper/cache/blups_demographics.rds by 00_build.R's last section).
## Run:   Rscript paper/figures/fig_demographics.R
source(here::here("paper", "figures", "_common.R"))
source(here("studies", "cross_sectional_demographics", "composite_figure.R"))

xsec_fits <- readRDS(here("studies", "cross_sectional_demographics", "cache", "fits.rds"))
p <- make_demographics_composite(xsec_fits)

## Tall stacked forests: PNAS caps height at 22.5 cm, so this runs at the
## full double-column width and the maximum height (base_size 9 inside
## composite_figure.R keeps type >= 6 pt at this scale).
save_fig(p, "demographics", width = PNAS_2COL, height = PNAS_MAXH)
