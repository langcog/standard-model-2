## _setup_shared.R -- everything both the main manuscript and the standalone
## supplement need: packages, caches, palettes and the inline values quoted in prose.
##
## Extracted from the main .qmd setup chunk so the supplement can be rendered on its own.
## Science requires the supplementary materials as a SEPARATE file, and the SI cannot be
## built alone while the objects it quotes (n_excl, pct_excl, qc_mar, qc_nor and the
## caches) exist only inside the main document's setup chunk -- that is why rendering
## supplemental.qmd by itself has always failed with "object 'n_excl' not found".
##
## Sourced by both documents. Do not duplicate its contents into either.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork); library(latex2exp)
  library(knitr); library(scales); library(stringr); library(forcats)
  # png/grid dropped: they were used only by the plate-diagram raster removed from Fig 1.
  # [MCF: gridExtra appears unused too, but it predates this change so it is left alone.]
  library(broom); library(gridExtra)
  library(ggrepel)
})
source(here("paper", "_helpers.R"))



# Cached data for the three main-text figures — produced by
# `Rscript paper/build_cache_short.R`.
CACHE <- here("paper", "cache")
fig1_fan   <- readRDS(file.path(CACHE, "fig1_fan.rds"))        # Fig 1: M0-M3 schematic + per-dataset M3 fan
fig2       <- readRDS(file.path(CACHE, "fig2_efficiency.rds")) # Fig 2: per-word exposures-to-learn (M3)
llm_slopes <- readRDS(file.path(CACHE, "fig6_llm_slopes.rds")) # Fig 3: child vs LLM acceleration
# QC cache
TR <- readRDS(file.path(CACHE, "qc_spaghetti_data.rds"))
# for percent CV 
si_cvh    <- readRDS(here("paper", "cache", "si_cv_headline.rds"))  # forward CV: M2 vs kappa=1


# Inline text values (raw numbers in si_inline.rds, built by build_cache_short
# section 7); formatted here to the precision they render at in the text.
inl <- readRDS(file.path(CACHE, "si_inline.rds"))
.k1  <- function(x) sprintf("%.1f", x)                        # kappa: one decimal
.pct <- function(x) if (x >= 10) sprintf("%.0f", x) else sprintf("%.1f", x)
age_lo   <- inl$age_lo; age_hi <- inl$age_hi                  # integer months
loo_min  <- format(signif(inl$loo_min, 2), big.mark = ",")   # 2 sig figs
kap_lo   <- .k1(inl$kappa_lo); kap_lo_q5 <- .k1(inl$kappa_lo_q5); kap_lo_q95 <- .k1(inl$kappa_lo_q95)
kap_hi   <- .k1(inl$kappa_hi); kap_hi_q5 <- .k1(inl$kappa_hi_q5); kap_hi_q95 <- .k1(inl$kappa_hi_q95)
en_kappa <- .k1(inl$en_kappa); no_kappa <- .k1(inl$no_kappa)
en_sd    <- .k1(inl$en_sd);    no_sd    <- .k1(inl$no_sd)

# Random-effect summaries for main-text reporting (see tbl-ranef). sigma_b = between-
# child SD of acceleration; rho = efficiency-acceleration correlation; pct_gt1 = share
# of children above the pure-accumulator value kappa = 1; kappa_between_sd = SD of
# kappa ACROSS the five samples, the honest uncertainty scale given that CDI items are
# not locally independent (it is ~13x the narrowest within-sample interval).
sb_lo    <- .k1(inl$sb_lo);  sb_hi  <- .k1(inl$sb_hi)
rho_lo   <- sprintf("%+.2f", inl$rho_lo); rho_hi <- sprintf("%+.2f", inl$rho_hi)
pct_gt1  <- sprintf("%.1f", inl$pct_kappa_gt1)
kappa_between_sd <- .k1(inl$kappa_between_sd)

# Longitudinal sample + QC exclusions (bayes_long_sample.rds; built from the bundles by
# studies/bayes_long/qc_exclusion_report.R, so it renders independently of the fit caches).
bl_samp  <- readRDS(file.path(CACHE, "bayes_long_sample.rds"))
.bl      <- function(s) bl_samp$sample[bl_samp$sample$dataset == s, ]
qc_mar   <- .pct(.bl("marchman")$pct_kids_excluded)   # % children excluded, Marchman
qc_nor   <- .pct(.bl("norwegian")$pct_kids_excluded)  # % children excluded, Norwegian
n_total  <- format(.bl("total")$n_kids, big.mark = ",")      # kids retained, all datasets
n_excl   <- .bl("total")$kids_excluded                       # kids excluded, all datasets
pct_excl <- .pct(.bl("total")$pct_kids_excluded)

# dev = "cairo_pdf": the default pdf() device cannot encode non-ASCII glyphs. It fails
# with an mbcsToSbcs conversion error and substitutes ".", so theta, kappa and the
# up-arrow were rendering as "Latent Ability (.)", "Pure (.=1)" and ". Efficiency" in
# Fig 1B (and kappa/gamma in the LM figure). cairo_pdf handles UTF-8; it needs XQuartz
# for its X11 libraries, which is now installed. Verify with capabilities("cairo").
# The device has to follow the output format. cairo_pdf emits vector PDF, which is right
# for the LaTeX build but which Word cannot render -- a docx built with it embeds three
# .pdf files that open as blank placeholders. ragg_png handles UTF-8 as well as cairo does
# and satisfies Science's 300 dpi requirement for initial submission.
knitr::opts_chunk$set(
  echo = FALSE, message = FALSE, warning = FALSE, cache=TRUE,
  fig.align = "center", out.width = "100%", dpi = 300,
  dev = if (knitr::is_latex_output()) "cairo_pdf" else "ragg_png"
)
