## Build the cache for the revised 2-panel io figure (fig-io-imputed-proc):
##   Panel A "io-imputed": sigma_r sensitivity of the implied input-efficiency
##     share for EN/NO D fits (analytic curve sigma_r^2/sigma_xi^2) + the real
##     refit anchors (EN/NO at pinned sigma_r) with 95% CIs; sigma_r band +
##     meta band.
##   Panel B "io-proc": the joint input+processing COEFFICIENTS (per-SD effects,
##     the EN+count +proc fit joint_io_proc_lean_d2_enct -- all channels free) --
##     input vs processing effect on efficiency (xi) and acceleration (kappa),
##     with 90% CIs. (Variance-explained partition is reported in text instead.)
##
## RUN LOCALLY (needs the gitignored fits/summaries). Writes the committed
## cache paper/cache/fig_io_imputed_proc.rds so the manuscript renders without
## the raw fits.
suppressPackageStartupMessages({ library(here); library(dplyr); library(posterior) })
SUMM <- here("fits", "summaries")
CACHE <- here("paper", "cache")

## ---- Panel A: io-imputed sigma_r sensitivity ----
## sigma_xi (total efficiency SD) is data-identified and stable across the
## sigma_r pins, so read it from the 0.44-pin fits; the un-suffixed baseline
## (bundle pin 0.53, Sperry-era) is retired -- do NOT reference it, or a
## stale pre-cleanup summary on disk will silently feed the curves.
sx <- function(tag) {
  s <- as.data.frame(readRDS(file.path(SUMM, paste0(tag, ".summary.rds"))))
  s$mean[s$variable == "sigma_xi"]
}
## Tag lineage: the June anchors were tagged long_no_freq_slopes_* (the GCP
## runner passed the variant with its long_ prefix); the 2026-08/09 clean-data
## convergence refits are the SAME variant (variant_hyperpriors() strips the
## prefix) run as no_freq_slopes with STAN_TAG_SFX=_2k, 2000/1000 iters,
## adapt_delta 0.97. Read those; the long_* files on disk are pre-cleanup.
tag2k <- function(lang, sr) sprintf("no_freq_slopes%s_sigmaR_%s_2k",
                                    if (lang == "no") "_norwegian" else "", sr)
sigma_xi <- c(`English (D)` = sx(tag2k("en", "0p44")),
              `Norwegian (D)` = sx(tag2k("no", "0p44")))
sr_grid <- seq(0.25, 0.80, by = 0.005)
panelA_curves <- bind_rows(lapply(names(sigma_xi), function(m)
  data.frame(model = m, sigma_r = sr_grid, share = sr_grid^2 / sigma_xi[[m]]^2)))

## real refit anchors: share = 1 - pi_alpha (mean + 95% CI). Include each only
## if its summary exists (NO sigma_r pins arrive later -> re-run to add them).
anchor <- function(tag, sr, model) {
  fp <- file.path(SUMM, paste0(tag, ".summary.rds"))
  if (!file.exists(fp)) return(NULL)
  s <- as.data.frame(readRDS(fp)); r <- s[s$variable == "pi_alpha", ]
  data.frame(model = model, sigma_r = sr,
             share = 1 - r$mean, lo = 1 - r$q975, hi = 1 - r$q025)
}
## Anchor set (2026-08-15 decision): {0.35, 0.44, 0.58}.
##   0.35 = the DIRECT joint/io fits' own sigma_r estimate on observed
##          recordings (enct/io_count: 0.35 [0.32, 0.39]; study-centered, so
##          a lower bound on population sigma_r);
##   0.44 = channel-matched all-adult population re-anchor (main pin);
##   0.58 = literature upper. 0.53 (Sperry CDS) retired as legacy.
panelA_anchors <- bind_rows(
  anchor(tag2k("en", "0p35"), 0.35, "English (D)"),
  anchor(tag2k("en", "0p44"), 0.44, "English (D)"),
  anchor(tag2k("en", "0p58"), 0.58, "English (D)"),
  anchor(tag2k("no", "0p35"), 0.35, "Norwegian (D)"),
  anchor(tag2k("no", "0p44"), 0.44, "Norwegian (D)"),
  anchor(tag2k("no", "0p58"), 0.58, "Norwegian (D)"))
stopifnot(nrow(panelA_anchors) == 6)   # all six _2k summaries must be present

## sigma_r band: lower edge = the direct estimate's 95% low (0.32), upper =
## literature (0.58); anchors now span the band. Dashed main pin at 0.44.
panelA <- list(curves = panelA_curves, anchors = panelA_anchors,
               sr_band = c(0.32, 0.58), sr_main = 0.44, meta = c(0.04, 0.07))

## ---- Panel B: joint io-proc COEFFICIENTS (per-SD effects, EN+count +proc fit) ----
## The four channel x trait per-SD effects (log-odds), same numbers the schematic draws:
##   Efficiency:   Input = sigma_r,        Processing = |eff_proc_xi|  (faster -> higher xi)
##   Acceleration: Input = eff_input_k,    Processing = |eff_proc_k|   (faster -> higher kappa)
## Processing flipped to "faster processing -> +" (eff_proc_* are negative on rt0). The qmd
## rescales the acceleration row x log(30/21) to level-equivalent theta units for display.
d3 <- as.data.frame(readRDS(file.path(SUMM, "joint_io_proc_lean_d2_enct_2k.summary.rds")))
gv <- function(v) d3[d3$variable == v, ]
panelB <- tibble::tribble(
  ~channel,       ~factor,      ~med,                    ~lo,                     ~hi,
  "Efficiency",   "Input",       gv("sigma_r")$mean,      gv("sigma_r")$q5,        gv("sigma_r")$q95,
  "Efficiency",   "Processing", -gv("eff_proc_xi")$mean, -gv("eff_proc_xi")$q95,  -gv("eff_proc_xi")$q5,
  "Acceleration", "Input",       gv("eff_input_k")$mean,  gv("eff_input_k")$q5,    gv("eff_input_k")$q95,
  "Acceleration", "Processing", -gv("eff_proc_k")$mean,  -gv("eff_proc_k")$q95,   -gv("eff_proc_k")$q5) |>
  mutate(channel = factor(channel, levels = c("Efficiency", "Acceleration")),
         factor  = factor(factor,  levels = c("Input", "Processing")))

## ---- schematic params (panel A): the SAME EN+count +proc fit, so the trajectories are
## drawn from the estimates rather than hand-set example values. a0/log_H/mu_r are data
## constants from the bundle; everything else is read from the fit summary above.
schem <- list(
  a0 = 21, log_H = 5.90, mu_r = 7.34,
  kpop     = 1 + gv("delta")$mean,        # population acceleration 1 + delta
  sigma_xi = gv("sigma_xi")$mean,         # between-child efficiency SD (curve separation)
  d_in_xi  = gv("sigma_r")$mean,          # input  -> efficiency  (coeff-1 identity = sigma_r)
  d_in_k   = gv("eff_input_k")$mean,      # input  -> acceleration
  d_pr_xi  = -gv("eff_proc_xi")$mean,     # processing(faster) -> efficiency
  d_pr_k   = -gv("eff_proc_k")$mean)      # processing(faster) -> acceleration

saveRDS(list(panelA = panelA, panelB = panelB, schem = schem),
        file.path(CACHE, "fig_io_imputed_proc.rds"))
cat(sprintf("Wrote fig_io_imputed_proc.rds\n  Panel A anchors: %d (%s)\n",
            nrow(panelA_anchors), paste(unique(panelA_anchors$model), collapse=", ")))
cat("  Panel B (EN+count +proc) coefficients:\n"); print(as.data.frame(panelB), digits = 2, row.names = FALSE)
