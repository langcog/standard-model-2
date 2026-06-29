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
sx <- function(tag) {
  s <- as.data.frame(readRDS(file.path(SUMM, paste0(tag, ".summary.rds"))))
  s$mean[s$variable == "sigma_xi"]
}
sigma_xi <- c(`English (D)` = sx("long_no_freq_slopes"),
              `Norwegian (D)` = sx("long_no_freq_slopes_norwegian"))
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
panelA_anchors <- bind_rows(
  anchor("long_no_freq_slopes_sigmaR_0p44", 0.44, "English (D)"),
  anchor("long_no_freq_slopes",             0.53, "English (D)"),
  anchor("long_no_freq_slopes_sigmaR_0p58", 0.58, "English (D)"),
  anchor("long_no_freq_slopes_norwegian_sigmaR_0p44", 0.44, "Norwegian (D)"),
  anchor("long_no_freq_slopes_norwegian",             0.53, "Norwegian (D)"),
  anchor("long_no_freq_slopes_norwegian_sigmaR_0p58", 0.58, "Norwegian (D)"))

## sigma_r band = channel-matched apples-to-apples per-sample range; dashed at 0.44.
panelA <- list(curves = panelA_curves, anchors = panelA_anchors,
               sr_band = c(0.36, 0.58), sr_main = 0.44, meta = c(0.04, 0.07))

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
