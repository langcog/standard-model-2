## Build the cache for the revised 2-panel io figure (fig-io-imputed-proc):
##   Panel A "io-imputed": sigma_r sensitivity of the implied input-efficiency
##     share for EN/NO D fits (analytic curve sigma_r^2/sigma_xi^2) + the real
##     refit anchors (EN/NO at pinned sigma_r) with 95% CIs; sigma_r band +
##     meta band.
##   Panel B "io-proc": the joint input+processing variance partition (D'3,
##     all channels free) -- share of efficiency (xi) and acceleration (kappa)
##     variance carried by input vs processing, with 95% CIs.
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

## ---- Panel B: joint io-proc partition (D'3) ----
d3 <- as.data.frame(readRDS(file.path(SUMM, "joint_io_proc_d3.summary.rds")))
gv <- function(v) d3[d3$variable == v, ]
dr <- as_draws_df(readRDS(file.path(SUMM, "joint_io_proc_d3.draws.rds")))
vk <- dr$var_input_k + dr$var_proc_k + dr$var_resid_k          # total kappa variance/draw
q3 <- function(x) c(med = median(x), lo = quantile(x, .05, names = FALSE),
                    hi = quantile(x, .95, names = FALSE))
si_k <- q3(dr$var_input_k / vk); sp_k <- q3(dr$var_proc_k / vk)
panelB <- tibble::tribble(
  ~channel,       ~factor,      ~med,                     ~lo,                    ~hi,
  "Efficiency",   "Input",      gv("share_input_xi")$mean, gv("share_input_xi")$q5, gv("share_input_xi")$q95,
  "Efficiency",   "Processing", gv("share_proc_xi")$mean,  gv("share_proc_xi")$q5,  gv("share_proc_xi")$q95,
  "Acceleration", "Input",      si_k["med"],               si_k["lo"],              si_k["hi"],
  "Acceleration", "Processing", sp_k["med"],               sp_k["lo"],              sp_k["hi"]) |>
  mutate(channel = factor(channel, levels = c("Efficiency", "Acceleration")),
         factor  = factor(factor,  levels = c("Input", "Processing")))

saveRDS(list(panelA = panelA, panelB = panelB),
        file.path(CACHE, "fig_io_imputed_proc.rds"))
cat(sprintf("Wrote fig_io_imputed_proc.rds\n  Panel A anchors: %d (%s)\n",
            nrow(panelA_anchors), paste(unique(panelA_anchors$model), collapse=", ")))
cat("  Panel B (D'3) partition:\n"); print(as.data.frame(panelB), digits = 2, row.names = FALSE)
