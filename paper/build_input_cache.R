## build_input_cache.R — data for the input-experiments triptych (Fig 3).
## Produces paper/cache/fig3_input.rds with three panels, all on one
## y-axis ("input share of between-child variance, 1 - pi"):
##   (A) input share vs AGE (continuous; intercept channel)
##   (B) ANALYTIC sensitivity of the input share to the sigma_r pin, with
##       refit anchors overlaid (see Supplemental Methods derivation)
##   (C) input share per model (EN D, NO D, IO wide-delta) vs the
##       meta-analytic band
##
## INTERCEPT CHANNEL ONLY for now: the slope/gamma channel (pi_zeta) needs
## the D' fits. EN/NO read the post-dedup D summaries; NO is currently the
## OLD pre-dedup summary as a placeholder until the salvage lands. Rebuild
## when new NO D + the D' fits arrive.

suppressPackageStartupMessages({ library(dplyr); library(here); library(tibble) })
CACHE <- here("paper", "cache"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

## ---- model scalars from the slim summaries ----
scal <- function(tag) {
  s <- as.data.frame(readRDS(here("fits", "summaries", paste0(tag, ".summary.rds"))))
  g <- function(v, col = "median") {
    if (!col %in% names(s)) return(NA_real_)
    i <- which(s$variable == v); if (length(i)) s[[col]][i] else NA_real_
  }
  sx <- g("sigma_xi"); sa <- g("sigma_alpha")
  list(sigma_xi = sx, sigma_alpha = sa, sigma_zeta = g("sigma_zeta"),
       rho = g("rho_xi_zeta"), sigma_r = sqrt(max(sx^2 - sa^2, 1e-6)),
       pi_alpha = g("pi_alpha"), pi_lo = g("pi_alpha", "q025"),
       pi_hi = g("pi_alpha", "q975"))
}
en <- scal("long_no_freq_slopes")             # post-dedup EN D
no <- scal("long_no_freq_slopes_norwegian")   # OLD pre-dedup placeholder

## ---- IO wide-delta (from the fit draws) ----
iofit <- readRDS(here("fits", "io_pooled_widedelta.rds"))
iod <- iofit$draws(format = "df")
io_pi <- iod$sigma_alpha^2 / (iod$sigma_alpha^2 + iod$sigma_r^2)
io <- list(sigma_xi = sqrt(median(iod$sigma_alpha)^2 + median(iod$sigma_r)^2),
           sigma_alpha = median(iod$sigma_alpha), sigma_zeta = median(iod$sigma_zeta),
           rho = 0, sigma_r = median(iod$sigma_r),
           pi_alpha = median(io_pi),
           pi_lo = unname(quantile(io_pi, 0.025)),
           pi_hi = unname(quantile(io_pi, 0.975)))

MODELS <- list(`English (D)` = en, `Norwegian (D)` = no,
               `Input-observed` = io)
MODEL_LEVELS <- names(MODELS)

## ---- Panel A: input share vs age (intercept channel) ----
## Var(theta_i(t)) = sigma_xi^2 + 2 rho sigma_xi sigma_zeta L + sigma_zeta^2 L^2
## input variance (intercept channel) = sigma_r^2.  share = sigma_r^2 / Var.
a0 <- 18; ages <- seq(8, 30, by = 0.5); L <- log(ages / a0)
panelA <- bind_rows(lapply(names(MODELS), function(nm) {
  p <- MODELS[[nm]]
  vtheta <- p$sigma_xi^2 + 2 * p$rho * p$sigma_xi * p$sigma_zeta * L + p$sigma_zeta^2 * L^2
  tibble(model = nm, age = ages, input_share = p$sigma_r^2 / vtheta)
})) |> mutate(model = factor(model, levels = MODEL_LEVELS))

## ---- Panel B: analytic sensitivity (EN, NO) + refit anchors ----
## Holding sigma_xi^2 at its fitted (data-determined) value, the input
## share is the algebraic split  sigma_r^2 / sigma_xi^2.  Refit anchors
## show the small second-order drift (sigma_xi^2 shrinks slightly as the
## pin grows) -- see Supplemental Methods.
sr_grid <- seq(0.2, 1.0, by = 0.01)
panelB_curve <- bind_rows(lapply(c("English (D)", "Norwegian (D)"), function(nm) {
  p <- MODELS[[nm]]
  tibble(model = nm, sigma_r = sr_grid, input_share = sr_grid^2 / p$sigma_xi^2)
})) |> mutate(model = factor(model, levels = MODEL_LEVELS))
anch <- function(tag, sr, model, note) {
  fp <- here("fits", "summaries", paste0(tag, ".summary.rds"))
  if (!file.exists(fp)) return(NULL)
  s <- as.data.frame(readRDS(fp)); pa <- s$median[s$variable == "pi_alpha"]
  tibble(model = model, sigma_r = sr, input_share = 1 - pa, note = note)
}
panelB_anchors <- bind_rows(
  anch("long_no_freq_slopes", en$sigma_r, "English (D)", "post-dedup (on pin)"),
  anch("long_no_freq_slopes_sigmaR_0p80", 0.80, "English (D)", "pre-dedup refit"))

## ---- Panel C: input share per model + meta band ----
panelC <- bind_rows(lapply(names(MODELS), function(nm) {
  p <- MODELS[[nm]]
  tibble(source = nm, kind = "model",
         input_share = 1 - p$pi_alpha, lo = 1 - p$pi_hi, hi = 1 - p$pi_lo)
}))
## Meta-analytic input-quantity share of vocabulary variance, Coffey 2026:
## 4-7%. (Anderson 2021 give a "modest proportion" -- MCF to fill the number.)
meta <- tibble(source = "Coffey 2026 (meta)", kind = "meta",
               input_share = 0.055, lo = 0.04, hi = 0.07)
panelC <- bind_rows(panelC, meta)

saveRDS(list(panelA = panelA, panelB_curve = panelB_curve,
             panelB_anchors = panelB_anchors, panelC = panelC,
             sr_main = 0.53, sr_range = c(0.40, 0.70),  # plausible band: MCF to set
             a0 = a0, intercept_only = TRUE,
             no_is_placeholder = TRUE),
        file.path(CACHE, "fig3_input.rds"))
cat("Wrote paper/cache/fig3_input.rds\n\n")
cat("Panel C (input share 1 - pi_alpha):\n"); print(as.data.frame(panelC))
cat(sprintf("\nPanel A input share at a0 (age 18): EN=%.3f NO=%.3f IO=%.3f\n",
            panelA$input_share[panelA$model=="English (D)" & panelA$age==18],
            panelA$input_share[panelA$model=="Norwegian (D)" & panelA$age==18],
            panelA$input_share[panelA$model=="Input-observed" & panelA$age==18]))
