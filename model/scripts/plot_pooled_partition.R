## Variance-partition + coefficient plot for the pooled IO ladder.
##
## Top row : Stacked variance partition (input vs efficiency / residual)
##           for INTERCEPT and SLOPE (baseline / gamma-add / gamma-mult).
##           This is the r:alpha (intercept) and gamma:zeta (slope) story.
##
## Bottom row: Coefficient plot with 90% CIs for the key parameters that
##             produced the partition above (sigma_r, sigma_alpha, gamma,
##             sigma_zeta).
##
## Output: figs/io/pooled_partition.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

cat("Loading fits ...\n")
fits <- list(
  baseline = readRDS(file.path(PATHS$fits_dir, "io_pooled.rds")),
  add      = readRDS(file.path(PATHS$fits_dir, "io_pooled_gamma_add.rds")),
  mult     = readRDS(file.path(PATHS$fits_dir, "io_pooled_gamma_mult.rds"))
)

# ---- pull draws --------------------------------------------------------
get_draws <- function(fit, gamma = FALSE) {
  d <- fit$draws(format = "df")
  out <- list(
    sigma_r     = d$sigma_r,
    sigma_alpha = d$sigma_alpha,
    sigma_zeta  = d$sigma_zeta,
    delta       = d$delta
  )
  if (gamma) out$gamma <- d$gamma
  # kappa_study to get A (base slope without zeta) — average across studies
  ks <- as.matrix(d[, grep("^kappa_study\\[", colnames(d), value = TRUE), drop = FALSE])
  out$A <- rowMeans(ks)
  out
}
dr <- list(
  baseline = get_draws(fits$baseline, gamma = FALSE),
  add      = get_draws(fits$add,      gamma = TRUE),
  mult     = get_draws(fits$mult,     gamma = TRUE)
)

# ---- per-draw variance partition --------------------------------------
# Intercept: V_input = sigma_r^2, V_eff = sigma_alpha^2 (same across models;
#   use the baseline draws as reference)
# Slope:
#   baseline : V_input = 0, V_eff = sigma_zeta^2
#   add      : V_input = (gamma * sigma_r)^2, V_eff = sigma_zeta^2
#   mult     : V_input = (gamma*sigma_r)^2 * (A^2 + sigma_zeta^2)
#              + an interaction term sigma_zeta^2 * (gamma*sigma_r)^2
#              Total slope variance:
#                Var(slope) = gamma^2 sigma_r^2 (A^2 + sigma_zeta^2) + sigma_zeta^2
#              Following the same r:alpha-style partition we keep V_eff = sigma_zeta^2
#              and call everything else input-driven.

slope_part <- function(d, model) {
  if (model == "baseline") {
    list(V_input = rep(0, length(d$sigma_zeta)),
         V_eff   = d$sigma_zeta^2)
  } else if (model == "add") {
    list(V_input = (d$gamma * d$sigma_r)^2,
         V_eff   = d$sigma_zeta^2)
  } else { # mult
    g2sr2 <- (d$gamma * d$sigma_r)^2
    list(V_input = g2sr2 * (d$A^2 + d$sigma_zeta^2),
         V_eff   = d$sigma_zeta^2)
  }
}

partition_rows <- list()
# intercept (shared structure — baseline draws as canonical)
V_in_int  <- dr$baseline$sigma_r^2
V_ef_int  <- dr$baseline$sigma_alpha^2
share_int <- V_in_int / (V_in_int + V_ef_int)
partition_rows[["intercept"]] <- data.frame(
  channel = "intercept",
  model   = "shared",
  share_lo = quantile(share_int, 0.05),
  share_md = quantile(share_int, 0.50),
  share_hi = quantile(share_int, 0.95),
  V_in_md = median(V_in_int),
  V_ef_md = median(V_ef_int),
  total_md = median(V_in_int + V_ef_int)
)
for (m in c("baseline","add","mult")) {
  sp <- slope_part(dr[[m]], m)
  share <- sp$V_input / (sp$V_input + sp$V_eff)
  partition_rows[[paste0("slope_", m)]] <- data.frame(
    channel = "slope", model = m,
    share_lo = quantile(share, 0.05),
    share_md = quantile(share, 0.50),
    share_hi = quantile(share, 0.95),
    V_in_md = median(sp$V_input),
    V_ef_md = median(sp$V_eff),
    total_md = median(sp$V_input + sp$V_eff)
  )
}
parts <- bind_rows(partition_rows)
print(parts)

# Long form for stacked bars (input share + eff share, summing to 100%)
bars <- parts |>
  mutate(label = ifelse(channel == "intercept", "intercept (shared)",
                  paste0("slope (", model, ")"))) |>
  select(label, V_in_md, V_ef_md, share_md, share_lo, share_hi) |>
  mutate(input_pct = 100 * V_in_md / (V_in_md + V_ef_md),
         eff_pct   = 100 - input_pct)
bars$label <- factor(bars$label, levels = c(
  "intercept (shared)",
  "slope (baseline)", "slope (add)", "slope (mult)"))

# Stacked-bar dataframe
bars_long <- bars |>
  select(label, input_pct, eff_pct) |>
  pivot_longer(c(input_pct, eff_pct), names_to = "component", values_to = "pct") |>
  mutate(component = factor(component, levels = c("eff_pct","input_pct"),
                            labels = c("efficiency / residual", "input-driven")))

# Pre-compute input-share label text (with 90% CI)
bar_labs <- bars |>
  mutate(txt = sprintf("input %.1f%%\n[%.1f, %.1f]",
                       100*share_md, 100*share_lo, 100*share_hi))

# ---- panel A: stacked variance partition ------------------------------
panelA <- ggplot(bars_long, aes(label, pct, fill = component)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(data = bar_labs, aes(label, y = 50, label = txt),
            inherit.aes = FALSE, size = 3.2, lineheight = 0.95) +
  scale_fill_manual(values = c("efficiency / residual" = "#888888",
                               "input-driven"          = "#D55E00")) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
  labs(x = NULL, y = "% of channel variance", fill = NULL,
       title = "Variance partition: input-driven vs efficiency-side",
       subtitle = "intercept (r:α) and slope (γ·log_r_dev : ζ) share the same input latent log_r_dev") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())

# ---- panel B: coefficient plot ----------------------------------------
coef_df <- bind_rows(
  data.frame(param = "σ_r",          model = "shared", draws = dr$baseline$sigma_r),
  data.frame(param = "σ_α",          model = "shared", draws = dr$baseline$sigma_alpha),
  data.frame(param = "σ_ζ",          model = "baseline", draws = dr$baseline$sigma_zeta),
  data.frame(param = "σ_ζ",          model = "add",      draws = dr$add$sigma_zeta),
  data.frame(param = "σ_ζ",          model = "mult",     draws = dr$mult$sigma_zeta),
  data.frame(param = "γ (additive)", model = "add",      draws = dr$add$gamma),
  data.frame(param = "γ·κ (mult, effective)", model = "mult",
             draws = dr$mult$gamma * dr$mult$A)
)
coef_summary <- coef_df |>
  group_by(param, model) |>
  summarise(md = median(draws), lo = quantile(draws, 0.05),
            hi = quantile(draws, 0.95), .groups = "drop") |>
  mutate(label = ifelse(model == "shared", param, paste0(param, " [", model, "]")),
         label = factor(label, levels = rev(c(
           "σ_r", "σ_α",
           "σ_ζ [baseline]", "σ_ζ [add]", "σ_ζ [mult]",
           "γ (additive) [add]", "γ·κ (mult, effective) [mult]"))))

panelB <- ggplot(coef_summary, aes(md, label)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  labs(x = "posterior median (90% CI)", y = NULL,
       title = "Underlying parameter posteriors",
       subtitle = "γ·κ rescales the multiplicative coef to slope-units (comparable to additive γ)") +
  theme_minimal(base_size = 11)

out_path <- file.path("figs", "io", "pooled_partition.png")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(out_path, panelA / panelB + plot_layout(heights = c(1.0, 1.1)),
       width = 11, height = 9, dpi = 150)
cat(sprintf("Wrote %s\n", out_path))
