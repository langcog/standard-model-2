## Imputed population input share vs the meta-analytic range (single column).
## Factored out of the io-proc figure so the population-level imputation is
## not confused with the direct io-proc study: analytic curves
## share = sigma_r^2 / sigma_xi^2 for the EN and NO D fits, the six real
## refit anchors at sigma_r = {0.35, 0.44, 0.58} with 95% CIs, the plausible
## sigma_r band [0.32, 0.58] and the meta-analytic 4-7% band.
##
## Reads: paper/cache/fig_io_imputed_proc.rds (panelA; built by
##        paper/build_fig_io_cache.R from the _2k anchor fits).
## Run:   Rscript paper/figures/fig_imputed_meta.R
source(here::here("paper", "figures", "_common.R"))

io_fig <- readRDS(file.path(CACHE, "fig_io_imputed_proc.rds"))
a <- io_fig$panelA
pct <- scales::percent_format(accuracy = 1)

p <- ggplot() +
  annotate("rect", xmin = a$sr_band[1], xmax = a$sr_band[2], ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.45) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = a$meta[1], ymax = a$meta[2],
           fill = "grey80", alpha = 0.45) +
  geom_vline(xintercept = a$sr_main, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_line(data = a$curves, aes(sigma_r, share, color = model), linewidth = 0.8) +
  geom_pointrange(data = a$anchors, aes(sigma_r, share, ymin = lo, ymax = hi, color = model),
                  fatten = 2.2, linewidth = 0.6, shape = 21, fill = "white", stroke = 0.9) +
  annotate("text", x = 0.255, y = mean(a$meta), label = "meta-analytic\nrange 4-7%",
           hjust = 0, size = 2.3, color = "grey30", lineheight = 0.85) +
  annotate("text", x = mean(a$sr_band), y = 0.205, label = "plausible population σ_r",
           size = 2.3, color = "grey30") +
  scale_color_manual(values = MODEL_PAL, name = NULL) +
  scale_y_continuous(labels = pct, limits = c(0, 0.22), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = seq(0.3, 0.8, 0.1)) +
  labs(x = expression(sigma[r] ~ "(population SD of log input rate)"),
       y = "implied input share of efficiency variance") +
  theme_fig(8) +
  theme(legend.position = c(0.97, 0.03), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

save_fig(p, "imputed_meta", width = PNAS_1COL, height = PNAS_1COL * 0.9)
