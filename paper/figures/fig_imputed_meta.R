## Imputed population input share vs the meta-analytic range (single column).
## Factored out of the io-proc figure so the population-level imputation is
## not confused with the direct io-proc study: analytic curves
## share = sigma_r^2 / sigma_xi^2 for the EN and NO D fits, the six real
## refit anchors at sigma_r = {0.35, 0.44, 0.58} with 95% CIs, the plausible
## sigma_r band [0.32, 0.58] and the meta-analytic 4-7% band. Both curves are
## in the input hue (CHANNEL_PAL) -- the quantity is input's share -- and are
## direct-labelled by language (solid English, dashed Norwegian). No panel
## letter: single panel.
##
## Reads: paper/cache/fig_io_imputed_proc.rds (panelA; built by
##        paper/build_fig_io_cache.R from the _2k anchor fits).
## Run:   Rscript paper/figures/fig_imputed_meta.R
source(here::here("paper", "figures", "_common.R"))

io_fig <- readRDS(file.path(CACHE, "fig_io_imputed_proc.rds"))
a <- io_fig$panelA
pct <- scales::percent_format(accuracy = 1)
LANG <- c(`English (D)` = "English", `Norwegian (D)` = "Norwegian")
curves  <- a$curves  |> mutate(model = LANG[model])
anchors <- a$anchors |> mutate(model = LANG[model])
ends    <- curves |> group_by(model) |> slice_max(sigma_r, n = 1) |> ungroup()
INP <- CHANNEL_PAL[["Input"]]

p <- ggplot() +
  annotate("rect", xmin = a$sr_band[1], xmax = a$sr_band[2], ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = a$meta[1], ymax = a$meta[2],
           fill = "grey85", alpha = 0.5) +
  geom_vline(xintercept = a$sr_main, linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_line(data = curves, aes(sigma_r, share, linetype = model), colour = INP, linewidth = 0.7) +
  geom_errorbar(data = anchors, aes(sigma_r, ymin = lo, ymax = hi), colour = INP,
                width = 0, linewidth = 0.5) +
  geom_point(data = anchors, aes(sigma_r, share), colour = INP, shape = 21, fill = "white",
             size = 1.8, stroke = 0.8) +
  geom_text(data = ends, aes(sigma_r + 0.012, share, label = model), hjust = 0, size = TXT,
            colour = "grey20") +
  annotate("text", x = 0.255, y = mean(a$meta), label = "meta-analytic\nrange 4–7%",
           hjust = 0, size = TXT, color = "grey30", lineheight = 0.9) +
  annotate("text", x = mean(a$sr_band), y = 0.207, label = "plausible~population~sigma[r]",
           parse = TRUE, size = TXT, color = "grey30") +
  scale_linetype_manual(values = c(English = "solid", Norwegian = "22"), guide = "none") +
  scale_y_continuous(labels = pct, limits = c(0, 0.22), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = seq(0.3, 0.8, 0.1), expand = expansion(mult = c(0.02, 0.27))) +
  labs(x = expression(sigma[r] ~ "(population SD of log input rate)"),
       y = "implied input share of efficiency variance") +
  theme_fig()

save_fig(p, "imputed_meta", width = PNAS_1COL, height = PNAS_1COL * 0.9)
