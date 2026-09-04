## io-proc figure (two column): the direct input + processing study.
##   (A) schematic: how a +/-2 SD change in input vs processing moves an
##       average child's trajectory, drawn from the JOINT fit's own median
##       estimates (io_fig$schem), in latent-ability and words-produced space.
##       Direct-labelled; input green, processing purple (CHANNEL_PAL).
##   (B) coefficients: each channel x trait effect shown twice -- fit ALONE
##       (io-only / proc-only arm, open) and in the JOINT model (filled). All
##       three arms fit the same english_count bundle, so the comparison is
##       exact. Per-SD theta units; the acceleration facet is rescaled by
##       log(30/21) to level-equivalent units at 30 months. 90% CIs.
## The imputed-population sweep that used to be panel C lives in
## fig_imputed_meta.R now.
##
## Reads: paper/cache/fig_io_imputed_proc.rds (schem + panelB_arms,
##        built by paper/build_fig_io_cache.R).
## Run:   Rscript paper/figures/fig_ioproc.R
source(here::here("paper", "figures", "_common.R"))

io_fig <- readRDS(file.path(CACHE, "fig_io_imputed_proc.rds"))

## ---- (A) schematic: per-SD input vs processing perturbations ----
sc <- io_fig$schem
SD <- 2
A_of  <- function(xi, k, age) xi + sc$log_H + k * log(age / sc$a0)
vocab <- function(A) 680 * plogis((A - 15) / 2.5)
ages  <- seq(12, 36, length.out = 120)
conds <- tibble(
  cond = c("average child", "+2 SD input", "−2 SD input", "+2 SD processing", "−2 SD processing"),
  chan = c("base", "Input", "Input", "Processing", "Processing"),
  sgn  = c(0, 1, -1, 1, -1),
  dxi  = sgn * SD * c(0, sc$d_in_xi, sc$d_in_xi, sc$d_pr_xi, sc$d_pr_xi),
  dk   = sgn * SD * c(0, sc$d_in_k,  sc$d_in_k,  sc$d_pr_k,  sc$d_pr_k))
QL <- c("latent ability (θ)", "words produced")
dfS <- tidyr::expand_grid(conds, age = ages) |>
  mutate(th = A_of(sc$mu_r + dxi, sc$kpop + dk, age), wd = vocab(th)) |>
  tidyr::pivot_longer(c(th, wd), names_to = "quantity", values_to = "value") |>
  mutate(quantity = factor(ifelse(quantity == "th", QL[1], QL[2]), levels = QL),
         chan = factor(chan, levels = c("base", "Input", "Processing")),
         lty  = ifelse(sgn < 0, "minus", "solid"))
## direct labels at the right end of each line, pushed apart vertically so
## none collide (deterministic; the labels keep the lines' vertical order)
spread_labels <- function(y, gap, lo, hi) {
  o <- order(y); ys <- y[o]
  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
  ys <- ys - (mean(ys) - mean(y))            # block centred on the line ends ...
  ys <- ys - max(0, max(ys) - hi)            # ... but kept inside the panel
  ys <- ys + max(0, lo - min(ys))
  y[o] <- ys; y
}
rng  <- dfS |> group_by(quantity) |> summarise(lo = min(value), hi = max(value), .groups = "drop")
ends <- dfS |> filter(age == max(age)) |> left_join(rng, by = "quantity") |> group_by(quantity) |>
  mutate(ylab = spread_labels(value, gap = 0.13 * (hi[1] - lo[1]), lo[1], hi[1])) |>
  ungroup()
pA <- ggplot(dfS, aes(age, value, colour = chan, linetype = lty, group = cond)) +
  geom_line(linewidth = 0.55) +
  geom_segment(data = ends, aes(x = 36.2, xend = 37.6, y = value, yend = ylab),
               inherit.aes = FALSE, colour = "grey60", linewidth = 0.25) +
  geom_text(data = ends, aes(x = 38, y = ylab, label = cond), inherit.aes = FALSE,
            hjust = 0, size = TXT, colour = "grey20") +
  facet_wrap(~ quantity, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(base = VOCAB_COL, CHANNEL_PAL), guide = "none") +
  scale_linetype_manual(values = c(solid = "solid", minus = "22"), guide = "none") +
  scale_x_continuous(breaks = c(12, 18, 24, 30, 36), expand = expansion(mult = c(0.03, 0.4))) +
  labs(x = "age (months)", y = NULL, title = "±2 SD of input vs processing, average child") +
  theme_fig()

## ---- (B) coefficients: separate arm (open) vs joint (filled) ----
KACC <- log(30 / 21)
OUT  <- c("Effect on efficiency", "Effect on acceleration (by 30 mo)")
pb <- io_fig$panelB_arms |>
  mutate(acc = channel == "Acceleration",
         e = ifelse(acc, med * KACC, med), l = ifelse(acc, lo * KACC, lo), h = ifelse(acc, hi * KACC, hi),
         outcome = factor(ifelse(acc, OUT[2], OUT[1]), levels = OUT),
         model = factor(model, levels = c("Separate model", "Joint model")),
         factor = factor(factor, levels = rev(names(CHANNEL_PAL))))
pB <- ggplot(pb, aes(e, factor, colour = factor, shape = model)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_errorbar(aes(xmin = l, xmax = h), position = position_dodge(width = 0.6),
                width = 0, linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 1.9, fill = "white", stroke = 0.8) +
  ggh4x::facet_wrap2(~ outcome, ncol = 1, strip = trait_strips()) +
  scale_shape_manual(values = c("Separate model" = 21, "Joint model" = 16), name = NULL) +
  guides(shape = guide_legend(override.aes = list(colour = "grey20", fill = "white"))) +
  scale_colour_manual(values = CHANNEL_PAL, guide = "none") +
  labs(x = "effect of +1 SD (θ units, log-odds), 90% CI", y = NULL,
       title = "Each channel alone vs jointly") +
  theme_fig() +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0), panel.grid.major.y = element_blank())

p <- (pA | pB) + plot_layout(widths = c(1.45, 1))
save_fig(tag_panels(p), "ioproc", width = PNAS_2COL, height = 3.4)
