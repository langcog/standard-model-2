## Trajectory demo = panel B of the graphical-model figure: what the model's
## two child-level traits do to a trajectory. One reference child (0.5 SD
## below the average efficiency), then the same child with +1 SD efficiency
## (xi: a parallel shift of theta at every age) and with +1 SD acceleration
## (kappa: pivots on the anchor age a0 -- identical at a0, below the reference
## before it, above after), in latent-ability and words-produced space. The
## population parameters (mu_xi, kappa, sigma_xi, sigma_kappa, log H, a0) are
## the joint io-proc fit's (io_fig$schem). Direct-labelled in the trait hues
## (TRAIT_PAL); the reference child neutral.
## Writes out/model_trajectories.pdf, a COMPONENT that fig_graphical_model.tex
## places beside the TikZ diagram (no panel letter here: the .tex adds A / B).
## Sized to the diagram's natural height (8.6 cm) at single-column width.
##
## Reads: paper/cache/fig_io_imputed_proc.rds (schem; built by
##        paper/build_fig_io_cache.R).
## Run:   Rscript paper/figures/fig_model_trajectories.R   (then the .tex)
source(here::here("paper", "figures", "_common.R"))

sc <- readRDS(file.path(CACHE, "fig_io_imputed_proc.rds"))$schem
A_of  <- function(xi, k, age) xi + sc$log_H + k * log(age / sc$a0)
vocab <- function(A) 680 * plogis((A - 15) / 2.5)
ages  <- seq(12, 36, length.out = 120)
xi0   <- sc$mu_r - 0.5 * sc$sigma_xi                   # reference child: slightly below average
conds <- tibble(
  cond = c("reference child", "+1 SD efficiency (ξ)", "+1 SD acceleration (κ)"),
  trait = c("base", "Efficiency", "Acceleration"),
  xi = xi0 + c(0, sc$sigma_xi, 0),
  k  = sc$kpop + c(0, 0, sc$sigma_kappa))
QL <- c("latent ability (θ)", "words produced")
dfS <- tidyr::expand_grid(conds, age = ages) |>
  mutate(th = A_of(xi, k, age), wd = vocab(th)) |>
  tidyr::pivot_longer(c(th, wd), names_to = "quantity", values_to = "value") |>
  mutate(quantity = factor(ifelse(quantity == "th", QL[1], QL[2]), levels = QL),
         trait = factor(trait, levels = c("base", names(TRAIT_PAL))))

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
a0_lab <- tibble(quantity = factor(QL[1], levels = QL), x = sc$a0, y = rng$lo[1],
                 label = sprintf("a[0] == %g~mo", sc$a0))

p <- ggplot(dfS, aes(age, value, colour = trait, group = cond)) +
  geom_vline(xintercept = sc$a0, linetype = "13", colour = "grey55", linewidth = 0.35) +
  geom_text(data = a0_lab, aes(x, y, label = label), inherit.aes = FALSE, parse = TRUE,
            hjust = -0.1, vjust = 0, size = TXT, colour = "grey40") +
  geom_line(linewidth = 0.6) +
  geom_segment(data = ends, aes(x = 36.2, xend = 37.6, y = value, yend = ylab),
               inherit.aes = FALSE, colour = "grey60", linewidth = 0.25) +
  geom_text(data = ends, aes(x = 38, y = ylab, label = cond), inherit.aes = FALSE,
            hjust = 0, size = TXT, colour = "grey20") +
  facet_wrap(~ quantity, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(base = VOCAB_COL, TRAIT_PAL), guide = "none") +
  scale_x_continuous(breaks = c(12, 18, 24, 30, 36), expand = expansion(mult = c(0.03, 0.55))) +
  labs(x = "age (months)", y = NULL, title = "One child, +1 SD of efficiency vs acceleration") +
  theme_fig()

save_fig(p, "model_trajectories", width = PNAS_1COL, height = 8.6 / 2.54)
