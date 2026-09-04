## io-proc figure (two column): the direct input + processing study.
##   (A) schematic: how a +/-2 SD change in input vs processing moves a
##       lower- and a higher-ability child's trajectory, drawn from the JOINT
##       fit's own estimates (io_fig$schem), in latent-ability and words-produced space.
##   (B) coefficients: each channel x trait effect shown twice -- fit ALONE
##       (io-only / proc-only arm, open) and in the JOINT model (filled). All
##       three arms fit the same english_count bundle, so the comparison is
##       exact. Per-SD theta units; the acceleration row is rescaled by
##       log(30/21) to level-equivalent units at 30 months. 90% CIs.
## The imputed-population sweep that used to be panel C lives in
## fig_imputed_meta.R now.
##
## Reads: paper/cache/fig_io_imputed_proc.rds (schem + panelB + panelB_arms,
##        built by paper/build_fig_io_cache.R).
## Run:   Rscript paper/figures/fig_ioproc.R
source(here::here("paper", "figures", "_common.R"))

io_fig <- readRDS(file.path(CACHE, "fig_io_imputed_proc.rds"))

## ---- (A) schematic: per-SD input vs processing perturbations ----
sc <- io_fig$schem
a0 <- sc$a0; logH <- sc$log_H; mu_r <- sc$mu_r; kpop <- sc$kpop; sig_xi <- sc$sigma_xi
SEP <- 1.35; SD <- 2
d_in_xi <- sc$d_in_xi * SD; d_in_k <- sc$d_in_k * SD     # input: level AND slope
d_pr_xi <- sc$d_pr_xi * SD; d_pr_k <- sc$d_pr_k * SD     # processing (faster): level (+ small slope)
mu_d <- 15; scl <- 2.5; Nitems <- 680
A_of  <- function(xi, k, age) xi + logH + k * log(age / a0)
vocab <- function(A) Nitems * plogis((A - mu_d) / scl)
ages  <- seq(12, 36, length.out = 120)
QL <- c("latent ability (θ)", "words produced")
kids  <- tibble(kid = names(KID_PAL),
                xi0 = c(mu_r - SEP * sig_xi, mu_r + SEP * sig_xi))
conds <- tibble(
  cond = c("base", "input +2SD", "input -2SD", "processing +2SD", "processing -2SD"),
  chan = c("base", "input", "input", "processing", "processing"),
  dxi  = c(0, d_in_xi, -d_in_xi, d_pr_xi, -d_pr_xi),
  dk   = c(0, d_in_k,  -d_in_k,  d_pr_k,  -d_pr_k))
dfS <- tidyr::expand_grid(kids, conds) |> rowwise() |>
  mutate(d = list(tibble(age = ages,
                         th = A_of(xi0 + dxi, kpop + dk, ages),
                         wd = vocab(A_of(xi0 + dxi, kpop + dk, ages))))) |>
  ungroup() |> tidyr::unnest(d) |>
  tidyr::pivot_longer(c(th, wd), names_to = "quantity", values_to = "value") |>
  mutate(chan = factor(chan, levels = c("base", "input", "processing")),
         quantity = factor(ifelse(quantity == "th", QL[1], QL[2]), levels = QL))
pA <- ggplot(dfS, aes(age, value, color = kid, linetype = chan,
                      group = interaction(kid, cond, quantity))) +
  geom_line(linewidth = 0.55) +
  facet_wrap(~ quantity, scales = "free_y") +
  scale_linetype_manual(values = c(base = "solid", input = "dashed", processing = "dotted"),
                        labels = c("base", "±2 SD input", "±2 SD processing"), name = NULL) +
  scale_color_manual(values = KID_PAL, name = NULL) +
  labs(x = "age (months)", y = NULL,
       title = "A. A ±2 SD change in input vs processing") +
  theme_fig(8) +
  theme(legend.position = "right")

## ---- (B) coefficients: separate arm (open) vs joint (filled) ----
KACC <- log(30 / 21)
pb <- io_fig$panelB_arms |>
  mutate(acc = channel == "Acceleration",
         e = ifelse(acc, med * KACC, med), l = ifelse(acc, lo * KACC, lo), h = ifelse(acc, hi * KACC, hi),
         effect = factor(paste(factor, "→", ifelse(acc, "acceleration\n(by 30 mo)", "efficiency")),
                         levels = rev(c("Input → efficiency", "Processing → efficiency",
                                        "Input → acceleration\n(by 30 mo)",
                                        "Processing → acceleration\n(by 30 mo)"))),
         model = factor(model, levels = c("Separate model", "Joint model")),
         factor = factor(factor, levels = c("Input", "Processing")))
pB <- ggplot(pb, aes(x = e, y = effect, color = factor, shape = model)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_errorbar(aes(xmin = l, xmax = h), position = position_dodge(width = 0.55),
                width = 0, linewidth = 0.6, orientation = "y") +
  geom_point(position = position_dodge(width = 0.55), size = 2.4, fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c("Separate model" = 21, "Joint model" = 16), name = NULL) +
  scale_color_manual(values = FACTOR_PAL, guide = "none") +
  labs(x = "effect of +1 SD (θ units, log-odds), 90% CI", y = NULL,
       title = "B. Coefficients: each channel alone vs jointly") +
  theme_fig(8) +
  theme(legend.position = "top", legend.justification = "left",
        panel.grid.major.y = element_blank())

p <- pA / pB + plot_layout(heights = c(1, 1))
save_fig(p, "ioproc", width = PNAS_2COL, height = PNAS_2COL * 0.72)
