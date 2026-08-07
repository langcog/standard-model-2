#!/usr/bin/env Rscript
# CANDIDATE FIGURE (standalone; not wired into the paper): the CHILDES LMs show DIMINISHING
# returns to experience -- their effective kappa DECLINES as they accumulate data, the mirror
# image of children (kappa ~10, ~flat). Artifact-free: measured from improvement-per-decade
# directly (NOT per-word sigmoid midpoints, which are window-confounded; see ONSET_KAPPA_NOTE.md).
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork)})
ROOT <- "/Users/mcfrank/Projects/standard_model_2"; OUT <- file.path(ROOT,"studies/llm")

lad <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"), show_col_types=FALSE) |>
  mutate(lx = log10(words))
# mean competence (surprisal) per seed x budget, pooled over words
comp <- lad |> group_by(seed, words, lx) |> summarise(surp = mean(surprisal), .groups="drop")
budg <- comp |> group_by(words, lx) |> summarise(surp = mean(surp), .groups="drop") |> arrange(lx)

# effective kappa = 0.434 * (nats improvement per log10-decade), local between consecutive rungs.
# (logit-competence ~ -surprisal for these low-prob CDI words; 0.434 = 1/ln10 -> natural-log kappa)
loc <- function(d) d |> arrange(lx) |> mutate(
  k_eff = 0.434 * -(surp - lag(surp)) / (lx - lag(lx)),
  lx_mid = (lx + lag(lx))/2, words_mid = 10^lx_mid)
per_seed <- comp |> group_by(seed) |> group_modify(~loc(.x)) |> ungroup() |> filter(is.finite(k_eff))
band <- per_seed |> group_by(words_mid, lx_mid) |>
  summarise(k = mean(k_eff), lo = quantile(k_eff,.1), hi = quantile(k_eff,.9), .groups="drop")

# smooth quadratic on pooled competence -> analytic derivative -> kappa curve
qfit <- lm(surp ~ lx + I(lx^2), budg); b <- coef(qfit)
kcurve <- tibble(lx = seq(min(budg$lx), max(budg$lx), length=100)) |>
  mutate(words = 10^lx, k_eff = 0.434 * -(b[2] + 2*b[3]*lx))
b2 <- summary(qfit)$coef[3,]
k_early <- with(subset(band, words_mid==min(words_mid)), k); k_late <- with(subset(band, words_mid==max(words_mid)), k)
cat(sprintf("effective kappa: early=%.2f  late=%.2f   quad b2=%.3f (p=%.2g)\n", k_early, k_late, b2[1], b2[4]))

pal_lm <- "#2c7fb8"
# Panel A: competence flattens (raw phenomenon)
pA <- ggplot(budg, aes(words, surp)) +
  geom_line(data=tibble(words=10^kcurve$lx, surp=predict(qfit, kcurve)), colour=pal_lm, linewidth=.7) +
  geom_point(size=2, colour=pal_lm) +
  scale_x_log10(breaks=c(5e5,1e6,3e6,1e7,2.4e7), labels=c("0.5M","1M","3M","10M","24M")) +
  labs(title="(A) Competence gains flatten with data",
       x="training data (distinct words)", y="mean surprisal (lower = better)") +
  theme_minimal(base_size=11)

# Panel B: effective kappa declines below unit accumulator; children ~10 flat
pB <- ggplot(band, aes(words_mid, k)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill=pal_lm, alpha=.15) +
  geom_line(data=kcurve, aes(words, k_eff), colour=pal_lm, linewidth=.7) +
  geom_point(size=2, colour=pal_lm) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey55") +
  annotate("text", x=6e5, y=1.06, hjust=0, size=3, colour="grey40",
           label="unit accumulator (kappa = 1)") +
  annotate("text", x=6e5, y=0.35, hjust=0, size=3.1, colour=pal_lm, fontface="bold",
           label=sprintf("model decelerates:\nkappa %.2f -> %.2f", k_early, k_late)) +
  annotate("segment", x=8e6, xend=8e6, y=1.9, yend=2.15,
           arrow=arrow(length=unit(.15,"cm")), colour="#c41e37") +
  annotate("text", x=8e6, y=1.95, vjust=0, size=3.1, colour="#c41e37",
           label="Children: kappa ~ 10,\n~flat across development") +
  scale_x_log10(breaks=c(5e5,1e6,3e6,1e7,2.4e7), labels=c("0.5M","1M","3M","10M","24M")) +
  coord_cartesian(ylim=c(0.2, 2.2)) +
  labs(title="(B) Returns to experience DECLINE (diminishing)",
       x="training data (distinct words)", y="effective kappa (returns to experience)") +
  theme_minimal(base_size=11)

g <- pA + pB + plot_annotation(
  title="CHILDES language models show diminishing returns to experience",
  subtitle="one learner over its own accumulating data (10 seeds); the mirror image of children's accelerating returns",
  theme=theme(plot.title=element_text(face="bold", size=12)))
ggsave(file.path(OUT,"fig_kappa_decline.png"), g, width=10, height=4.3, dpi=140)
cat("wrote fig_kappa_decline.png\n")
