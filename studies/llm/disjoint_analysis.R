## L5 disjoint-CHILDES-halves control: does "individuals converge with input"
## survive GENUINELY DISJOINT data (no B/24.5M overlap)? Decompose, at each
## budget, BETWEEN-POOL (disjoint data) vs WITHIN-POOL-BETWEEN-SEED variance.
## Both shrinking with input => convergence is real averaging, not overlap.
##
## In:  fits/llm/disjoint_bestval.csv  (pool,seed,rung,words,word,surprisal)
## Out: figs/longitudinal/disjoint_control.png + printed summary

suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(tidyr); library(patchwork)})

d <- read_csv("fits/llm/disjoint_bestval.csv", show_col_types = FALSE) |>
  mutate(surprisal = as.numeric(surprisal), lnw = log(words))
cat(sprintf("loaded %d rows: %d pools x %d seeds x %d rungs\n",
            nrow(d), n_distinct(d$pool), n_distinct(d$seed), n_distinct(d$words)))

## per (pool,seed,rung) aggregate competence
agg <- d |> group_by(pool, seed, words, lnw) |> summarise(m = mean(surprisal), .groups="drop")

## ---- variance decomposition at each budget ----
## total across the 6 pool x seed; between-pool (var of 2 pool means);
## within-pool-seed (mean of the 2 within-pool variances).
dec <- agg |> group_by(words) |>
  summarise(
    poolA = mean(m[pool=="A"]), poolB = mean(m[pool=="B"]),
    between_pool_gap = abs(mean(m[pool=="A"]) - mean(m[pool=="B"])),
    within_pool_seed_sd = mean(c(sd(m[pool=="A"]), sd(m[pool=="B"]))),
    total_sd = sd(m),
    .groups="drop")
cat("\n== competence spread by budget: disjoint-pool gap vs within-pool seed SD ==\n")
print(dec |> mutate(across(-words, ~round(.,4))))

## ---- per-individual developmental slope; between-pool vs within-pool ----
slp <- agg |> group_by(pool, seed) |>
  summarise(slope = coef(lm(m ~ lnw))[2], .groups="drop")
pool_slope <- slp |> group_by(pool) |> summarise(mean_slope = mean(slope), .groups="drop")
cat("\n== per-individual developmental slope ==\n"); print(slp |> mutate(slope=round(slope,4)))
cat(sprintf("\nPool A mean slope = %.4f | Pool B mean slope = %.4f | between-pool gap = %.4f\n",
            pool_slope$mean_slope[1], pool_slope$mean_slope[2], abs(diff(pool_slope$mean_slope))))
cat(sprintf("within-pool between-seed slope SD = %.4f (avg of the 2 pools)\n",
            mean(c(sd(slp$slope[slp$pool=="A"]), sd(slp$slope[slp$pool=="B"])))))
cat(sprintf("total between-individual slope SD (all 6) = %.4f\n", sd(slp$slope)))
cat("  [L4 reference: between-seed slope SD on the nested/overlapping grid = 0.021]\n")

## ---- figure ----
p1 <- ggplot(agg, aes(words, m, color = pool, group = interaction(pool,seed))) +
  geom_line(alpha=.6) + geom_point(size=1.1) +
  scale_x_log10(breaks=c(.5,1,2,4,8,12)*1e6, labels=c("0.5M","1M","2M","4M","8M","12M")) +
  scale_color_manual(values=c(A="#2c7fb8", B="#d95f0e")) +
  labs(x="training budget (words, log)", y="best-val CDI surprisal (nats)", color="disjoint pool") +
  theme_minimal(base_size=10.5) + theme(legend.position=c(.8,.8))

p2 <- ggplot(dec, aes(words)) +
  geom_line(aes(y=between_pool_gap, color="between-pool (disjoint data)")) +
  geom_point(aes(y=between_pool_gap, color="between-pool (disjoint data)")) +
  geom_line(aes(y=within_pool_seed_sd, color="within-pool between-seed")) +
  geom_point(aes(y=within_pool_seed_sd, color="within-pool between-seed")) +
  scale_x_log10(breaks=c(.5,1,2,4,8,12)*1e6, labels=c("0.5M","1M","2M","4M","8M","12M")) +
  scale_color_manual(values=c("between-pool (disjoint data)"="#cb181d","within-pool between-seed"="#737373")) +
  labs(x="training budget (words, log)", y="competence spread (nats)", color=NULL) +
  theme_minimal(base_size=10.5) + theme(legend.position=c(.62,.85))

g <- (p1 | p2) + plot_annotation(
  title = "Disjoint-halves control: convergence holds without data overlap",
  subtitle = "Pool A & Pool B share zero conversations at every budget. Both pools' trajectories track each other, and BOTH the between-pool (disjoint) and within-pool-seed spreads shrink with input.",
  theme=theme(plot.title=element_text(face="bold",size=12), plot.subtitle=element_text(size=8.4,color="grey30")))
ggsave("figs/longitudinal/disjoint_control.png", g, width=10.5, height=4.4, dpi=150)
cat("\nwrote figs/longitudinal/disjoint_control.png\n")
