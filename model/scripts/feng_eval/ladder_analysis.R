## Developmental-ladder analysis: competence (best-val held-out CDI surprisal)
## vs distinct-input budget, across seeds. Each seed is an "individual"; the
## spread across seeds of the developmental slope is the LM analog of children's
## between-child sigma_kappa (~3.5).
##
## Input:  outputs/feng_eval/ladder_bestval_partial.csv  (seed,rung,words,word,surprisal)
## Output: outputs/figs/longitudinal/ladder_development.png + printed summary.

suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(tidyr)})

d <- read_csv("outputs/feng_eval/ladder_bestval_partial.csv", show_col_types = FALSE) |>
  mutate(surprisal = as.numeric(surprisal), lnw = log(words))
cat(sprintf("loaded %d rows: %d seeds x %d rungs x ~%d words\n",
            nrow(d), n_distinct(d$seed), n_distinct(d$words),
            round(nrow(d)/n_distinct(d$seed)/n_distinct(d$words))))

## ---- aggregate developmental curve: mean competence over words, per (seed,rung) ----
agg <- d |> group_by(seed, words, lnw) |> summarise(mean_surp = mean(surprisal), .groups = "drop")

## per-seed developmental slope kappa_seed = d(surprisal)/d(ln words)  (negative = learning)
slopes <- agg |> group_by(seed) |>
  summarise(kappa = coef(lm(mean_surp ~ lnw))[2], .groups = "drop")
cat("\n== per-seed developmental slope (nats surprisal per e-fold of input) ==\n")
print(slopes |> mutate(kappa = round(kappa, 4)))
cat(sprintf("\nkappa_pop_LM (mean over seeds) = %.4f\n", mean(slopes$kappa)))
cat(sprintf("sigma_kappa_LM (SD across seeds) = %.4f   [CV = %.1f%%]\n",
            sd(slopes$kappa), 100*sd(slopes$kappa)/abs(mean(slopes$kappa))))

## between-seed SD of mean competence at each rung (do seeds converge?)
cat("\n== between-seed spread of competence at each rung ==\n")
bs <- agg |> group_by(words) |>
  summarise(mean = mean(mean_surp), sd_across_seeds = sd(mean_surp), .groups = "drop")
print(bs |> mutate(across(c(mean, sd_across_seeds), ~round(., 4))))

## ---- per-word developmental slopes (each seed): distribution + between-seed stability ----
pw <- d |> group_by(seed, word) |>
  summarise(slope = coef(lm(surprisal ~ lnw))[2], .groups = "drop")
pw_seed <- pw |> group_by(seed) |>
  summarise(median_slope = median(slope), iqr_lo = quantile(slope,.25), iqr_hi = quantile(slope,.75),
            .groups = "drop")
cat("\n== per-word slope distribution, per seed (median [IQR]) ==\n")
print(pw_seed |> mutate(across(-seed, ~round(., 3))))
cat(sprintf("\nbetween-seed SD of per-seed MEDIAN per-word slope = %.4f\n",
            sd(pw_seed$median_slope)))

## ---- figure: developmental curves by seed ----
seed_med <- mean(c(0.72,0.74,0.74))
p <- ggplot(agg, aes(words, mean_surp, color = factor(seed))) +
  geom_line(alpha = .8) + geom_point(size = 1.3) +
  scale_x_log10(breaks = c(.5,1,2,4,8,16)*1e6, labels = c("0.5M","1M","2M","4M","8M","16M")) +
  labs(title = "LM development: competence vs distinct input, 5 seeds (10 rungs; 24M pending)",
       subtitle = sprintf("Between-seed sigma of developmental slope = %.3f  (kids sigma_kappa ~ 3.5, on a steeper, rising axis)",
                          sd(slopes$kappa)),
       x = "training budget (words, log scale)",
       y = "best-val held-out CDI surprisal (nats; lower = more competent)",
       color = "seed") +
  theme_minimal(base_size = 11)
ggsave("outputs/figs/longitudinal/ladder_development.png", p, width = 9, height = 5, dpi = 150)
cat("\nwrote outputs/figs/longitudinal/ladder_development.png\n")
