## Model-over-development returns from the CHILDES ladder (converged models at
## increasing distinct-input budgets). The level-matched analogue of children's
## kappa_c: does ONE learner's competence ACCELERATE (increasing returns ->
## learning-to-learn) or DECELERATE (diminishing returns -> scaling law) as it
## accumulates DATA?
##
## Reproduce from saved runs (run from repo root):
##   Rscript studies/llm/autoresearch_pilot/ladder_returns.R
## Input : fits/llm/ladder_bestval_finer.csv   (10 seeds x 18 budgets, per-word
##         converged best-val surprisal; built by studies/llm/extract_ladder_bestval.py)
## Output: studies/llm/autoresearch_pilot/ladder_returns.png + console stats
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(here)})

d <- read_csv(here("fits/llm/ladder_bestval_finer.csv"), show_col_types = FALSE) |>
  mutate(words = as.numeric(words), surprisal = as.numeric(surprisal), lb = log10(words))
nb <- max(count(d, seed, word)$n)
d <- d |> group_by(seed, word) |> filter(n() == nb) |> ungroup()
cat("seeds:", n_distinct(d$seed), " budgets:", n_distinct(d$words), " words:", n_distinct(d$word), "\n")

# aggregate competence + marginal improvement per decade of data
agg <- d |> group_by(words, lb) |> summarise(mean_surp = mean(surprisal), .groups = "drop") |>
  arrange(lb) |> mutate(nats_per_decade = (mean_surp - lag(mean_surp)) / (lb - lag(lb)),
                        kappa_approx = -nats_per_decade / log(10))  # logit ~= -surprisal for low-P words
cat("\n== competence by budget (improvement per decade; effective kappa ~= -slope/ln10) ==\n")
print(as.data.frame(agg |> transmute(words, mean_surp = round(mean_surp, 3),
                                      nats_per_decade = round(nats_per_decade, 3),
                                      kappa_approx = round(kappa_approx, 3))))

qfit <- lm(mean_surp ~ lb + I(lb^2), agg); b2 <- coef(qfit)["I(lb^2)"]
cat(sprintf("\naggregate quadratic coef (lb^2) = %.3f -> %s\n", b2,
  ifelse(b2 > 0, "CONVEX: drop DECELERATES = DIMINISHING returns (kappa falls with data)",
                 "CONCAVE: drop ACCELERATES = INCREASING returns")))
pw <- d |> group_by(seed, word) |>
  summarise(b2 = tryCatch(coef(lm(surprisal ~ lb + I(lb^2)))[3], error = function(e) NA), .groups = "drop")
cat(sprintf("per-word quadratic: frac diminishing (b2>0) = %.2f, median b2 = %.3f\n",
            mean(pw$b2 > 0, na.rm = TRUE), median(pw$b2, na.rm = TRUE)))
cat(sprintf("effective kappa: early(0.5-1M) ~= %.2f, late(16-24M) ~= %.2f  | children kappa_c ~= 10 (flat)\n",
            -mean(agg$nats_per_decade[2:3]) / log(10), -mean(agg$nats_per_decade[17:18]) / log(10)))

g <- ggplot(agg, aes(words, mean_surp)) + geom_line(colour = "#2c7fb8") + geom_point(size = 1.8) +
  scale_x_log10(breaks = c(.5,1,2,4,8,16,24)*1e6, labels = c("0.5M","1M","2M","4M","8M","16M","24M")) +
  labs(title = "CHILDES ladder: converged competence vs data budget (diminishing returns)",
       subtitle = "per-decade drop shrinks -> effective kappa falls ~1.3->0.5; opposite of child acceleration (kappa_c~=10)",
       x = "distinct-input budget (words, log)", y = "mean held-out CDI surprisal (nats)") +
  theme_minimal(base_size = 11)
ggsave(here("studies/llm/autoresearch_pilot/ladder_returns.png"), g, width = 7.5, height = 4.4, dpi = 150)
cat("\nwrote studies/llm/autoresearch_pilot/ladder_returns.png\n")
