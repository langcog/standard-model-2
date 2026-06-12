## FINAL developmental-ladder analysis: 10 seeds x 18 rungs (0.5M..24M words),
## each cell = best-val competence of a GPT-2-small trained to convergence on
## that (nested) budget. Three outputs:
##   1. aggregate developmental curve + per-seed slope + between-seed sigma
##   2. per-word C&B 4-PL sigmoid on the DEVELOPMENT axis -> kappa-equivalent
##      (0.434/ParamScale on log10(words)), comparable to training-axis 0.73
##      and kids' ~12
##   3. between-seed sigma of the per-seed kappa -> the LM sigma_kappa analog
##
## Input:  fits/llm/ladder_bestval_finer.csv  (10 seeds x 18 budgets, complete)
## Output: figs/llm/ladder_development_final.png
##         fits/llm/ladder_kappa_summary.csv

suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(tidyr); library(patchwork)})

d <- read_csv("fits/llm/ladder_bestval_finer.csv", show_col_types = FALSE) |>
  mutate(surprisal = as.numeric(surprisal), lnw = log(words), l10w = log10(words))
NB <- n_distinct(d$words)   # full-ladder budget count (18)
cat(sprintf("loaded %d rows: %d seeds x %d rungs\n",
            nrow(d), n_distinct(d$seed), n_distinct(d$words)))

## ================= 1. aggregate curve + linear dev slope =================
agg <- d |> group_by(seed, words, lnw) |> summarise(mean_surp = mean(surprisal), .groups="drop")
slopes <- agg |> group_by(seed) |> summarise(lin_slope = coef(lm(mean_surp ~ lnw))[2], .groups="drop")
cat("\n== per-seed LINEAR dev slope (nats per e-fold) ==\n"); print(slopes |> mutate(across(-seed, ~round(.,4))))
cat(sprintf("mean = %.4f, SD = %.4f (CV %.1f%%)\n",
            mean(slopes$lin_slope), sd(slopes$lin_slope), 100*sd(slopes$lin_slope)/abs(mean(slopes$lin_slope))))

bs <- agg |> group_by(words) |> summarise(sd_seeds = sd(mean_surp), .groups="drop")
cat("\n== between-seed SD of competence by rung ==\n"); print(bs |> mutate(across(everything(), ~signif(.,3))))

## ================= 2. per-word 4-PL sigmoid on dev axis =================
four_pl_fit <- function(x, y) {
  ## C&B 4-parameter logistic: y = lo + (up-lo)/(1+exp((x-mid)/scale))
  st <- list(up = max(y), lo = min(y), mid = mean(x), scale = 0.5)
  fit <- tryCatch(
    nls(y ~ lo + (up - lo)/(1 + exp((x - mid)/scale)), start = st,
        control = nls.control(maxiter = 200, warnOnly = TRUE),
        lower = c(up = min(y)-5, lo = min(y)-5, mid = min(x)-3, scale = 1e-3),
        upper = c(up = max(y)+5, lo = max(y)+5, mid = max(x)+3, scale = 50),
        algorithm = "port"),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- coef(fit)
  tibble(up = co["up"], lo = co["lo"], mid = co["mid"], scale = co["scale"])
}

pw <- d |> group_by(seed, word) |> filter(n() == NB) |>
  group_modify(~{ r <- four_pl_fit(.x$l10w, .x$surprisal); if (is.null(r)) tibble() else r }) |>
  ungroup() |>
  mutate(kappa = 0.434 / scale, range = up - lo) |>
  filter(scale > 0.01, scale < 10, range > 1.0)   # C&B filters

kpseed <- pw |> group_by(seed) |>
  summarise(n_fit = n(), kappa_med = median(kappa),
            iqr_lo = quantile(kappa,.25), iqr_hi = quantile(kappa,.75), .groups="drop")
cat("\n== per-word DEV-axis kappa (0.434/ParamScale), per seed ==\n")
print(kpseed |> mutate(across(-c(seed,n_fit), ~round(.,3))))
cat(sprintf("\nkappa_pop_LM (dev axis, mean of per-seed medians) = %.3f\n", mean(kpseed$kappa_med)))
cat(sprintf("sigma_kappa_LM (SD of per-seed medians)            = %.4f\n", sd(kpseed$kappa_med)))
cat(sprintf("CV = %.1f%%   |   kids: kappa ~12, sigma_kappa ~3.5 (CV ~35%%)\n",
            100*sd(kpseed$kappa_med)/mean(kpseed$kappa_med)))
cat(sprintf("training-axis reference: GPT-2-CHILDES 24.5M full-training kappa ~0.74\n"))

write_csv(bind_rows(
  kpseed |> mutate(metric = "perword_kappa_dev_axis"),
  slopes |> rename(kappa_med = lin_slope) |> mutate(metric = "linear_dev_slope")),
  "fits/llm/ladder_kappa_summary.csv")

## ================= 3. figure =================
p1 <- ggplot(agg, aes(words, mean_surp, color = factor(seed))) +
  geom_line(alpha=.8) + geom_point(size=1.2) +
  scale_x_log10(breaks = c(.5,1,2,4,8,16,24)*1e6, labels=c("0.5M","1M","2M","4M","8M","16M","24M")) +
  labs(x = "training budget (words, log)", y = "held-out CDI surprisal (nats)", color = "seed") +
  theme_minimal(base_size = 10.5) + theme(legend.position = c(.85,.75))

p2 <- ggplot(pw, aes(kappa, color = factor(seed), fill = factor(seed))) +
  geom_density(alpha=.08) +
  geom_vline(xintercept = 0.74, linetype = "dotted") +
  annotate("text", x=.78, y=Inf, label="training-axis LM kappa 0.74", hjust=0, vjust=1.4, size=2.7) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(x = "per-word developmental kappa (0.434/scale)", y = "density",
       color = "seed", fill = "seed") +
  theme_minimal(base_size = 10.5) + theme(legend.position = "none")

g <- (p1 | p2) + plot_annotation(
  title = "LM development across 10 individuals (seeds): near-identical, decelerating, shallow",
  subtitle = sprintf(paste0("Left: competence vs distinct input, 18 budgets to convergence; between-seed slope SD %.3f (CV %.0f%%). ",
                            "Right: per-word dev-axis kappa, median %.2f — vs children kappa ~12, sigma_kappa ~3.5 (CV ~35%%)."),
                     sd(slopes$lin_slope), 100*sd(slopes$lin_slope)/abs(mean(slopes$lin_slope)),
                     mean(kpseed$kappa_med)),
  theme = theme(plot.title = element_text(face="bold", size=12),
                plot.subtitle = element_text(size=8.6, color="grey30")))
ggsave("figs/llm/ladder_development_final.png", g, width = 10.5, height = 4.6, dpi = 150)
cat("\nwrote figs/llm/ladder_development_final.png\n")
