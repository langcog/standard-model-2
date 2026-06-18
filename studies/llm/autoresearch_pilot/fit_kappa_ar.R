## Fit per-CDI-word acquisition sigmoids from the autoresearch probe trajectory.
## Input: word_surprisal.csv (step,epoch,word,n_occurrences,mean_nll) -- mean NLL
## of each CDI target token at log-spaced training steps.
## Same 4-PL fit as the ladder: y = lo + (up-lo)/(1+exp((x-mid)/sc)) on x=log10(step);
## kappa = 0.434/sc (logit per natural-log experience), C&B convention.
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2)})

args <- commandArgs(trailingOnly = TRUE)
csv <- ifelse(length(args) >= 1, args[1], "/tmp/ar_pilot/word_surprisal.csv")
tag <- ifelse(length(args) >= 2, args[2], "5min")

d <- read_csv(csv, show_col_types = FALSE) |>
  mutate(step = as.integer(step), mean_nll = as.numeric(mean_nll), l10 = log10(step))
cat(sprintf("loaded %d rows, %d words, steps: %s\n", nrow(d), n_distinct(d$word),
            paste(sort(unique(d$step)), collapse = " ")))

four_pl <- function(x, y) tryCatch({
  f <- nls(y ~ lo + (up - lo) / (1 + exp((x - mid) / sc)),
           start = list(up = max(y), lo = min(y), mid = mean(x), sc = 0.5),
           lower = c(up = min(y) - 5, lo = min(y) - 5, mid = min(x) - 3, sc = 1e-3),
           upper = c(up = max(y) + 5, lo = max(y) + 5, mid = max(x) + 3, sc = 50),
           algorithm = "port", control = nls.control(maxiter = 200, warnOnly = TRUE))
  co <- coef(f); c(sc = unname(co["sc"]), rng = unname(co["up"] - co["lo"]))
}, error = function(e) c(sc = NA_real_, rng = NA_real_))

res <- d |> group_by(word) |> filter(n() >= 6) |>
  group_modify(~{ p <- four_pl(.x$l10, .x$mean_nll)
                  tibble(sc = p["sc"], rng = p["rng"], npts = nrow(.x),
                         drop = max(.x$mean_nll) - min(.x$mean_nll)) }) |>
  ungroup() |>
  mutate(kappa = 0.434 / sc)

pass <- res |> filter(is.finite(kappa), sc > 0.01, sc < 10, rng > 1)
cat(sprintf("\n=== %s ===\n", tag))
cat(sprintf("words fit: %d / %d ; passing C&B filter (range>1 nat): %d (%.0f%%)\n",
            sum(is.finite(res$kappa)), nrow(res), nrow(pass), 100 * nrow(pass) / nrow(res)))
cat(sprintf("median surprisal drop across words: %.2f nats\n", median(res$drop, na.rm = TRUE)))
if (nrow(pass) >= 5) {
  cat("kappa quantiles (5/25/50/75/95):\n"); print(round(quantile(pass$kappa, c(.05,.25,.5,.75,.95)), 3))
  cat(sprintf("median kappa = %.3f ; frac kappa>1 = %.2f ; frac>2 = %.2f\n",
              median(pass$kappa), mean(pass$kappa > 1), mean(pass$kappa > 2)))
}
write_csv(res, sub("\\.csv$", paste0("_kappa_", tag, ".csv"), csv))

if (nrow(pass) >= 5) {
  g <- ggplot(pass, aes(kappa)) +
    geom_density(fill = "#41ab5d", alpha = 0.4, colour = "#2c7fb8") +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    coord_cartesian(xlim = c(0, max(4, quantile(pass$kappa, .97)))) +
    labs(title = sprintf("ClimbMix GPT-2 per-CDI-word kappa (%s, n=%d words)", tag, nrow(pass)),
         subtitle = "training-step axis; kappa=1 dashed = unit accumulator",
         x = "kappa (0.434 / sigmoid scale)", y = "density") +
    theme_minimal(base_size = 11)
  ggsave(sub("\\.csv$", paste0("_kappa_", tag, ".png"), csv), g, width = 7, height = 4, dpi = 150)
  cat("wrote density plot\n")
}
