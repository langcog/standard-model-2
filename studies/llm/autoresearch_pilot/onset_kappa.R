## Flavor-(ii) test: does a word's ONSET (acquisition step) predict its kappa
## (steepness)?  Positive onset->kappa, controlling frequency = "later words
## learned faster" = learning-to-learn signature.
## CAVEAT: single-pass/undertrained data -> late-onset words are more TRUNCATED
## (transition near the end, lower asymptote unobserved) which can upward-bias
## their kappa and manufacture a spurious positive slope. Read accordingly; the
## clean test is a converged late-held-out-lexicon run.
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2)})

d <- read_csv("/tmp/ar_pilot/trajectory_singlepass.csv", show_col_types = FALSE) |>
  mutate(step = as.integer(step), mean_nll = as.numeric(mean_nll), l10 = log10(step))

fp <- function(x, y) tryCatch({
  f <- nls(y ~ lo + (up - lo) / (1 + exp((x - mid) / sc)),
           start = list(up = max(y), lo = min(y), mid = mean(x), sc = 0.5),
           lower = c(up = min(y)-5, lo = min(y)-5, mid = min(x)-3, sc = 1e-3),
           upper = c(up = max(y)+5, lo = max(y)+5, mid = max(x)+3, sc = 50),
           algorithm = "port", control = nls.control(maxiter = 200, warnOnly = TRUE))
  co <- coef(f); c(sc = unname(co["sc"]), mid = unname(co["mid"]), rng = unname(co["up"] - co["lo"]))
}, error = function(e) c(sc = NA, mid = NA, rng = NA))

res <- d |> group_by(word) |> filter(n() >= 6) |>
  group_modify(~{ p <- fp(.x$l10, .x$mean_nll); tibble(sc = p["sc"], mid = p["mid"], rng = p["rng"]) }) |>
  ungroup() |> mutate(kappa = 0.434 / sc)

fr <- read_csv("/tmp/ar_pilot/freq.csv", show_col_types = FALSE) |> mutate(logfreq = log10(count))
m <- res |> inner_join(fr, by = "word") |>
  filter(is.finite(kappa), sc > 0.01, sc < 10, rng > 1, count > 0)

cat("n words:", nrow(m), " | onset(log10 step) range:", round(range(m$mid), 2), "\n\n")
cat("cor(kappa, onset)      :", round(cor(m$kappa, m$mid), 3), "\n")
cat("cor(onset, logfreq)    :", round(cor(m$mid, m$logfreq), 3), "  (freq drives onset earlier)\n")
cat("cor(kappa, logfreq)    :", round(cor(m$kappa, m$logfreq), 3), "\n")
r_ko <- cor(resid(lm(kappa ~ logfreq, m)), resid(lm(mid ~ logfreq, m)))
cat("partial cor(kappa,onset | logfreq):", round(r_ko, 3), "  <- learning-to-learn test\n\n")
cat("== kappa ~ onset + logfreq ==\n"); print(round(summary(lm(kappa ~ mid + logfreq, m))$coef, 4))

g <- ggplot(m, aes(mid, kappa, color = logfreq)) +
  geom_point(alpha = .6) + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = .6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_color_viridis_c(name = "log10 freq") +
  labs(title = "Does later onset -> steeper kappa? (single-pass pilot; truncation-confounded)",
       subtitle = "positive slope = 'later words learned faster' = learning-to-learn signature",
       x = "onset = sigmoid midpoint (log10 training step)", y = "kappa (0.434 / scale)") +
  theme_minimal(base_size = 11)
ggsave("/tmp/ar_pilot/onset_kappa.png", g, width = 7.5, height = 4.5, dpi = 150)
write_csv(m, "/tmp/ar_pilot/onset_kappa_fits.csv")
cat("\nwrote onset_kappa.png + fits\n")
