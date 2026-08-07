#!/usr/bin/env Rscript
# EARLY vs LATE-LEARNED WORDS on the CHILDES models: does a word's ONSET (when it is
# acquired) predict its KAPPA (how fast it is acquired)? Positive onset->kappa, controlling
# frequency = "later-learned words are learned faster" = the child acceleration signature.
# Children show this strongly; the question is whether the LMs do. Reuses the per-word sigmoid
# fits already behind Fig 6 (checkpoint axis) + fits the ladder (development axis). No training.
#
# onset = sigmoid midpoint (log10 of experience: training step OR data budget)
# kappa = 0.434 / scale   (logit per natural-log unit, paper convention)
# Truncation control: late-onset words whose transition runs off the end of the window have an
# unobserved lower asymptote -> inflated kappa -> spurious positive slope. We flag words whose
# observed floor sits >1 nat above the fitted asymptote and report results with them dropped.
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2); library(purrr)})
ROOT <- "/Users/mcfrank/Projects/standard_model_2"
OUT  <- file.path(ROOT, "studies/llm")

# ---- frequency + lexical class (CHILDES) ----
freq <- read_csv(file.path(ROOT, "figs/longitudinal/psi_freq_regression_per_word.csv"),
                 show_col_types = FALSE) |>
  transmute(word = item, log_freq = log_p, lexical_class)

report <- function(df, axis) {
  # df: word, onset, kappa, range, gap, log_freq
  base <- df |> filter(is.finite(kappa), is.finite(onset), is.finite(log_freq), range > 1)
  trim <- base |> filter(gap <= 1)                       # truncation-controlled
  stat <- function(d, tag) {
    n <- nrow(d); rc <- cor(d$kappa, d$onset)
    pc <- cor(resid(lm(kappa ~ log_freq, d)), resid(lm(onset ~ log_freq, d)))  # partial | freq
    co <- summary(lm(kappa ~ onset + log_freq, d))$coef
    cat(sprintf("  [%s | %s] n=%d  cor(k,onset)=%+.3f  partial(k,onset|freq)=%+.3f  onset_beta=%+.3f (p=%.3g)  median_kappa=%.2f\n",
                axis, tag, n, rc, pc, co["onset","Estimate"], co["onset","Pr(>|t|)"], median(d$kappa)))
    tibble(axis=axis, set=tag, n=n, cor_k_onset=rc, partial_k_onset_freq=pc,
           onset_beta=co["onset","Estimate"], onset_p=co["onset","Pr(>|t|)"], median_kappa=median(d$kappa))
  }
  bind_rows(stat(base,"all"), stat(trim,"asymptote-reached")) |> mutate(.plot = list(base))
}

# ===================== (1) CHECKPOINT AXIS (training steps) =====================
ckpt <- map_dfr(c(0,42,123), function(s){
  read.delim(file.path(ROOT, sprintf("fits/llm/sigmoids/gpt2_childes_seed%d_sigmoids.txt", s))) |>
    transmute(word = Token, seed = s, onset = ParamXmid, scale = ParamScale,
              range = ParamUpper - ParamLower, gap = MinSurprisal - ParamLower) |>
    filter(scale > 0.01, scale < 10)
}) |> mutate(kappa = 0.434/scale) |> inner_join(freq, by = "word")

# ===================== (2) DEVELOPMENT AXIS (data-budget ladder) =====================
fit1 <- function(x, y) tryCatch({
  f <- nls(y ~ lo + (up-lo)/(1+exp((x-mid)/sc)),
           start=list(up=max(y), lo=min(y), mid=mean(x), sc=0.3),
           lower=c(up=min(y)-5, lo=min(y)-5, mid=min(x)-2, sc=1e-3),
           upper=c(up=max(y)+5, lo=max(y)+5, mid=max(x)+2, sc=50),
           algorithm="port", control=nls.control(maxiter=200, warnOnly=TRUE))
  co <- coef(f); tibble(onset=co["mid"], scale=co["sc"], range=co["up"]-co["lo"], gap=min(y)-co["lo"])
}, error=function(e) tibble(onset=NA,scale=NA,range=NA,gap=NA))

ladder <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"), show_col_types=FALSE) |>
  mutate(lx = log10(words))
dev <- ladder |> group_by(seed, word) |> filter(n() >= 6) |>
  group_modify(~fit1(.x$lx, .x$surprisal)) |> ungroup() |>
  filter(is.finite(scale), scale > 0.01, scale < 10) |>
  mutate(kappa = 0.434/scale) |> inner_join(freq, by = "word")

cat("=== onset -> kappa (later-learned words learned faster?) : POSITIVE = kid-like acceleration ===\n")
res <- bind_rows(report(ckpt,"checkpoint (training step)"), report(dev,"development (data budget)"))
write_csv(res |> select(-.plot), file.path(OUT,"onset_kappa_childes_stats.csv"))

# ---- figure: kappa vs onset, both axes ----
pd <- bind_rows(
  ckpt |> filter(range>1, gap<=1) |> transmute(axis="checkpoint (log10 training step)", onset, kappa, log_freq),
  dev  |> filter(range>1, gap<=1) |> transmute(axis="development (log10 data budget)",  onset, kappa, log_freq))
g <- ggplot(pd, aes(onset, kappa)) +
  geom_point(aes(colour=log_freq), alpha=.5, size=1) +
  geom_smooth(method="lm", se=TRUE, colour="black", linewidth=.6) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  scale_colour_viridis_c(name="log CHILDES\nfrequency") +
  facet_wrap(~axis, scales="free_x") + coord_cartesian(ylim=c(0,3)) +
  labs(title="Do CHILDES LMs learn their later-acquired words faster?",
       subtitle="onset = when a word is acquired; kappa = how fast. POSITIVE slope = acceleration (children); flat/negative = none/diminishing",
       x="onset (sigmoid midpoint)", y="kappa = 0.434 / scale") +
  theme_minimal(base_size=11)
ggsave(file.path(OUT,"onset_kappa_childes.png"), g, width=9, height=4.4, dpi=140)
cat("\nwrote onset_kappa_childes.png + onset_kappa_childes_stats.csv\n")
