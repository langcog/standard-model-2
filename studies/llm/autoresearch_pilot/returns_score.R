## Returns-to-experience scorer for a single autoresearch run's probe trajectory.
## Objective J = (improvement rate late) - (improvement rate early): how the
## learner's per-decade gain CHANGES with experience. J<0 = diminishing returns
## (decelerating; the default LM regime); J->0 or >0 = the kid-like, learning-to-
## learn direction we hunt for. Also reports the aggregate quadratic curvature and,
## as a (demoted) diagnostic, the median per-word slope kappa_w.
##   Rscript returns_score.R <word_surprisal.csv> [label] [append_tsv]
suppressPackageStartupMessages({library(dplyr); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
csv <- args[1]
label <- if (length(args) >= 2) args[2] else basename(dirname(normalizePath(csv)))
out <- if (length(args) >= 3) args[3] else ""

d <- read_csv(csv, show_col_types = FALSE) |>
  mutate(step = as.integer(step), mean_nll = as.numeric(mean_nll), lx = log10(step))
agg <- d |> group_by(step, lx) |> summarise(s = mean(mean_nll), .groups = "drop") |> arrange(lx)
mid <- median(agg$lx)
r_e <- as.numeric(-coef(lm(s ~ lx, agg |> filter(lx <= mid)))[2])   # nats improvement / decade, early
r_l <- as.numeric(-coef(lm(s ~ lx, agg |> filter(lx > mid)))[2])    # late
J   <- r_l - r_e
b2  <- as.numeric(coef(lm(s ~ lx + I(lx^2), agg))[3])
final_s <- tail(agg$s, 1)
pw <- d |> group_by(word) |> filter(n() >= 6) |>
  summarise(sl = as.numeric(coef(lm(mean_nll ~ lx))[2]), .groups = "drop")
kappa_w_med <- median(-pw$sl * 0.434, na.rm = TRUE)   # logit/ln-unit, diagnostic only

hdr <- "label\tr_early\tr_late\tJ\tquad_b2\tfinal_surp\tkappa_w_med\tn_steps"
line <- sprintf("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d",
                label, r_e, r_l, J, b2, final_s, kappa_w_med, nrow(agg))
cat(hdr, "\n", line, "\n", sep = "")
cat(sprintf("  -> %s returns (J=%.3f, b2=%.3f): %s\n", label, J, b2,
            ifelse(J < -0.05 | b2 > 0.05, "DIMINISHING (decelerating)",
                   ifelse(J > 0.05, "INCREASING (accelerating!)", "~flat"))))
if (nchar(out) > 0) {
  if (!file.exists(out)) cat(hdr, "\n", file = out)
  cat(line, "\n", file = out, append = TRUE)
}
