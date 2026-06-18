#!/usr/bin/env Rscript
# Finalize plots for the acceleration campaign. Reads bank_results.tsv, writes 3 PNGs.
suppressWarnings(suppressMessages(library(ggplot2)))
d <- read.delim("bank_results.tsv", stringsAsFactors = FALSE)
names(d)[names(d) == "final_surp"] <- "surp"
d$control <- d$r_early >= 0.4 & d$surp <= 6.1   # genuine = passes both artifact controls

# ---- family assignment ----
fam <- function(l) {
  if (grepl("^lr_0p[0-9]+$|^constlr$|^lr_[24]$", l)) return("LR-flat (magnitude)")
  if (grepl("init0p3", l)) return("LR x small-init")
  if (grepl("^rep[0-9]|^s1_rep", l)) return("repetition")
  if (grepl("^gf", l)) return("grokfast")
  if (grepl("^arch_", l)) return("depth")
  if (grepl("^batch_", l)) return("batch")
  if (grepl("^matlr_", l)) return("matrix-LR")
  if (grepl("^decaylow_", l)) return("decay x low")
  if (grepl("^warm|^finallr|^wsd", l)) return("schedule")
  if (grepl("^betas|^window", l)) return("optim/attn")
  if (grepl("^init_", l)) return("init-scale")
  if (grepl("^wd_", l)) return("weight-decay")
  if (grepl("^baseline", l)) return("baseline")
  "other"
}
d$family <- vapply(d$label, fam, "")

# ---- (1) J vs LR : the money plot ----
mult <- c(lr_0p03=.03, lr_0p05=.05, lr_0p07=.07, lr_0p1=.1, lr_0p15=.15,
          lr_0p25=.25, lr_0p35=.35, lr_0p5=.5, constlr=1, lr_2=2, lr_4=4)
lr <- d[d$label %in% names(mult), ]
lr$mult <- mult[lr$label]
lr$regime <- ifelse(lr$surp > 6.1, "undertrain (artifact)",
              ifelse(lr$J > 0, "genuine increasing", "diminishing"))
p1 <- ggplot(lr, aes(mult, J)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line(colour = "grey70") +
  geom_point(aes(colour = surp, shape = regime), size = 4) +
  scale_x_log10(breaks = c(.03,.05,.1,.25,.5,1,2,4)) +
  scale_colour_viridis_c(option = "C", name = "final\nsurprisal") +
  labs(title = "Flat-LR dose-response: lower LR -> increasing returns (J)",
       subtitle = "genuine frontier x0.07-0.1 (surp ~ baseline); below x0.07 competence breaks",
       x = "learning-rate multiplier (flat, log scale)", y = "J = late - early gain-rate") +
  theme_minimal(base_size = 12)
ggsave("fig_J_vs_LR.png", p1, width = 8, height = 5, dpi = 130)

# ---- (2) J vs final_surp : artifact-diagnostic quadrant ----
p2 <- ggplot(d, aes(surp, J, colour = control)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 6.1, linetype = 2, colour = "grey50") +
  geom_point(alpha = .8, size = 2) +
  scale_colour_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02"),
                      name = "passes controls\n(r_early>=.4 & surp<=6.1)") +
  annotate("text", x = 5.75, y = max(d$J), label = "genuine hits", colour = "#1b9e77", hjust = 0) +
  labs(title = "Artifact diagnostic: J vs competence",
       subtitle = "genuine increasing returns = upper-left (J>0 AND surp<=6.1); J>0 with surp>6.1 = slow-start/undertrain artifact",
       x = "final mean surprisal (lower = more competent; baseline 5.88)", y = "J") +
  theme_minimal(base_size = 12)
ggsave("fig_J_vs_finalsurp.png", p2, width = 8, height = 5, dpi = 130)

# ---- (3) best within-control J per family ----
agg <- do.call(rbind, lapply(split(d, d$family), function(g) {
  gc <- g[g$control, ]
  bestc <- if (nrow(gc)) max(gc$J) else NA_real_
  data.frame(family = g$family[1], bestJ_control = bestc, bestJ_raw = max(g$J))
}))
agg <- agg[order(-replace(agg$bestJ_control, is.na(agg$bestJ_control), -99)), ]
agg$family <- factor(agg$family, levels = agg$family)
agg$plotJ <- ifelse(is.na(agg$bestJ_control), agg$bestJ_raw, agg$bestJ_control)
agg$passes <- ifelse(is.na(agg$bestJ_control), "no control-passing variant", "best genuine (control-passing)")
p3 <- ggplot(agg, aes(family, plotJ, fill = passes)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_col() +
  scale_fill_manual(values = c(`best genuine (control-passing)` = "#1b9e77",
                               `no control-passing variant` = "grey70"), name = NULL) +
  coord_flip() +
  labs(title = "Best J per family (within artifact controls where possible)",
       subtitle = "only the LR-flat and LR x small-init families produce genuine increasing returns",
       x = NULL, y = "best J") +
  theme_minimal(base_size = 12)
ggsave("fig_bestJ_by_family.png", p3, width = 8, height = 5, dpi = 130)

cat("plots written: fig_J_vs_LR.png, fig_J_vs_finalsurp.png, fig_bestJ_by_family.png\n")
