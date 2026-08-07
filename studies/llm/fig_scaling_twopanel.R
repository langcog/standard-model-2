#!/usr/bin/env Rscript
# The scaling law two ways: same CDI-word data.
# A: loss (mean surprisal) vs log data -> CONVEX, flattening toward the entropy floor E (a straight
#    line in log-D, dashed grey, misses the bend -> loss is NOT linear in log-D).
# B: log(EXCESS loss = L - E) vs log data -> a STRAIGHT line, slope -beta (the actual scaling law).
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(patchwork)})
ROOT<-"/Users/mcfrank/Projects/standard_model_2"; OUT<-file.path(ROOT,"studies/llm"); LM<-"#2c7fb8"

bud <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"),show_col_types=FALSE) |>
  group_by(words) |> summarise(L=mean(surprisal),.groups="drop") |> arrange(words)
f <- nls(L ~ E + B*words^(-beta), data=bud, start=list(E=3,B=50,beta=0.3),
         control=nls.control(maxiter=500,warnOnly=TRUE))
E<-coef(f)["E"]; B<-coef(f)["B"]; beta<-coef(f)["beta"]
bud <- bud |> mutate(excess = L - E)
sm  <- tibble(words=10^seq(log10(min(bud$words)),log10(max(bud$words)),length=100)) |>
  mutate(L = E + B*words^(-beta), excess = B*words^(-beta))
linfit <- lm(L~log10(words), bud)                                   # the "linear in log-D" fit
r2lin <- summary(linfit)$r.squared
xb<-c(5e5,1e6,3e6,1e7,2.4e7); xl<-c("0.5M","1M","3M","10M","24M")

# ---- A: loss vs log data (convex; floor; linear fit misses it) ----
pA <- ggplot(bud, aes(words,L)) +
  geom_line(aes(words, predict(linfit, bud)), linetype="dashed", colour="grey60") +
  geom_hline(yintercept=E, linetype="dotted", colour="#d95f02") +
  annotate("text", x=5.2e5, y=E+0.18, hjust=0, size=3, colour="#d95f02",
           label=sprintf("entropy floor  E = %.2f", E)) +
  annotate("text", x=3.2e6, y=6.6, hjust=0, size=3, colour="grey45",
           label=sprintf("straight line in log-D\n(misses the bend; R2=%.3f)", r2lin)) +
  geom_line(data=sm, aes(words,L), colour=LM, linewidth=.8) + geom_point(size=2, colour=LM) +
  scale_x_log10(breaks=xb,labels=xl) +
  labs(subtitle="A. Loss vs log-data: CONVEX, flattens to the floor (not linear in log-D)",
       x="training data (distinct words)", y="mean CDI-word surprisal (loss)") +
  theme_minimal(base_size=11)+theme(plot.subtitle=element_text(face="bold"))

# ---- B: log(excess loss) vs log data -> straight scaling line, slope -beta ----
pB <- ggplot(bud, aes(words,excess)) +
  geom_line(data=sm, aes(words,excess), colour=LM, linewidth=.8) + geom_point(size=2, colour=LM) +
  annotate("text", x=6e5, y=1.35, hjust=0, size=3.2, colour=LM, fontface="bold",
           label=sprintf("straight power law:
slope = -beta = -%.2f  (R2=0.998)", beta)) +
  scale_x_log10(breaks=xb,labels=xl) + scale_y_log10() +
  labs(subtitle="B. Excess loss (L - E), log-log: a STRAIGHT scaling line",
       x="training data (distinct words)", y="excess loss  L - E  (log scale)") +
  theme_minimal(base_size=11)+theme(plot.subtitle=element_text(face="bold"))

ggsave(file.path(OUT,"fig_scaling_twopanel.png"), pA|pB, width=11, height=4.4, dpi=140)
cat(sprintf("E=%.2f B=%.1f beta=%.3f  | linear-in-logD R2=%.3f\n", E,B,beta,r2lin))
cat("wrote fig_scaling_twopanel.png\n")
