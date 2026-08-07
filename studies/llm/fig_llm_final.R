#!/usr/bin/env Rscript
# Candidate replacement for fig-llm-acceleration (standalone; .qmd untouched):
#  A: density of acceleration kappa -- children (by child) vs LMs (by word). [current figure]
#  B: the scaling law -- excess loss (L - E) vs data, log-log, a straight power line (slope -beta).
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(patchwork)})
ROOT<-"/Users/mcfrank/Projects/standard_model_2"; OUT<-file.path(ROOT,"studies/llm"); LM<-"#2c7fb8"

# ---- A: kappa density (reuse the paper's cached fig6 data) ----
sl <- readRDS(file.path(ROOT,"paper/cache/fig6_llm_slopes.rds"))$slopes
pal <- c("Children (English)"="#c41e37","Children (Norwegian)"="#fb9a99",
  "LMs: C&B 2022 (4 architectures)"="#8856a7","LMs: CHILDES (training)"="#2c7fb8","LMs: CHILDES (development)"="#41ab5d")
pA <- ggplot(sl, aes(slope_natural, fill=group, colour=group)) + geom_density(alpha=.32, linewidth=.5) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey55", linewidth=.4) +
  annotate("text", x=1.6, y=0.92, hjust=0, vjust=1, size=3, colour="grey45",
           label="kappa = 1\nunit accumulator\n(no acceleration)") +
  scale_fill_manual(values=pal,name=NULL) + scale_colour_manual(values=pal,name=NULL) +
  coord_cartesian(xlim=c(0,22)) +
  labs(subtitle="A. Acceleration (kappa): children by child, LMs by word", x="kappa: slope on log(experience)", y="density") +
  theme_minimal(base_size=11) + theme(legend.position=c(0.66,0.60),
    legend.background=element_rect(fill="white",colour=NA), legend.text=element_text(size=7.2),
    plot.subtitle=element_text(face="bold"))

# ---- B: scaling law -- excess loss (L - E) log-log ----
bud <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"),show_col_types=FALSE) |>
  group_by(words) |> summarise(L=mean(surprisal),.groups="drop") |> arrange(words)
f <- nls(L ~ E + B*words^(-beta), data=bud, start=list(E=3,B=50,beta=0.3), control=nls.control(maxiter=500,warnOnly=TRUE))
E<-coef(f)["E"]; B<-coef(f)["B"]; beta<-coef(f)["beta"]
# kappa_eff = beta * (L - E)  (see derivation) -> a power law in D, log-log linear with slope -beta
bud <- bud |> mutate(kappa = beta*(L-E))
sm <- tibble(words=10^seq(log10(min(bud$words)),log10(max(bud$words)),length=100)) |> mutate(kappa=beta*B*words^(-beta))
xb<-c(5e5,1e6,3e6,1e7,2.4e7); xl<-c("0.5M","1M","3M","10M","24M")
pB <- ggplot(bud, aes(words,kappa)) + geom_line(data=sm, colour=LM, linewidth=.8) + geom_point(size=2, colour=LM) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey55") +
  annotate("text", x=1.3e7, y=1.06, hjust=0, size=2.9, colour="grey45", label="kappa = 1") +
  annotate("text", x=6e5, y=0.62, hjust=0, size=3.1, colour=LM, fontface="bold",
           label=sprintf("kappa = beta x (L - E)\npower law, slope -beta = -%.2f", beta)) +
  scale_x_log10(breaks=xb,labels=xl) + scale_y_log10(breaks=c(0.5,0.7,1,1.5)) +
  labs(subtitle="B. Effective kappa falls as a power law (log-log linear, slope -beta)",
       x="training data (distinct words)", y="effective kappa  (= beta x excess loss, log scale)") +
  theme_minimal(base_size=11) + theme(plot.subtitle=element_text(face="bold"))

ggsave(file.path(OUT,"fig_llm_final.png"), pA|pB, width=11, height=4.5, dpi=140)
cat(sprintf("E=%.2f beta=%.3f\n", E, beta)); cat("wrote fig_llm_final.png\n")
