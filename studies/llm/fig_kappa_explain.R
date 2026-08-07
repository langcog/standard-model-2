#!/usr/bin/env Rscript
# Pedagogical figure: HOW the data-axis kappa-decline is measured.
# A: the model's competence curve (mean CDI-word surprisal) vs training data. We measure the LOCAL
#    slope (improvement per 10x data) early vs late -- it shrinks (2.8 -> 1.1 nats/decade) = the curve
#    flattens = diminishing returns. A straight dashed line = a constant-kappa=1 accumulator for reference.
# B: the same slope expressed as effective kappa (=0.434 x nats/decade), declining 1.2 -> 0.5, with the
#    unit-accumulator line (kappa=1) and children (kappa~10, constant) for contrast.
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(patchwork)})
ROOT<-"/Users/mcfrank/Projects/standard_model_2"; OUT<-file.path(ROOT,"studies/llm"); LM<-"#2c7fb8"

lad <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"),show_col_types=FALSE)
bud <- lad |> group_by(words) |> summarise(surp=mean(surprisal),.groups="drop") |>
  mutate(lx=log10(words)) |> arrange(lx)
q <- lm(surp~lx+I(lx^2), bud); b<-coef(q)
kap <- function(x) 0.434*-(b[2]+2*b[3]*x)          # smooth effective kappa (for the curve in B)
lo<-min(bud$lx); hi<-max(bud$lx)
# ROBUST regional slopes: local linear fit over the low third and high third of the budget range
early_s <- coef(lm(surp~lx, bud |> filter(lx <= quantile(lx,1/3))))[2]   # nats per decade (neg)
late_s  <- coef(lm(surp~lx, bud |> filter(lx >= quantile(lx,2/3))))[2]
kap_e <- 0.434*-early_s; kap_l <- 0.434*-late_s
cat(sprintf("early (<=~1.5M): drop/decade=%.2f nats  kappa=%.2f\n", -early_s, kap_e))
cat(sprintf("late (>=~8M):    drop/decade=%.2f nats  kappa=%.2f\n", -late_s, kap_l))
cat(sprintf("quad b2=%.3f (p=%.2g); per-word frac b2>0 (diminishing)=0.84\n", summary(q)$coef[3,1], summary(q)$coef[3,4]))
slope <- function(x) b[2]+2*b[3]*x

xb<-c(5e5,1e6,3e6,1e7,2.4e7); xl<-c("0.5M","1M","3M","10M","24M")
pred <- tibble(lx=seq(lo,hi,length=100)) |> mutate(words=10^lx, surp=predict(q,across()))
# slope segments drawn with the ROBUST regional slopes, anchored on the smooth curve
tri <- function(x0,s,half=0.45){ y0<-b[1]+b[2]*x0+b[3]*x0^2
  tibble(x=c(x0-half,x0+half), y=c(y0-s*half, y0+s*half)) }
te<-tri(lo+0.35, early_s); tl<-tri(hi-0.35, late_s)
# constant kappa=1 reference: slope -1/0.434 = -2.30 nats/decade, anchored at the first point
kref <- tibble(lx=c(lo,hi)) |> mutate(surp=bud$surp[1] + (-1/0.434)*(lx-lo), words=10^lx)

pA <- ggplot(bud, aes(words,surp)) +
  geom_line(data=kref, aes(words,surp), linetype="dashed", colour="grey60") +
  annotate("text", x=1.1e6, y=kref$surp[1]-2.0, angle=-26, size=2.7, colour="grey55",
           label="constant returns (kappa = 1)") +
  geom_line(data=pred, aes(words,surp), colour=LM, linewidth=.7) + geom_point(size=1.8, colour=LM) +
  geom_line(data=te, aes(10^x,y), colour="#d95f02", linewidth=1) +
  geom_line(data=tl, aes(10^x,y), colour="#d95f02", linewidth=1) +
  annotate("text", x=10^(lo+0.35), y=b[1]+b[2]*(lo+0.35)+b[3]*(lo+0.35)^2-0.9, size=3, colour="#d95f02",
           label=sprintf("early: -%.1f nats\nper 10x data", -early_s)) +
  annotate("text", x=10^(hi-0.55), y=b[1]+b[2]*(hi-0.35)+b[3]*(hi-0.35)^2+1.1, size=3, colour="#d95f02",
           label=sprintf("late: -%.1f nats\nper 10x data", -late_s)) +
  scale_x_log10(breaks=xb,labels=xl) +
  labs(subtitle="A. The model's competence curve flattens: each 10x of data buys less",
       x="training data (distinct words)", y="mean CDI-word surprisal\n(lower = more competent)") +
  theme_minimal(base_size=11)+theme(plot.subtitle=element_text(face="bold"))

kd <- tibble(lx=seq(lo,hi,length=100)) |> mutate(words=10^lx, k=kap(lx))
pB <- ggplot(kd, aes(words,k)) + geom_line(colour=LM, linewidth=.9) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey55") +
  annotate("text", x=6e5, y=1.08, hjust=0, size=2.9, colour="grey45", label="unit accumulator (kappa=1)") +
  annotate("text", x=6e5, y=0.35, hjust=0, size=3.1, colour=LM, fontface="bold",
           label=sprintf("kappa falls %.1f -> %.1f", kap_e, kap_l)) +
  annotate("segment", x=7e6, xend=7e6, y=2.1, yend=2.35, colour="#c41e37",
           arrow=arrow(length=unit(.12,"cm"))) +
  annotate("text", x=7e6, y=2.0, vjust=0, size=3, colour="#c41e37", label="children: kappa~10, constant") +
  scale_x_log10(breaks=xb,labels=xl) + coord_cartesian(ylim=c(0.2,2.4)) +
  labs(subtitle="B. Same thing as effective kappa: it DECLINES (diminishing returns)",
       x="training data (distinct words)", y="effective kappa\n(returns per e-fold of experience)") +
  theme_minimal(base_size=11)+theme(plot.subtitle=element_text(face="bold"))

ggsave(file.path(OUT,"fig_kappa_explain.png"), pA|pB, width=11, height=4.4, dpi=140)
cat("wrote fig_kappa_explain.png\n")
