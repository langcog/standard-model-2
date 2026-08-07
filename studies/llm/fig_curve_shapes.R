#!/usr/bin/env Rscript
# Does the Chang&Bergen sigmoid hide within-curve acceleration (e.g. a "grokking" slow-then-sudden
# shape)? Diagnostic on the TRAINING-STEP axis (73 log-spaced checkpoints/word -- where grokking would
# live). For selected words we overlay the raw data, the fitted 4-PL SIGMOID, and a shape-agnostic
# LOESS. If the loess reveals a slow-then-sudden shape the sigmoid misses, we should fit richer shapes.
# Words are chosen to INCLUDE the most acceleration-like curves (most concave / negative curvature).
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(tidyr)})
ROOT <- "/Users/mcfrank/Projects/standard_model_2"; OUT <- file.path(ROOT,"studies/llm")

d <- read_csv(file.path(ROOT,"fits/llm/surprisal_gpt2_childes_seed42.csv"), show_col_types=FALSE) |>
  filter(n_occurrences>=20) |> mutate(lx=log10(step))
# per-word: quadratic curvature b2 (on surprisal ~ lx + lx^2; b2<0 = concave = accelerating improvement)
# + total drop, to keep words that actually learn something.
wm <- d |> group_by(word) |> filter(n()>=10) |>
  summarise(b2 = coef(lm(mean_nll~lx+I(lx^2)))[3], drop = first(mean_nll[order(lx)])-last(mean_nll[order(lx)]),
            .groups="drop") |> filter(is.finite(b2), drop>2)
cat(sprintf("CHECKPOINT axis: n=%d words  frac decelerating(b2>0)=%.3f  accelerating(b2<0)=%.3f\n",
            nrow(wm), mean(wm$b2>0), mean(wm$b2<0)))

# select: 4 most 'accelerating' (b2 most negative) + 4 most 'decelerating' (b2 most positive)
acc <- wm |> arrange(b2) |> slice(1:4) |> mutate(grp="most acceleration-like (b2<0)")
dec <- wm |> arrange(desc(b2)) |> slice(1:4) |> mutate(grp="typical: diminishing (b2>0)")
sel <- bind_rows(acc,dec)
ds <- d |> inner_join(sel, by="word") |>
  mutate(lab=sprintf("%s  (b2=%+.2f)", word, b2),
         grp=factor(grp, levels=c("most acceleration-like (b2<0)","typical: diminishing (b2>0)")))

# fit 4-PL sigmoid per selected word for the overlay
sigfit <- function(df){ x<-df$lx; y<-df$mean_nll
  f<-tryCatch(nls(y~lo+(up-lo)/(1+exp((x-mid)/sc)),
    start=list(up=max(y),lo=min(y),mid=mean(x),sc=.5),lower=c(min(y)-5,min(y)-5,min(x)-3,1e-3),
    upper=c(max(y)+5,max(y)+5,max(x)+3,50),algorithm="port",control=nls.control(maxiter=300,warnOnly=TRUE)),error=function(e)NULL)
  if(is.null(f)) return(NULL)
  gx<-seq(min(x),max(x),length=120); tibble(lx=gx, fit=predict(f,newdata=list(x=gx)))}
sigc <- ds |> group_by(word,lab,grp) |> group_modify(~sigfit(.x)) |> ungroup()

fitpal <- c("4-PL sigmoid"="#d95f02", "loess (shape-agnostic)"="#1b9e77")
p <- ggplot(ds, aes(10^lx, mean_nll)) +
  geom_point(colour="grey55", size=.9, alpha=.8) +
  geom_line(data=sigc, aes(10^lx, fit, colour="4-PL sigmoid"), linewidth=.8) +
  geom_smooth(aes(colour="loess (shape-agnostic)"), method="loess", se=FALSE, span=.6, linewidth=.8) +
  facet_wrap(~grp+lab, ncol=4, scales="free_y") +
  scale_x_log10() + scale_colour_manual(values=fitpal, name=NULL) +
  labs(title="Do individual word-learning curves hide a 'grokking' shape the sigmoid misses?",
       subtitle="CHILDES LM, training-step axis (73 checkpoints). Raw points; 4-PL sigmoid vs shape-agnostic loess. Top row = the MOST acceleration-like words.",
       x="training step (log)", y="surprisal (nats)") +
  theme_minimal(base_size=10) + theme(legend.position="top", strip.text=element_text(size=7.5))
ggsave(file.path(OUT,"fig_curve_shapes.png"), p, width=11, height=6, dpi=140)
cat("wrote fig_curve_shapes.png\n")
