#!/usr/bin/env Rscript
# Intermediate / intuition figures for the LM acceleration story (standalone; not wired to paper):
#  (1) what do individual-word surprisal curves look like? (development/ladder axis)
#  (2) how does LM word-difficulty relate to children's word difficulty (delta_j, bayes_long M3)?
suppressPackageStartupMessages({library(dplyr);library(readr);library(tidyr);library(ggplot2);library(purrr);library(patchwork)})
ROOT <- "/Users/mcfrank/Projects/standard_model_2"; OUT <- file.path(ROOT,"studies/llm")
LM_BLUE <- "#2c7fb8"; CHILD_RED <- "#c41e37"

# ---- data ----
lad <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"), show_col_types=FALSE) |> mutate(lx=log10(words))
curve <- lad |> group_by(word, words, lx) |> summarise(surp=mean(surprisal), .groups="drop")   # mean over 10 seeds
cd <- read_csv(file.path(ROOT,"fits/bayes_long/summaries/marchman_a3_m3_psi.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item))) |> filter(!grepl(" ",word)) |> select(word, delta_j)
fq <- read_csv(file.path(ROOT,"data/intermediates/cdi_master_item_key.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item)), logfreq=log10(prob)) |>
  filter(is.finite(logfreq)) |> select(word, logfreq, lexical_class=lexical_category) |>
  distinct(word, .keep_all=TRUE)

fitmid <- function(x,y) tryCatch(coef(nls(y~lo+(up-lo)/(1+exp((x-mid)/sc)),
  start=list(up=max(y),lo=min(y),mid=mean(x),sc=.3),lower=c(min(y)-5,min(y)-5,min(x)-2,1e-3),
  upper=c(max(y)+5,max(y)+5,max(x)+2,50),algorithm="port",control=nls.control(maxiter=200,warnOnly=TRUE)))["mid"],
  error=function(e)NA_real_)
lmw <- curve |> group_by(word) |> summarise(final_surp=surp[which.max(lx)], onset=fitmid(lx,surp), .groups="drop")
d <- cd |> inner_join(lmw,by="word") |> left_join(fq,by="word") |> filter(is.finite(final_surp), is.finite(delta_j))

xbrk <- c(5e5,1e6,3e6,1e7,2.4e7); xlab <- c("0.5M","1M","3M","10M","24M")

# ========== FIG 1: individual word surprisal curves ==========
# pick 9 words spanning children's difficulty (delta_j), ordered easy->hard for kids
sel <- d |> arrange(delta_j) |> slice(round(seq(1, n(), length.out=9))) |> pull(word)
c1 <- curve |> filter(word %in% sel) |>
  left_join(d |> select(word, delta_j, logfreq), by="word") |>
  mutate(lab = sprintf("%s  (child δ=%+.1f, logfreq=%.1f)", word, delta_j, logfreq),
         lab = factor(lab, levels = d |> filter(word%in%sel) |> arrange(delta_j) |>
                        mutate(lab=sprintf("%s  (child δ=%+.1f, logfreq=%.1f)",word,delta_j,logfreq)) |> pull(lab)))
p1 <- ggplot(c1, aes(words, surp)) +
  geom_line(colour=LM_BLUE, linewidth=.7) + geom_point(colour=LM_BLUE, size=1.1) +
  facet_wrap(~lab, ncol=3) +
  scale_x_log10(breaks=xbrk, labels=xlab) +
  labs(title="(1) Individual word surprisal curves (CHILDES LM, development axis)",
       subtitle="mean over 10 seeds; panels ordered by children's difficulty (δ). surprisal falls as the model sees more distinct words.",
       x="training data (distinct words)", y="surprisal (nats)") +
  theme_minimal(base_size=10) + theme(strip.text=element_text(size=8))
ggsave(file.path(OUT,"fig_word_curves.png"), p1, width=9.5, height=6.2, dpi=140)

# ========== FIG 2: LM vs child difficulty, and the frequency dissociation ==========
rraw <- cor(d$final_surp, d$delta_j)
df <- d |> filter(is.finite(logfreq))     # words with a CHILDES frequency (near-full)
rpar <- cor(resid(lm(final_surp~logfreq,df)), resid(lm(delta_j~logfreq,df)))
pA <- ggplot(df, aes(delta_j, final_surp)) +
  geom_point(aes(colour=logfreq), alpha=.7, size=1.4) +
  geom_smooth(method="lm", se=TRUE, colour="grey25", linewidth=.5) +
  scale_colour_viridis_c(option="D", name="log CHILDES\nfrequency") +
  labs(title="(2A) LM vs child word difficulty",
       subtitle=sprintf("r = %.2f raw;  %.2f controlling frequency", rraw, rpar),
       x="child difficulty  δ  (bayes_long M3; higher = harder for children)",
       y="LM difficulty (surprisal at 24M words)") + theme_minimal(base_size=11)

# frequency dissociation: z-scored difficulty vs frequency, LM vs child (one axis, two series)
z <- function(v) as.numeric(scale(v))
dd <- df |> transmute(logfreq, LM=z(final_surp), Children=z(delta_j)) |>
  pivot_longer(c(LM,Children), names_to="who", values_to="difficulty")
rlm <- cor(df$logfreq,df$final_surp); rch <- cor(df$logfreq,df$delta_j)
pB <- ggplot(dd, aes(logfreq, difficulty, colour=who)) +
  geom_point(alpha=.25, size=.9) + geom_smooth(method="lm", se=FALSE, linewidth=1) +
  scale_colour_manual(values=c(LM=LM_BLUE, Children=CHILD_RED), name=NULL) +
  annotate("text", x=-3.1, y=2.0, hjust=0, size=3, colour=LM_BLUE, fontface="bold",
           label=sprintf("LM: r=%.2f\n(frequency drives it)", rlm)) +
  annotate("text", x=-3.1, y=-1.4, hjust=0, size=3, colour=CHILD_RED, fontface="bold",
           label=sprintf("Children: r=%.2f\n(~orthogonal)", rch)) +
  labs(title="(2B) Why they only partly agree: frequency",
       subtitle="LM difficulty tracks word frequency; children's does not",
       x="log CHILDES frequency", y="word difficulty (z-scored)") + theme_minimal(base_size=11)

ggsave(file.path(OUT,"fig_lm_child_difficulty.png"), pA + pB, width=11, height=4.5, dpi=140)
cat(sprintf("n=%d  r(LM,child)=%.2f raw, %.2f partial | r(freq,LM)=%.2f r(freq,child)=%.2f\n",
            nrow(d), rraw, rpar, rlm, rch))
cat("wrote fig_word_curves.png + fig_lm_child_difficulty.png\n")
