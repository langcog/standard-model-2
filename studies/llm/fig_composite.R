#!/usr/bin/env Rscript
# CANDIDATE COMPOSITE FIGURE (standalone; .qmd untouched) for the LM section.
#  A: per-word acquisition curves, children vs LM, both rising (child steep, LM shallow+negative)
#  B: LM vs child word difficulty, by lexical class
#  C: the kappa density (current fig-llm-acceleration / Fig 3)
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(tidyr);library(patchwork)})
ROOT <- "/Users/mcfrank/Projects/standard_model_2"; OUT <- file.path(ROOT,"studies/llm")
CHILD_RED <- "#c41e37"; LM_BLUE <- "#2c7fb8"

# ============ shared data ============
jjmap <- read_csv(file.path(ROOT,"fits/bayes_long/summaries/marchman_a3_m3_psi.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item))) |> select(jj, word, delta_j)
key <- read_csv(file.path(ROOT,"data/intermediates/cdi_master_item_key.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item))) |>
  transmute(word, lexical_class=lexical_category, logfreq=log10(prob)) |> distinct(word,.keep_all=TRUE)
lad <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"), show_col_types=FALSE) |>
  group_by(word, words) |> summarise(surp=mean(surprisal), .groups="drop") |> mutate(ln=log(words))

# ---------- Panel A ----------
b <- readRDS(file.path(ROOT,"fits/bayes_long/bundle_marchman_a3.rds")); sd <- b$stan_data
obs <- tibble(jj=sd$jj, y=sd$y, age=sd$admin_age[sd$aa])
cand <- obs |> group_by(jj) |> summarise(pmin=mean(y[age<=min(age)+3]), pmax=mean(y[age>=max(age)-3]), .groups="drop") |>
  inner_join(jjmap,by="jj") |> inner_join(key,by="word") |> filter(word %in% lad$word, pmax-pmin>0.35, pmax>0.6)
sel <- cand |> filter(lexical_class %in% c("nouns","predicates","function_words","other")) |>
  group_by(lexical_class) |> slice_max(pmax-pmin, n=1) |> ungroup()
r01 <- function(x)(x-min(x))/(max(x)-min(x)); rows<-list()
for(i in seq_len(nrow(sel))){ w<-sel$word[i]; jj<-sel$jj[i]
  ch <- obs |> filter(jj==!!jj) |> mutate(agb=round(age/2)*2) |> group_by(agb) |> summarise(p=mean(y),n=n(),.groups="drop") |>
    filter(n>=8) |> mutate(p=pmin(pmax(p,.02),.98), logit=log(p/(1-p)), ln=log(agb))
  kc<-coef(lm(logit~ln,ch))[2]; md<-lad |> filter(word==!!w) |> mutate(ev=-surp); km<-coef(lm(ev~ln,md))[2]
  rows[[i]] <- bind_rows(
    ch |> transmute(word=w, who="Children", x=r01(ln), evidence=logit),
    md |> transmute(word=w, who="LM", x=r01(ln), evidence=ev)) |>
    mutate(lab=sprintf("%s  (%s)\nkappa: child %.0f  vs  LM %.1f", w, sub("_"," ",sel$lexical_class[i]), kc, km)) }
dA <- bind_rows(rows)
pA <- ggplot(dA, aes(x, evidence, colour=who)) + geom_hline(yintercept=0, linewidth=.3, colour="grey85") +
  geom_point(size=1, alpha=.7) + geom_smooth(method="lm", se=FALSE, linewidth=.9) +
  facet_wrap(~lab, nrow=1, scales="free_y") +
  scale_colour_manual(values=c(Children=CHILD_RED, LM=LM_BLUE), name=NULL) +
  labs(subtitle="A. Per-word evidence accumulation: children (logit P produced, +) rise steeply; LMs (-surprisal, -) rise gradually",
       x="experience (rescaled log: child age | model data)", y="accumulated evidence") +
  theme_minimal(base_size=10) + theme(legend.position="right", strip.text=element_text(size=7.5),
    plot.subtitle=element_text(face="bold"))

# ---------- Panel B ----------
dB <- jjmap |> inner_join(lad |> group_by(word) |> summarise(final_surp=surp[which.max(ln)],.groups="drop"),by="word") |>
  inner_join(key,by="word") |> filter(lexical_class %in% c("nouns","predicates","function_words","other")) |>
  mutate(lexical_class=recode(lexical_class, function_words="function words"))
labB <- dB |> group_by(lexical_class) |> summarise(txt=sprintf("r=%.2f",cor(final_surp,delta_j)),.groups="drop")
palB <- c("nouns"="#1b9e77","predicates"="#d95f02","function words"="#7570b3","other"="#e7298a")
pB <- ggplot(dB, aes(delta_j, final_surp, colour=lexical_class)) +
  geom_point(alpha=.5,size=1) + geom_smooth(method="lm",se=FALSE,linewidth=.8) +
  geom_text(data=labB, aes(x=-Inf,y=Inf,label=txt), hjust=-0.1,vjust=1.4,size=3,colour="grey25",inherit.aes=FALSE) +
  facet_wrap(~lexical_class,nrow=1) + scale_colour_manual(values=palB,guide="none") +
  labs(subtitle="B. LM vs child word difficulty, by lexical class (positive = agree)",
       x="child difficulty  delta", y="LM difficulty (surprisal, 24M)") +
  theme_minimal(base_size=10) + theme(strip.text=element_text(size=8,face="bold"),
    plot.subtitle=element_text(face="bold"))

# ---------- Panel C: kappa density (current Fig 3) ----------
sl <- readRDS(file.path(ROOT,"paper/cache/fig6_llm_slopes.rds"))$slopes
pal <- c("Children (English)"="#c41e37","Children (Norwegian)"="#fb9a99",
  "LMs: C&B 2022 (4 architectures)"="#8856a7","LMs: CHILDES (training)"="#2c7fb8","LMs: CHILDES (development)"="#41ab5d")
pC <- ggplot(sl, aes(slope_natural, fill=group, colour=group)) + geom_density(alpha=.32, linewidth=.5) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey55", linewidth=.4) +
  annotate("text", x=1.4, y=0.9, hjust=0, size=2.8, colour="grey45", label="kappa=1\n(no acceleration)") +
  scale_fill_manual(values=pal,name=NULL) + scale_colour_manual(values=pal,name=NULL) +
  coord_cartesian(xlim=c(0,22)) +
  labs(subtitle="C. Distribution of acceleration (kappa): children by child, LMs by word", x="kappa", y="density") +
  theme_minimal(base_size=10) + theme(legend.position="right", legend.text=element_text(size=7),
    plot.subtitle=element_text(face="bold"))

comp <- pA / pB / pC + plot_layout(heights=c(1,1,1))
ggsave(file.path(OUT,"fig_composite.png"), comp, width=11, height=10.5, dpi=140)
cat("wrote fig_composite.png\n")
