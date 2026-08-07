#!/usr/bin/env Rscript
# Composite PANEL A: per-word acquisition curves, children vs LM, BOTH RISING on a shared
# "log-evidence" axis. Child = logit P(produce) vs log-age (positive region, steep: kappa~10).
# Model = -surprisal vs log-data (negative region -> "models negative", shallow: kappa~1).
# Experience rescaled to [0,1] per learner so the SLOPES (kappa) are visually comparable;
# actual kappa annotated so the magnitude isn't lost.
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(tidyr)})
ROOT <- "/Users/mcfrank/Projects/standard_model_2"; OUT <- file.path(ROOT,"studies/llm")
CHILD_RED <- "#c41e37"; LM_BLUE <- "#2c7fb8"

# ---- child data: empirical P(produce word j) by age, from the Marchman bundle ----
b <- readRDS(file.path(ROOT,"fits/bayes_long/bundle_marchman_a3.rds")); sd <- b$stan_data
jjmap <- read_csv(file.path(ROOT,"fits/bayes_long/summaries/marchman_a3_m3_psi.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item))) |> select(jj, word, delta_j)
key <- read_csv(file.path(ROOT,"data/intermediates/cdi_master_item_key.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item))) |> transmute(word, lexical_class=lexical_category) |> distinct(word,.keep_all=TRUE)
obs <- tibble(jj=sd$jj, y=sd$y, age=sd$admin_age[sd$aa])          # per-observation age + produced

# ---- model data: -surprisal vs data budget (ladder, mean over seeds) ----
lad <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"), show_col_types=FALSE) |>
  group_by(word, words) |> summarise(surp=mean(surprisal), .groups="drop") |> mutate(ln=log(words))

# pick 4 words: one per class, acquired IN-window (child P rises through 0.5), present in ladder
cand <- obs |> group_by(jj) |> summarise(pmin=mean(y[age<=min(age)+3]), pmax=mean(y[age>=max(age)-3]), .groups="drop") |>
  inner_join(jjmap, by="jj") |> inner_join(key, by="word") |> filter(word %in% lad$word, pmax-pmin > 0.35, pmax>0.6)
set.seed(1)
sel <- cand |> group_by(lexical_class) |> slice_max(pmax-pmin, n=1) |> ungroup() |>
  filter(lexical_class %in% c("nouns","predicates","function_words","other"))

rescale01 <- function(x) (x-min(x))/(max(x)-min(x))
rows <- list()
for(i in seq_len(nrow(sel))){
  w <- sel$word[i]; jj <- sel$jj[i]
  # child: P(produce) by 2-month age bin -> logit; kappa_child = slope of logit on ln(age)
  ch <- obs |> filter(jj==!!jj) |> mutate(agb=round(age/2)*2) |> group_by(agb) |>
    summarise(p=mean(y), n=n(), .groups="drop") |> filter(n>=8) |>
    mutate(p=pmin(pmax(p,.02),.98), logit=log(p/(1-p)), ln=log(agb))
  kc <- coef(lm(logit~ln, ch))[2]
  # model: -surprisal vs ln(data); kappa_model = slope
  md <- lad |> filter(word==!!w) |> mutate(ev=-surp)
  km <- coef(lm(ev~ln, md))[2]
  rows[[i]] <- bind_rows(
    ch |> transmute(word=w, lexical_class=sel$lexical_class[i], who="Children", x=rescale01(ln), evidence=logit, kap=kc),
    md |> transmute(word=w, lexical_class=sel$lexical_class[i], who="LM", x=rescale01(ln), evidence=ev, kap=km))
}
dd <- bind_rows(rows)
# per-word facet label with both kappas
labs <- dd |> distinct(word, lexical_class, who, kap) |> pivot_wider(names_from=who, values_from=kap) |>
  mutate(lab=sprintf("%s (%s):  child kappa=%.0f vs LM kappa=%.1f", word, sub("_"," ",lexical_class), Children, LM))
dd <- dd |> left_join(labs |> select(word, lab), by="word")

p <- ggplot(dd, aes(x, evidence, colour=who)) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey80") +
  geom_point(size=1.1, alpha=.7) + geom_smooth(method="lm", se=FALSE, linewidth=.9) +
  facet_wrap(~lab, nrow=1, scales="free_y") +
  scale_colour_manual(values=c(Children=CHILD_RED, LM=LM_BLUE), name=NULL) +
  labs(title="How children vs language models accumulate evidence for individual words",
       subtitle="both rising = knowledge accumulating; children (logit P produced, +) rise steeply, LMs (-surprisal, -) rise gradually",
       x="experience (rescaled log: child age | model data)", y="accumulated evidence") +
  theme_minimal(base_size=10) + theme(legend.position="top", strip.text=element_text(size=8))
ggsave(file.path(OUT,"fig_panelA_curves.png"), p, width=11, height=3.6, dpi=140)
cat("selected words:\n"); print(as.data.frame(labs[,c("word","lexical_class","Children","LM")]), row.names=FALSE)
cat("wrote fig_panelA_curves.png\n")
