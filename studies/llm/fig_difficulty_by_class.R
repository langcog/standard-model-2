#!/usr/bin/env Rscript
# Candidate composite PANEL B: LM word-difficulty vs children's difficulty (delta_j), broken down by
# lexical class (matches fig-efficiency's Dark2 classes). Tests where LM & children (dis)agree.
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2)})
ROOT <- "/Users/mcfrank/Projects/standard_model_2"; OUT <- file.path(ROOT,"studies/llm")

cd <- read_csv(file.path(ROOT,"fits/bayes_long/summaries/marchman_a3_m3_psi.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item))) |> filter(!grepl(" ",word)) |> select(word, delta_j)
lad <- read_csv(file.path(ROOT,"fits/llm/ladder_bestval_finer.csv"), show_col_types=FALSE) |> mutate(lx=log10(words))
lmw <- lad |> group_by(word,words,lx) |> summarise(s=mean(surprisal),.groups="drop") |>
  group_by(word) |> summarise(final_surp=s[which.max(lx)], .groups="drop")
key <- read_csv(file.path(ROOT,"data/intermediates/cdi_master_item_key.csv"), show_col_types=FALSE) |>
  mutate(word=tolower(sub("^(id|ul):","",item))) |>
  transmute(word, lexical_class=lexical_category, logfreq=log10(prob)) |> distinct(word,.keep_all=TRUE)
d <- cd |> inner_join(lmw,by="word") |> inner_join(key,by="word") |>
  filter(lexical_class %in% c("nouns","predicates","function_words","other")) |>
  mutate(lexical_class=recode(lexical_class, function_words="function words"))

# per-class correlations (LM-child, and LM-freq to show the mechanism)
lab <- d |> group_by(lexical_class) |> summarise(
  n=n(), r_lmchild=cor(final_surp,delta_j), r_lmfreq=cor(final_surp,logfreq),
  r_childfreq=cor(delta_j,logfreq), .groups="drop") |>
  mutate(txt=sprintf("r(LM,child)=%.2f  (n=%d)", r_lmchild, n))
cat("per-class: LM-child, LM-freq, child-freq correlations\n"); print(as.data.frame(lab[,1:5]), row.names=FALSE)

pal <- c("nouns"="#1b9e77","predicates"="#d95f02","function words"="#7570b3","other"="#e7298a")
p <- ggplot(d, aes(delta_j, final_surp)) +
  geom_point(aes(colour=lexical_class), alpha=.55, size=1.2) +
  geom_smooth(aes(colour=lexical_class), method="lm", se=FALSE, linewidth=.8) +
  geom_text(data=lab, aes(x=-Inf, y=Inf, label=txt), hjust=-0.05, vjust=1.4, size=3, colour="grey25") +
  facet_wrap(~lexical_class, nrow=1) +
  scale_colour_manual(values=pal, guide="none") +
  labs(title="LM vs child word difficulty, by lexical class",
       subtitle="each point a word; positive slope = LM & children agree on difficulty",
       x="child difficulty  δ  (bayes_long M3; higher = harder for children)",
       y="LM difficulty (surprisal at 24M words)") +
  theme_minimal(base_size=11) + theme(strip.text=element_text(face="bold"))
ggsave(file.path(OUT,"fig_difficulty_by_class.png"), p, width=11, height=3.6, dpi=140)
cat("wrote fig_difficulty_by_class.png\n")
