## ONE axis for everything: all glmer + Bayesian estimates of the two dissociation channels,
## input->acceleration and processing->efficiency, in common per-1-SD units (the Bayesian
## eff_input_k = gamma_in*sigma_r and eff_proc_xi are on the same scale as the glmer input_z /
## proc_z coefficients -- the io_partition prototype established this overlay is exact).
## Lets us compare carefully: where the clean English signal sits, where pooling lands it, and
## how noisy each new cohort is. RUN LOCALLY. -> figs/all_estimates_forest.png
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here); library(readr); library(tidyr)})

gl  <- read_csv(here("studies/io_proc_glmer/results/glmer_ladder_coefs.csv"),  show_col_types = FALSE)
coh <- read_csv(here("studies/io_proc_glmer/results/cohort_glmer_coefs.csv"),   show_col_types = FALSE)
ds  <- read_csv(here("studies/io_proc_glmer/results/input_accel_by_dataset.csv"),show_col_types = FALSE)
sm2 <- read_csv(here("studies/io_proc_glmer/results/sm2_overlay.csv"),          show_col_types = FALSE)

mk <- function(channel, label, method, language, est, lo, hi)
  tibble(channel, label, method, language, est, lo, hi)

## ---------- input -> acceleration ----------
IA <- bind_rows(
  # glmer pooled English (the clean reference)
  gl  %>% filter(spec=="input_full", term=="acceleration") %>% transmute(channel="input → acceleration",
          label="glmer · pooled EN (input+CDI, 193)", method="glmer", language="English", est, lo, hi),
  gl  %>% filter(spec=="both_common", term=="acceleration", channel=="input") %>% transmute(channel="input → acceleration",
          label="glmer · pooled EN both,142", method="glmer", language="English", est, lo, hi),
  # glmer per English dataset + SLENA
  ds  %>% transmute(channel="input → acceleration",
          label=sprintf("glmer · %s (%s, %d)", dataset, substr(language,1,2), n_kids),
          method="glmer", language, est, lo, hi),
  # glmer new sumscore cohorts
  coh %>% filter(term=="acceleration", grepl("HABLA|ELENA", spec)) %>% transmute(channel="input → acceleration",
          label=sprintf("glmer · %s", sub(" \\(.*","",spec)), method="glmer",
          language=ifelse(grepl("ES",spec),"Spanish","English"), est, lo, hi),
  # Bayesian pooled fits
  sm2 %>% filter(spec=="en_d2", channel=="input", term=="acceleration") %>% transmute(channel="input → acceleration",
          label="Bayesian · io-proc-lean EN (403)", method="Bayesian", language="English", est, lo, hi),
  sm2 %>% filter(spec=="bi_d2", channel=="input", term=="acceleration") %>% transmute(channel="input → acceleration",
          label="Bayesian · bi-lean (558)", method="Bayesian", language="Bilingual", est, lo, hi))

## ---------- processing -> efficiency (level) ----------
slp <- read_csv(here("studies/io_proc_glmer/results/slena_proc.csv"), show_col_types = FALSE)
PE <- bind_rows(
  gl  %>% filter(spec=="proc_full", term=="level") %>% transmute(channel="processing → efficiency",
          label="glmer · pooled EN (proc+CDI, 326)", method="glmer", language="English", est, lo, hi),
  gl  %>% filter(spec=="both_common", term=="level", channel=="processing") %>% transmute(channel="processing → efficiency",
          label="glmer · pooled EN (both, 142)", method="glmer", language="English", est, lo, hi),
  slp %>% filter(term=="level") %>% transmute(channel="processing → efficiency",
          label=sprintf("glmer · SLENA (ES, %d)", n_kids), method="glmer", language="Spanish", est, lo, hi),
  sm2 %>% filter(spec=="en_d2", channel=="processing", term=="level") %>% transmute(channel="processing → efficiency",
          label="Bayesian · io-proc-lean EN (403)", method="Bayesian", language="English", est, lo, hi),
  sm2 %>% filter(spec=="bi_d2", channel=="processing", term=="level") %>% transmute(channel="processing → efficiency",
          label="Bayesian · bi-lean (558)", method="Bayesian", language="Bilingual", est, lo, hi))

d <- bind_rows(IA, PE) %>%
  mutate(language = factor(language, c("English","Spanish","Bilingual")),
         method   = factor(method, c("glmer","Bayesian")),
         channel  = factor(channel, c("input → acceleration","processing → efficiency")))
# order rows within each facet by estimate (desc) so the forest reads cleanly
ord <- d %>% group_by(channel) %>% mutate(r = rank(est, ties.method="first")) %>% ungroup()
d$label <- factor(d$label, levels = unique(ord$label[order(ord$channel, ord$r)]))

PAL <- c(English="#009E73", Spanish="#D55E00", Bilingual="#7570b3")
p <- ggplot(d, aes(est, label, color=language, shape=method)) +
  geom_vline(xintercept=0, linewidth=0.3, color="grey55") +
  geom_linerange(aes(xmin=lo, xmax=hi), linewidth=0.7) +
  geom_point(size=2.7, fill="white", stroke=0.9) +
  facet_wrap(~channel, scales="free", ncol=1) +
  scale_color_manual(values=PAL, name=NULL) +
  scale_shape_manual(values=c(glmer=16, Bayesian=24), name=NULL) +
  labs(x="coefficient (per 1 SD input / processing; log-odds slope for accel, intercept for efficiency)", y=NULL,
       title="All estimates on one axis: input→acceleration vs processing→efficiency",
       subtitle="glmer (circles) + Bayesian (triangles); colour = language. Proc→efficiency agrees everywhere; input→accel is English-driven and noisy elsewhere") +
  theme_minimal(base_size=10) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        plot.title=element_text(face="bold", size=11), strip.text=element_text(face="bold", size=10.5),
        axis.text.y=element_text(size=8.5))
ggsave(here("studies/io_proc_glmer/figs/all_estimates_forest.png"), p, width=9.5, height=7.2, dpi=140)
cat("wrote figs/all_estimates_forest.png\n")
d %>% transmute(channel, label, method, language, est=round(est,2), ci=sprintf("[%.2f,%.2f]",lo,hi)) %>%
  as.data.frame() %>% print(row.names=FALSE)
