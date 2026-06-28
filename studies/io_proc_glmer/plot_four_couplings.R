## All FOUR couplings as a 2x2 panel grid -- input/processing x efficiency(level)/acceleration --
## with every glmer + Bayesian estimate on a shared estimator y-axis. Common per-1-SD units;
## processing is in "faster-processing -> outcome" units on both sides (glmer proc_z = -residual;
## Bayesian eff_proc_* sign-flipped in sm2_overlay). Lets us read the whole dissociation at once.
## RUN LOCALLY. -> figs/four_couplings_forest.png
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here); library(readr)})

gl  <- read_csv(here("studies/io_proc_glmer/results/glmer_ladder_coefs.csv"),  show_col_types = FALSE)
coh <- read_csv(here("studies/io_proc_glmer/results/cohort_glmer_coefs.csv"),   show_col_types = FALSE)
slp <- read_csv(here("studies/io_proc_glmer/results/slena_proc.csv"),           show_col_types = FALSE)
sm2 <- read_csv(here("studies/io_proc_glmer/results/sm2_overlay.csv"),          show_col_types = FALSE)

coupling <- function(channel, term) ifelse(channel == "input",
  ifelse(term == "level", "input → efficiency", "input → acceleration"),
  ifelse(term == "level", "processing → efficiency", "processing → acceleration"))
row_of <- function(df, est_lab, method, language) df %>%
  transmute(coupling = coupling(channel, term), estimator = est_lab, method, language, est, lo, hi)

d <- bind_rows(
  row_of(filter(gl, spec == "both_common"),                 "glmer · EN both (142)",      "glmer", "English"),
  row_of(filter(gl, spec == "input_full"),                  "glmer · EN 1-channel",       "glmer", "English"),
  row_of(filter(gl, spec == "proc_full"),                   "glmer · EN 1-channel",       "glmer", "English"),
  row_of(filter(coh, grepl("SLENA", spec)),                 "glmer · SLENA (ES, 29)",     "glmer", "Spanish"),
  row_of(slp,                                               "glmer · SLENA (ES, 29)",     "glmer", "Spanish"),
  row_of(filter(coh, grepl("HABLA", spec)),                 "glmer · HABLA (ES, 103)",    "glmer", "Spanish"),
  row_of(filter(coh, grepl("ELENA", spec)),                 "glmer · ELENA (EN, 24)",     "glmer", "English"),
  row_of(filter(sm2, spec == "en_d2"),                      "Bayes · io-proc-lean EN (403)","Bayesian","English"),
  row_of(filter(sm2, spec == "enct_proc"),                  "Bayes · EN+count +proc (413)","Bayesian","English"),
  row_of(filter(sm2, spec == "enct_noproc"),                "Bayes · EN+count no-proc (413)","Bayesian","English"),
  row_of(filter(sm2, spec == "bi_d2"),                      "Bayes · bi-lean (558)",      "Bayesian","Bilingual"))

EST_ORDER <- rev(c("glmer · EN 1-channel", "glmer · EN both (142)", "glmer · ELENA (EN, 24)",
                   "glmer · HABLA (ES, 103)", "glmer · SLENA (ES, 29)",
                   "Bayes · io-proc-lean EN (403)", "Bayes · EN+count +proc (413)",
                   "Bayes · EN+count no-proc (413)", "Bayes · bi-lean (558)"))
d <- d %>% mutate(estimator = factor(estimator, EST_ORDER),
                  language  = factor(language, c("English","Spanish","Bilingual")),
                  method    = factor(method, c("glmer","Bayesian")),
                  coupling  = factor(coupling, c("input → efficiency","input → acceleration",
                                                 "processing → efficiency","processing → acceleration")))
PAL <- c(English="#009E73", Spanish="#D55E00", Bilingual="#7570b3")
p <- ggplot(d, aes(est, estimator, color=language, shape=method)) +
  geom_vline(xintercept=0, linewidth=0.3, color="grey55") +
  geom_linerange(aes(xmin=lo, xmax=hi), linewidth=0.7) +
  geom_point(size=2.6, fill="white", stroke=0.9) +
  facet_wrap(~coupling, scales="free_x", ncol=2) +
  scale_color_manual(values=PAL, name=NULL) +
  scale_shape_manual(values=c(glmer=16, Bayesian=24), name=NULL) +
  labs(x="coefficient (per 1 SD input / faster processing; log-odds: intercept for efficiency, log-age slope for acceleration)",
       y=NULL, title="All four couplings: the full input × processing dissociation, glmer + Bayesian",
       subtitle="input→accel credible + stable across ALL English fits (403/413, ±count, ±proc); only +Spanish (bi-lean) drops to 0.37. proc→efficiency strong everywhere") +
  theme_minimal(base_size=9.5) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        plot.title=element_text(face="bold", size=11), strip.text=element_text(face="bold", size=10),
        axis.text.y=element_text(size=8))
ggsave(here("studies/io_proc_glmer/figs/four_couplings_forest.png"), p, width=11, height=6.4, dpi=140)
cat("wrote figs/four_couplings_forest.png\n")
d %>% arrange(coupling, estimator) %>%
  transmute(coupling, estimator, est=round(est,2), ci=sprintf("[%.2f,%.2f]",lo,hi)) %>% as.data.frame() %>% print(row.names=FALSE)
