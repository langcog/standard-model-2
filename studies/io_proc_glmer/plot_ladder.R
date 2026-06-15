## Figure: input & processing coefficients (level + acceleration) with 95% CI
## across glmer specs/samples -- the robustness/benchmark view (cf. paper Fig 3B).
## v1 = glmer specs only; v2 = + Bayesian SM2 D'3 overlay (input channel, approx
## theta~logit scale). Reads results/glmer_ladder_coefs.csv. RUN LOCALLY.
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here); library(readr)})
co <- read_csv(here("studies/io_proc_glmer/results/glmer_ladder_coefs.csv"), show_col_types = FALSE)

lab <- c(input_full="input · full (193)", input_common="input · common (142)",
         both_common="both (adjusted) · common (142)",
         proc_common="proc · common (142)", proc_full="proc · full (326)",
         SM2_bayes="SM2 D'3 · Bayesian")
ord <- rev(c("input · full (193)","input · common (142)","both (adjusted) · common (142)",
             "proc · common (142)","proc · full (326)","SM2 D'3 · Bayesian"))
prep <- function(d) d %>% mutate(
  spec_lab = factor(lab[spec], levels = ord),
  term = factor(term, levels = c("level","acceleration"),
                labels = c("Level (intercept @ 21mo)","Acceleration (x log age)")),
  channel = factor(channel, levels = c("input","processing")))

dodge <- position_dodge(width = 0.55)
make_fig <- function(d, subtitle) ggplot(d, aes(est, spec_lab, color = channel, shape = model_type, group = channel)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_linerange(aes(xmin = lo, xmax = hi), position = dodge, linewidth = 0.7) +
  geom_point(position = dodge, size = 2.7) +
  facet_wrap(~ term, scales = "free_x") +
  scale_color_manual(values = c(input = "#1f78b4", processing = "#e6701b")) +
  scale_shape_manual(values = c(glmer = 16, `Bayesian (SM2 D3)` = 17)) +
  labs(x = "coefficient (log-odds per 1 SD of channel)", y = NULL, color = NULL, shape = NULL,
       title = "Input vs processing: level & acceleration coefficients", subtitle = subtitle) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10.5), strip.text = element_text(face = "bold"))

p1 <- make_fig(prep(filter(co, model_type == "glmer")),
               "Double dissociation: input -> ACCELERATION, processing -> LEVEL (stable across samples & adjustment)")
ggsave(here("studies/io_proc_glmer/figs/io_proc_glmer_coefs.png"), p1, width = 10, height = 3.4, dpi = 150)

p2 <- make_fig(prep(co),
               "Benchmark vs Bayesian: SM2 (triangles, input only) under-credits acceleration, over-credits level (approx theta~logit scale)")
ggsave(here("studies/io_proc_glmer/figs/io_proc_glmer_coefs_vs_sm2.png"), p2, width = 10, height = 3.8, dpi = 150)
cat("wrote figs/io_proc_glmer_coefs.png and figs/io_proc_glmer_coefs_vs_sm2.png\n")
