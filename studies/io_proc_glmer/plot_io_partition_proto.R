## PROTOTYPE (interim, not for the paper): an io-proc COEFFICIENT panel toward a future
## fig-io-partition update, following the paper's fig-io-partition(B) plotting grammar:
## vertical continuous axis, efficiency vs acceleration as the categorical x split, the
## paper's channel palette (Input green / Processing pink), theme_minimal.
## Two estimates of the same thing:
##   - glmer "both" (joint input+processing, 142 kids w/ both channels, full longitudinal)
##   - io-proc-lean Bayesian (tau=1, 403 kids, ALL current data incl. one-timepoint kids)
## Efficiency (level) coefs as-is; acceleration coefs rescaled x log(t_ref/a0) so the two
## sit on ONE theta-units axis ("contribution to log-odds of production by 30 mo").
## Renders TWO variants: (1) facet-by-model, (2) single panel w/ model as point shape.
## RUN LOCALLY.  -> figs/io_partition_proto_facet.png, figs/io_partition_proto_panel.png
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here); library(readr)})

co  <- read_csv(here("studies/io_proc_glmer/results/glmer_ladder_coefs.csv"), show_col_types = FALSE)
sm2 <- read_csv(here("studies/io_proc_glmer/results/sm2_overlay.csv"),        show_col_types = FALSE)

A0 <- 21; T_REF <- 30; KACC <- log(T_REF / A0)   # acceleration -> level-equivalent theta units (0.357)
FACTOR_PAL <- c(Input = "#009E73", Processing = "#CC79A7")   # paper fig-io-partition channel palette
SUB <- sprintf("acceleration rescaled ×log(%d/%d)=%.2f to level-equivalent θ units; τ=1", T_REF, A0, KACC)
YLAB <- sprintf("contribution to log-odds of production by %d mo (θ units)", T_REF)

d <- bind_rows(filter(co, spec == "both_common"), filter(sm2, spec == "lean_t1")) %>% mutate(
  acc  = term == "acceleration",
  est2 = ifelse(acc, est * KACC, est),
  lo2  = ifelse(acc, lo  * KACC, lo),
  hi2  = ifelse(acc, hi  * KACC, hi),
  term_x  = factor(term, levels = c("level", "acceleration"),
                   labels = c("efficiency", paste0("acceleration\n(by ", T_REF, " mo)"))),
  channel = factor(channel, levels = c("input", "processing"), labels = c("Input", "Processing")),
  model   = factor(recode(spec, both_common = "glmer", lean_t1 = "io-proc-lean"),
                   levels = c("glmer", "io-proc-lean")),
  model_lab = factor(recode(spec, both_common = "glmer (both channels, n=142)",
                                  lean_t1     = "io-proc-lean (Bayesian, all data, n=403)"),
                     levels = c("glmer (both channels, n=142)",
                                "io-proc-lean (Bayesian, all data, n=403)")))

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold", size = 10.5),
        plot.margin = margin(5.5, 5.5, 5.5, 14),
        legend.position = "top")

## ---- Variant 1: facet by model ------------------------------------------------------
dg <- position_dodge(width = 0.5)
p_facet <- ggplot(d, aes(term_x, est2, color = channel, group = channel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_pointrange(aes(ymin = lo2, ymax = hi2), position = dg, fatten = 3, linewidth = 0.9) +
  facet_wrap(~ model_lab) +
  scale_color_manual(values = FACTOR_PAL, name = NULL) +
  labs(x = NULL, y = YLAB, title = "Input vs processing: efficiency & acceleration", subtitle = SUB) +
  base_theme
ggsave(here("studies/io_proc_glmer/figs/io_partition_proto_facet.png"), p_facet,
       width = 8.5, height = 4.0, dpi = 150)

## ---- Variant 2: single panel, model as point shape ----------------------------------
dg2 <- position_dodge(width = 0.6)
p_panel <- ggplot(d, aes(term_x, est2, color = channel, shape = model,
                         group = interaction(channel, model))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_linerange(aes(ymin = lo2, ymax = hi2), position = dg2, linewidth = 0.8) +
  geom_point(position = dg2, size = 3) +
  scale_color_manual(values = FACTOR_PAL, name = NULL) +
  scale_shape_manual(values = c(glmer = 16, `io-proc-lean` = 17), name = NULL) +
  labs(x = NULL, y = YLAB, title = "Input vs processing: efficiency & acceleration", subtitle = SUB) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  base_theme
ggsave(here("studies/io_proc_glmer/figs/io_partition_proto_panel.png"), p_panel,
       width = 7.5, height = 4.0, dpi = 150)

cat("wrote figs/io_partition_proto_{facet,panel}.png\n")
