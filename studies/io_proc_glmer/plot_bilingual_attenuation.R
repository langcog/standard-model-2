## Bilingual diagnostic: the input->acceleration coefficient across three estimators of the
## SAME quantity, ordered by how many weak-signal kids each includes:
##   glmer "both"      (142 kids: input AND longitudinal item-level CDI -- the clean identifiers)
##   io-proc-lean EN   (403 kids, CURRENT Rasch-lean + wide-delta production fit, en_d2)
##   bilingual bi-lean (558 kids: + count-only HABLA/ELENA kids via the sumscore branch)
## English input->accel is CREDIBLE (0.62, CI clears 0); the bilingual count cohorts both
## halve it (0.62->0.37) and collapse convergence (ess 134->22) via the structural kappa
## feedback. proc->efficiency stays put (the dense RT channel). Apples-to-apples = same model.
## RUN LOCALLY. -> figs/bilingual_input_accel_attenuation.png
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here); library(readr)})

co  <- read_csv(here("studies/io_proc_glmer/results/glmer_ladder_coefs.csv"), show_col_types = FALSE)
sm2 <- read_csv(here("studies/io_proc_glmer/results/sm2_overlay.csv"),        show_col_types = FALSE)
KACC <- log(30 / 21)                                       # acceleration -> level-equivalent theta units
PAL  <- c(Input = "#009E73", Processing = "#CC79A7")

d <- bind_rows(filter(co, spec == "both_common"),
               filter(sm2, spec %in% c("en_d2", "bi_d2"))) %>%
  mutate(acc  = term == "acceleration",
         est2 = ifelse(acc, est * KACC, est),
         lo2  = ifelse(acc, lo  * KACC, lo),
         hi2  = ifelse(acc, hi  * KACC, hi),
         channel = factor(channel, levels = c("input", "processing"), labels = c("Input", "Processing")),
         term_x  = factor(term, levels = c("level", "acceleration"),
                          labels = c("efficiency (ξ)", "acceleration (κ, by 30mo)")),
         model = factor(recode(spec, both_common = "glmer\n(142: input+long.CDI)",
                                     en_d2       = "io-proc-lean EN\n(403, credible)",
                                     bi_d2       = "bilingual\n(558: +count-only)"),
                        levels = c("glmer\n(142: input+long.CDI)",
                                   "io-proc-lean EN\n(403, credible)",
                                   "bilingual\n(558: +count-only)")))

dg <- position_dodge(width = 0.55)
p <- ggplot(d, aes(model, est2, color = channel, group = channel)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_line(position = dg, linewidth = 0.5, alpha = 0.5) +
  geom_linerange(aes(ymin = lo2, ymax = hi2), position = dg, linewidth = 0.8) +
  geom_point(position = dg, size = 2.6) +
  facet_wrap(~ term_x) +
  scale_color_manual(values = PAL, name = NULL) +
  labs(x = NULL, y = "contribution to log-odds of production (θ units)",
       title = "input→acceleration attenuates + destabilizes as weak-signal kids are added",
       subtitle = "proc→efficiency is stable across all three; only the κ-input channel dilutes (acceleration ×log(30/21))") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(size = 8), strip.text = element_text(face = "bold"))
ggsave(here("studies/io_proc_glmer/figs/bilingual_input_accel_attenuation.png"), p, width = 9, height = 4.6, dpi = 140)
cat("wrote figs/bilingual_input_accel_attenuation.png\n")
d %>% filter(acc) %>% transmute(model = gsub("\n.*", "", model), channel, est2 = round(est2, 3),
                                ci = sprintf("[%.2f, %.2f]", lo2, hi2)) %>% as.data.frame() %>% print(row.names = FALSE)
