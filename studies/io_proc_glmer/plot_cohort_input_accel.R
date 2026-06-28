## Input->acceleration by DATA SOURCE: a model-independent (glmer) look at which cohorts can
## identify the channel at all. The English both-channel sample (142) is the only informative
## one; every NEW cohort -- clean Spanish item-level AND the sumscore cohorts -- has a wide CI
## spanning 0. So the bilingual additions can't contribute input->accel signal, only noise +
## sigma_zeta inflation (4.16->5.60), which is why pooling them destabilizes gamma_in.
## RUN LOCALLY. -> figs/cohort_input_accel.png
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here); library(readr)})

en  <- read_csv(here("studies/io_proc_glmer/results/glmer_ladder_coefs.csv"), show_col_types = FALSE) %>%
  filter(spec == "both_common", channel == "input", term == "acceleration") %>%
  transmute(spec = "English \"both\" (EN, item+input)", n_kids, est, lo, hi, kind = "English reference")
co <- read_csv(here("studies/io_proc_glmer/results/cohort_glmer_coefs.csv"), show_col_types = FALSE) %>%
  filter(term == "acceleration") %>% transmute(spec, n_kids, est, lo, hi, kind = "NEW cohort (added to bilingual)")

d <- bind_rows(en, co) %>%
  mutate(lab = sprintf("%s  (n=%d)", spec, n_kids),
         lab = factor(lab, levels = rev(lab)))

p <- ggplot(d, aes(est, lab, color = kind)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_vline(xintercept = en$est, linetype = "dashed", linewidth = 0.4, color = "#009E73") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.9) +
  geom_point(size = 3) +
  scale_color_manual(values = c("English reference" = "#009E73", "NEW cohort (added to bilingual)" = "#D55E00"), name = NULL) +
  labs(x = "input → acceleration coefficient (glmer, per 1 SD within-cohort input)", y = NULL,
       title = "Only the English both-channel sample identifies input→acceleration",
       subtitle = "every new cohort (incl. clean Spanish item-level) has a wide CI spanning 0 → pooling adds noise, not signal") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9))
ggsave(here("studies/io_proc_glmer/figs/cohort_input_accel.png"), p, width = 9, height = 3.8, dpi = 140)
cat("wrote figs/cohort_input_accel.png\n")
