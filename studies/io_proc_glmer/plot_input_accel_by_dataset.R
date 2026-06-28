## input->acceleration per dataset, colored by language. All four ENGLISH datasets point
## positive (replicated, ~0.5-1.3, clustering near the pooled glmer ~0.85); SLENA Spanish is
## the lone negative -- suggestive of a possible language difference, but its CI is too wide to
## confirm (29 kids, input at one age). Motivates a language-varying gamma_in test.
## RUN LOCALLY. -> figs/input_accel_by_dataset.png
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here); library(readr)})

d <- read_csv(here("studies/io_proc_glmer/results/input_accel_by_dataset.csv"), show_col_types = FALSE) %>%
  mutate(lab = factor(sprintf("%s (%s, n=%d)", dataset, substr(language,1,2), n_kids)))
EN_POOL <- 0.845                                            # glmer both_common (English, 142) reference

d <- d %>% arrange(language, -est) %>% mutate(lab = factor(lab, levels = rev(lab)))
p <- ggplot(d, aes(est, lab, color = language)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_vline(xintercept = EN_POOL, linetype = "dashed", linewidth = 0.4, color = "#009E73") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.9) +
  geom_point(size = 3) +
  scale_color_manual(values = c(English = "#009E73", Spanish = "#D55E00"), name = NULL) +
  labs(x = "input → acceleration coefficient (glmer, per 1 SD within-dataset input)", y = NULL,
       title = "input→acceleration: replicated + positive across English datasets; Spanish is the lone negative",
       subtitle = "every English cohort points positive (0.5–1.3); SLENA −2.35 but CI [−6, 1.3] — suggestive, not conclusive") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold", size = 10.5), axis.text.y = element_text(size = 9))
ggsave(here("studies/io_proc_glmer/figs/input_accel_by_dataset.png"), p, width = 9, height = 3.8, dpi = 140)
cat("wrote figs/input_accel_by_dataset.png\n")
