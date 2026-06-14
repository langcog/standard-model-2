## Data-check plot: per-child input-rate estimates by dataset and method
## (LENA vs head-cam), to see how comparable the rate estimates are across
## measurement methods. Feeds the sigma_r measurement-model design.
## RUN LOCALLY.  Output: figs/data_checks/input_rate_by_dataset.png
suppressPackageStartupMessages({ library(here); library(dplyr); library(ggplot2) })
iob <- readRDS(here("fits", "io_pooled_subset_data.rds")); sd <- iob$stan_data
ii2 <- iob$df %>% distinct(ii, study)
pc <- tibble(logr = sd$log_r_obs, ii = sd$recording_to_child) %>%
  left_join(ii2, by = "ii") %>%
  group_by(study, ii) %>% summarise(mlog = mean(logr), n_rec = n(), .groups = "drop") %>%
  mutate(method = ifelse(study == "BabyView", "head-cam", "LENA"),
         study  = factor(study, levels = c("AM2018", "FMW2013", "SEEDLingS", "BabyView")))

MPAL <- c("LENA" = "#0072B2", "head-cam" = "#D55E00")
p <- ggplot(pc, aes(study, mlog, color = method)) +
  geom_violin(aes(fill = method), alpha = 0.15, color = NA, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = NA, linewidth = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 0.9) +
  scale_color_manual(values = MPAL, name = NULL) + scale_fill_manual(values = MPAL, guide = "none") +
  labs(x = NULL, y = "per-child mean log input rate", title = "Input-rate estimates by dataset & method") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figs", "data_checks", "input_rate_by_dataset.png"), p, width = 7, height = 4, dpi = 140)
cat("wrote figs/data_checks/input_rate_by_dataset.png\n")
## headline: LENA and head-cam rate estimates overlap (medians ~7.0-7.4 log);
## FMW2013 has only 1 recording/kid (no within-child replicate for noise).
