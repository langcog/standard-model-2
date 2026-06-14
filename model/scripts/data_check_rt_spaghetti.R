## Data-check plot: per-child log-RT (LWL) trajectories by age, faceted by the
## RT datasets in the joint io+proc bundle. Companion to the CDI spaghetti.
## RUN LOCALLY.  Output: figs/data_checks/rt_spaghetti.png
suppressPackageStartupMessages({ library(here); library(dplyr); library(ggplot2) })
b  <- readRDS(here("fits", "joint_io_proc_subset_data.rds"))
nm <- c("AM2018", "FM2012", "FMW2013", "fernald_totlot", "BabyView", "SEEDLingS")
lwl <- b$lwl %>%
  left_join(b$child_info %>% transmute(ii, ds = nm[study]), by = "ii") %>%
  filter(!is.na(ds))

p <- ggplot(lwl, aes(lwl_age, lwl_log_rt)) +
  geom_line(aes(group = ii), color = "grey70", alpha = 0.4, linewidth = 0.3) +
  geom_point(color = "grey55", alpha = 0.4, size = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "#c41e37", linewidth = 0.9, span = 1) +
  facet_wrap(~ factor(ds, levels = nm[1:4]), nrow = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = "RT (ms)", breaks = c(500, 1000, 2000))) +
  labs(x = "age (months)", y = "log RT", title = "LWL reaction time by age, per child") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10))

dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figs", "data_checks", "rt_spaghetti.png"), p, width = 11, height = 3.4, dpi = 140)
cat("wrote figs/data_checks/rt_spaghetti.png\n")
