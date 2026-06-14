## Data-check plot: per-child log-RT (LWL) trajectories by age, faceted by the
## RT datasets in the joint io+proc bundle. Companion to the CDI spaghetti.
## RUN LOCALLY.  Output: figs/data_checks/rt_spaghetti.png
suppressPackageStartupMessages({ library(here); library(dplyr); library(ggplot2) })
b  <- readRDS(here("fits", "joint_io_proc_subset_data.rds"))
nm <- c("AM2018", "FM2012", "FMW2013", "fernald_totlot", "BabyView", "SEEDLingS")
lwl <- b$lwl %>%
  left_join(b$child_info %>% transmute(ii, ds = nm[study]), by = "ii") %>%
  filter(!is.na(ds)) %>%
  transmute(ds, lwl_age, lwl_log_rt, gid = paste0(ds, "::", ii))

## SEEDLings derived RT (Zhu et al., the new channel) -- from the extractor output
se <- readr::read_csv(here("data/seedlings/seedlings_lwl_rt.csv"), show_col_types = FALSE) %>%
  transmute(ds = "SEEDLingS-LWL", lwl_age, lwl_log_rt, gid = paste0(ds, "::", lab_subject_id))
lwl <- bind_rows(lwl, se)
FACETS <- c("AM2018", "FM2012", "FMW2013", "fernald_totlot", "SEEDLingS-LWL")

p <- ggplot(lwl %>% filter(ds %in% FACETS), aes(lwl_age, lwl_log_rt)) +
  geom_line(aes(group = gid), color = "grey70", alpha = 0.4, linewidth = 0.3) +
  geom_point(color = "grey55", alpha = 0.4, size = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "#c41e37", linewidth = 0.9, span = 1) +
  facet_wrap(~ factor(ds, levels = FACETS), nrow = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = "RT (ms)", breaks = c(500, 1000, 2000))) +
  labs(x = "age (months)", y = "log RT",
       title = "LWL reaction time by age, per child (SEEDLingS-LWL = newly derived from Zhu et al. raw)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 9))

dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figs", "data_checks", "rt_spaghetti.png"), p, width = 13, height = 3.4, dpi = 140)
cat("wrote figs/data_checks/rt_spaghetti.png\n")
