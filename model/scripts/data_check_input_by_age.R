## Data-check: observed input (log rate) by age, faceted by study, to test whether
## input DRIFTS WITH AGE. This decides (a) whether the input measurement model
## needs an age slope, and (b) whether the within-child replicate spread that
## identifies sigma_meas is pure noise (flat) or contaminated by age drift.
## We annotate the WITHIN-CHILD age slope (child-demeaned) per study -- that's the
## sigma_meas-relevant quantity (between-child age differences don't contaminate it).
## RUN LOCALLY.  Output: figs/data_checks/input_by_age.png
suppressPackageStartupMessages({ library(here); library(dplyr); library(tidyr); library(ggplot2); library(readr) })

## ---- gather per-recording (study, child, age_mo, log_input) from each source ----
lena <- read_csv(here("data/raw/AM2018/lena_am2018_fmw2013.csv"), show_col_types = FALSE) %>%
  filter(Study == "TL3") %>%                           # AM2018 (TL3) only from old file
  transmute(study = "adams_marchman_2018", child = as.character(SubjectID1),
            a16 = AGE16M, r16 = AWCHr16M, a18 = AGE18M, r18 = AWCHr18M) %>%
  pivot_longer(c(a16, r16, a18, r18), names_to = c(".value", "tp"),
               names_pattern = "([ar])(16|18)") %>%
  transmute(study, child, age_mo = a, log_input = log(r)) %>% filter(is.finite(log_input), age_mo > 0)
## FMW2013 input now from the cleaned 2-timepoint file (TLO+ELENA, 18+24mo)
fmw <- read_csv(here("data/raw/FMW2013/TLOELENA_LENA_1824.csv"), show_col_types = FALSE) %>%
  transmute(study = "fmw_2013", child = as.character(SubjectID1),
            a18 = AGE18M, r18 = AWCHr18M, a24 = AGE24M, r24 = AWCHr24M) %>%
  pivot_longer(c(a18, r18, a24, r24), names_to = c(".value", "tp"),
               names_pattern = "([ar])(18|24)") %>%
  transmute(study, child, age_mo = a, log_input = log(r)) %>% filter(is.finite(log_input), age_mo > 0)

seed <- read_csv(here("data/raw/seedlings/lena_data.csv"), show_col_types = FALSE) %>%
  transmute(study = "SEEDLingS", child = as.character(subj), age_mo = month,
            log_input = log(awc_perhr)) %>%
  filter(is.finite(log_input), age_mo <= 30)   # drop the 4;6 (54mo) follow-up; model uses 6-17mo

bv <- readRDS(here("fits/babyview_subset_data.rds"))$videos %>%
  transmute(study = "BabyView", child = as.character(subject_id), age_mo, log_input = log_r_obs) %>%
  filter(is.finite(log_input))

dat <- bind_rows(lena, fmw, seed, bv) %>%
  mutate(study = factor(study, levels = c("adams_marchman_2018","fmw_2013","SEEDLingS","BabyView")))

## ---- within-child age slope (child-demeaned) per study ----
slopes <- dat %>% group_by(study, child) %>% filter(n() >= 2) %>%
  mutate(age_c = age_mo - mean(age_mo), log_c = log_input - mean(log_input)) %>%
  group_by(study) %>%
  summarise(b = coef(lm(log_c ~ age_c - 1))[1],
            se = summary(lm(log_c ~ age_c - 1))$coefficients[1, 2],
            n_kids = n_distinct(child), .groups = "drop") %>%
  mutate(lab = sprintf("within-child slope: %.3f ± %.3f log/mo\n(%s; %d kids w/ repeats)",
                       b, se, ifelse(abs(b) < 2 * se, "n.s.", "p<.05"), n_kids))
cat("within-child age slopes (log input units / month):\n"); print(as.data.frame(slopes[, c("study","b","se","n_kids")]))

## ---- plot ----
labpos <- dat %>% group_by(study) %>% summarise(x = min(age_mo), y = max(log_input), .groups = "drop") %>%
  left_join(slopes, by = "study")
p <- ggplot(dat, aes(age_mo, log_input)) +
  geom_line(aes(group = child), color = "grey75", alpha = 0.4, linewidth = 0.25) +
  geom_point(alpha = 0.35, size = 0.6, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#1f78b4", linewidth = 0.9) +
  geom_text(data = labpos, aes(x, y, label = lab), hjust = 0, vjust = 1, size = 2.6, color = "#1f78b4") +
  facet_wrap(~ study, scales = "fixed", nrow = 1) +   # shared x AND y limits for comparison
  labs(x = "age (months)", y = "log input rate (study-specific units)",
       title = "Observed input by age, per child (blue = pooled lm; text = within-child slope)",
       subtitle = "Flat within-child slope => replicate spread is ~pure noise => estimate sigma_meas freely; sloped => input model needs an age term") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())

dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figs", "data_checks", "input_by_age.png"), p, width = 13, height = 3.6, dpi = 140)
cat("wrote figs/data_checks/input_by_age.png\n")
