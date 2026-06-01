## Quick visual test of hypothesis 2 (input-channel age confound).
##
## If observed log_r covaries with RECORDING age within a study, the
## pooled IO model will attribute that age-related variation to the
## intercept (xi via log_r_dev) instead of the slope (slope_i),
## systematically deflating kappa_study. Scatter log_r_obs vs age per
## recording, faceted by study, with within-study Pearson correlation.
##
## Output: figs/io/input_by_age.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

bv <- readRDS(file.path(PATHS$fits_dir, "babyview_subset_data.rds"))$videos |>
  transmute(study = "BabyView", child = as.integer(child_ii),
            age = age_mo, log_r_obs)

ss <- readRDS(file.path(PATHS$fits_dir, "seedlings_subset_data.rds"))$recordings |>
  transmute(study = "SEEDLingS", child = as.integer(child_ii),
            age = as.numeric(month), log_r_obs)

parse_which <- function(w) as.numeric(sub("AWCHr(\\d+)M", "\\1", w))
am <- readRDS(file.path(PATHS$fits_dir, "io_am2018_subset_data.rds"))$recordings |>
  transmute(study = "AM2018", child = as.integer(child_ii),
            age = parse_which(which), log_r_obs)
fm <- readRDS(file.path(PATHS$fits_dir, "io_fmw2013_subset_data.rds"))$recordings |>
  transmute(study = "FMW2013", child = as.integer(child_ii),
            age = parse_which(which), log_r_obs)

rec <- bind_rows(bv, ss, am, fm) |>
  mutate(study = factor(study, levels = c("BabyView","SEEDLingS","AM2018","FMW2013")))

# Within-study correlations
corrs <- rec |>
  group_by(study) |>
  summarise(r = cor(age, log_r_obs),
            n_rec = n(),
            n_kids = n_distinct(child),
            age_range = sprintf("%.0f-%.0f", min(age), max(age)),
            .groups = "drop")
cat("\n=== Within-study age vs log_r_obs correlations ===\n")
print(corrs, digits = 3)

# Kid-mean correlations (more directly tied to the model bias mechanism)
kid_means <- rec |>
  group_by(study, child) |>
  summarise(age_kid = mean(age), log_r_kid = mean(log_r_obs), .groups = "drop")
kid_corrs <- kid_means |>
  group_by(study) |>
  summarise(r_kid = cor(age_kid, log_r_kid), n = n(), .groups = "drop")
cat("\n=== Kid-mean age vs kid-mean log_r_obs (between-kid only) ===\n")
print(kid_corrs, digits = 3)

# Plot
study_pal <- c(BabyView = "#E69F00", SEEDLingS = "#56B4E9",
               AM2018   = "#009E73", FMW2013   = "#D55E00")

ann <- corrs |> mutate(label = sprintf("r = %+.2f\n%d rec, %d kids",
                                       r, n_rec, n_kids))

p <- ggplot(rec, aes(age, log_r_obs)) +
  geom_point(aes(color = study), alpha = 0.45, size = 1.1) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              linewidth = 0.55, alpha = 0.18, formula = y ~ x) +
  geom_text(data = ann, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.08, vjust = 1.25,
            inherit.aes = FALSE, size = 3.2, lineheight = 0.95) +
  facet_wrap(~ study, scales = "free", nrow = 1) +
  scale_color_manual(values = study_pal) +
  labs(x = "age at recording (months)",
       y = "observed log input rate (log_r_obs)",
       title = "Within-study age vs observed input",
       subtitle = "If r is substantially > 0, age-driven variance leaks into xi (deflating slope_i)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

out <- "figs/io/input_by_age.png"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 13, height = 4.2, dpi = 150)
cat(sprintf("\nWrote %s\n", out))
