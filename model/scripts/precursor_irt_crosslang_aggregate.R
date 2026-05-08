## Aggregate cached per-language precursor fits into the panel + table.
##
## Reads everything in fits/precursor/*.rds (whatever has finished so
## far in precursor_irt_crosslang.R) and regenerates:
##   outputs/figs/longitudinal/precursor_crosslang_summary.csv
##   outputs/figs/longitudinal/precursor_crosslang_sd_by_age.csv
##   outputs/figs/longitudinal/precursor_crosslang_panel.png
##   outputs/figs/longitudinal/precursor_crosslang_summary.png
##
## Useful for rebuilding the deck while the long-running per-language
## script is still working through additional languages.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

OUT_DIR   <- file.path(PATHS$figs_dir, "longitudinal")
CACHE_DIR <- file.path(PATHS$fits_dir, "precursor")
MIN_OBS_PER_AGE <- 5

cache_files <- list.files(CACHE_DIR, pattern = "\\.rds$", full.names = TRUE)
cat(sprintf("Found %d cached language fits in %s\n",
            length(cache_files), CACHE_DIR))

results <- lapply(cache_files, readRDS)
results <- Filter(Negate(is.null), results)
cat(sprintf("Loaded %d non-null fits.\n", length(results)))

# Drop any fit that doesn't have the new format
results <- Filter(function(r) all(c("language", "pop_slope",
                                    "sigma_theta_within",
                                    "age_summary") %in% names(r)),
                  results)
cat(sprintf("After format filter: %d fits.\n", length(results)))
if (length(results) == 0) stop("No usable fits.")

# Per-language scalars
summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(language = r$language,
             n_admins = r$n_admins, n_items = r$n_items,
             age_min = r$age_range[1], age_max = r$age_range[2],
             pop_slope = r$pop_slope,
             sigma_theta_pooled = r$sigma_theta_pooled,
             sigma_theta_within = r$sigma_theta_within,
             months_per_sd_within = r$months_per_sd_within,
             stringsAsFactors = FALSE)
})) |> arrange(desc(n_admins))
write.csv(summary_df,
          file.path(OUT_DIR, "precursor_crosslang_summary.csv"),
          row.names = FALSE)
cat(sprintf("Wrote: %s (%d languages)\n",
            file.path(OUT_DIR, "precursor_crosslang_summary.csv"),
            nrow(summary_df)))
print(summary_df, digits = 3, row.names = FALSE)

# Per-age-bin SD
sd_by_age <- do.call(rbind, lapply(results, function(r) {
  data.frame(language = r$language, age_int = r$age_summary$age_int,
             n = r$age_summary$n,
             mean_theta = r$age_summary$mean_theta,
             sd_theta = r$age_summary$sd_theta,
             stringsAsFactors = FALSE)
}))
write.csv(sd_by_age,
          file.path(OUT_DIR, "precursor_crosslang_sd_by_age.csv"),
          row.names = FALSE)

# Faceted panel
plot_df <- sd_by_age |>
  group_by(language) |>
  filter(n() >= 3) |>
  ungroup() |>
  mutate(language = factor(language,
                           levels = summary_df$language))

p_panel <- ggplot(plot_df, aes(age_int, sd_theta, group = language)) +
  geom_line(linewidth = 0.5, colour = "grey25", alpha = 0.7) +
  geom_point(size = 0.9, colour = "grey25") +
  facet_wrap(~ language,
             ncol = max(3, ceiling(sqrt(length(unique(plot_df$language)))))) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Age (months)",
       y = expression("Within-age SD of " * theta[admin] * " (logits)"),
       title = sprintf("Cross-language precursor: within-age SD of ability grows with age (%d languages)",
                       length(unique(plot_df$language))),
       subtitle = "1-PL Rasch IRT (glmer) on Wordbank cross-sectional WS data; per-admin theta extracted; SD computed per integer age bin.") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "precursor_crosslang_panel.png"),
       p_panel, width = 13, height = 8, dpi = 150)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_crosslang_panel.png")))

# Summary scatter
p_summary <- ggplot(summary_df,
                    aes(pop_slope, months_per_sd_within, label = language)) +
  geom_point(aes(size = n_admins), alpha = 0.7, colour = "steelblue") +
  geom_text(size = 2.6, hjust = -0.15, vjust = 0.4, check_overlap = TRUE) +
  scale_size_continuous(range = c(2, 7), name = "N admins") +
  labs(x = "Population age slope (logits / month)",
       y = "Within-age SD / pop slope (months / SD)",
       title = "Cross-language: 'months of development' equivalents",
       subtitle = "How many months of typical development +/- 1 SD across kids corresponds to.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "precursor_crosslang_summary.png"),
       p_summary, width = 10, height = 6.5, dpi = 150)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_crosslang_summary.png")))

cat("\nDone.\n")
