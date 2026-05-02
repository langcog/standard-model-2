## Forest plot of per-class log-frequency slope beta_c from
## class_beta_slopes fits on English (and Norwegian when available).
##
## beta_c = 1 means frequency contributes per the unit-accumulator
## baseline; beta_c < 1 means frequency contributes less; beta_c = 0
## means frequency is uninformative once per-word psi_j is free.

source("model/R/config.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(posterior); library(rstan)
})

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

extract_beta_c <- function(variant, dataset_key, lang_label) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) {
    cat(sprintf("  missing: %s\n", path)); return(NULL)
  }
  fit <- readRDS(path)
  d <- as_draws_df(fit)
  bundle <- load_dataset_bundle(dataset_key)
  C <- bundle$stan_data$C
  class_names <- bundle$class_levels[seq_len(C)]

  bind_rows(lapply(seq_len(C), function(c) {
    x <- d[[sprintf("beta_c[%d]", c)]]
    tibble(language = lang_label,
           class_id = c,
           class    = class_names[c],
           median = median(x),
           lo95   = quantile(x, 0.025, names = FALSE),
           hi95   = quantile(x, 0.975, names = FALSE),
           lo80   = quantile(x, 0.10, names = FALSE),
           hi80   = quantile(x, 0.90, names = FALSE))
  }))
}

cat("Loading class_beta fits...\n")
res <- bind_rows(
  extract_beta_c("long_class_beta_slopes",            "english",   "English"),
  extract_beta_c("long_class_beta_slopes_norwegian",  "norwegian", "Norwegian")
)

if (nrow(res) == 0) stop("No class_beta fits found.")

# Print
res_print <- res %>%
  mutate(cell = sprintf("%.2f [%.2f, %.2f]", median, lo95, hi95))
cat("\nbeta_c posterior summary:\n")
print(as.data.frame(res_print %>% select(language, class, cell)),
      row.names = FALSE)

# ---- Forest plot ---- #
res <- res %>%
  mutate(class = factor(class,
                        levels = c("function_words", "predicates",
                                    "other", "nouns")))

p <- ggplot(res, aes(y = class, color = language)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  geom_pointrange(aes(x = median, xmin = lo95, xmax = hi95),
                  position = position_dodge(width = 0.4),
                  size = 0.4, linewidth = 0.5) +
  geom_pointrange(aes(x = median, xmin = lo80, xmax = hi80),
                  position = position_dodge(width = 0.4),
                  size = 0.6, linewidth = 1.2, fatten = 0) +
  geom_text(aes(x = median,
                label = sprintf("%.2f", median)),
            position = position_dodge(width = 0.4),
            color = "black", size = 3, hjust = -0.6, vjust = -0.6) +
  scale_color_manual(values = c("English" = "#1f77b4",
                                 "Norwegian" = "#d62728")) +
  scale_x_continuous(breaks = seq(-0.2, 1.2, by = 0.2),
                     limits = c(-0.3, 1.2)) +
  labs(x = expression("Per-class log-frequency slope " ~ beta[c]),
       y = NULL,
       title = expression("Class-specific frequency slope " ~ beta[c]),
       subtitle = paste0("Dashed line at 1 = unit-accumulator baseline.  ",
                          "Dotted at 0 = frequency irrelevant.  ",
                          "Wider pointrange = 95% CrI; thicker = 80%.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_FIGS, "class_beta_forest.png"),
       p, width = 9, height = 5, dpi = 200)
cat("\nWrote class_beta_forest.png\n")
