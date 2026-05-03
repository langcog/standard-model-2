## Targeted comparisons across the M0..M5 nested family on English.
##
## Mapping:
##   M0  long_m0                     no time, no freq, no slopes, no 2PL
##   M1  long_m1                     unit time + unit freq, delta pinned
##   M2  long_baseline               + free delta
##   M3  long_slopes                 + per-child zeta
##   M4  long_class_beta_slopes      + class-specific beta_c
##   M5  long_m5                     + 2PL (lambda_j)
##
## Outputs:
##   outputs/figs/longitudinal/nested_family_scalars.png
##   outputs/figs/longitudinal/nested_family_loo.csv
##
## Targeted comparisons:
##   M1 vs M0:  does adding the unit time + freq baseline help?
##   M2 vs M1:  is acceleration (delta) needed?              (RQ1)
##   M3 vs M2:  do per-child slopes (zeta_i) matter?
##   M4 vs M3:  does class-specific beta_c improve fit?      (RQ2)
##   M5 vs M4:  does per-word 2PL discrimination help?

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan); library(loo)
})

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

FAMILY <- list(
  list(label = "M0", variant = "long_m0",
       desc  = "no time, no freq"),
  list(label = "M1", variant = "long_m1",
       desc  = "unit time + freq, delta pinned"),
  list(label = "M2", variant = "long_baseline",
       desc  = "+ free delta"),
  list(label = "M3", variant = "long_slopes",
       desc  = "+ per-child slopes"),
  list(label = "M4", variant = "long_class_beta_slopes",
       desc  = "+ class-specific beta_c"),
  list(label = "M5", variant = "long_m5",
       desc  = "+ 2PL (lambda_j)")
)

PARAMS <- c("sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
            "pi_alpha", "delta", "sigma_lambda")

# ---- 1. Pull scalar posteriors ---- #
all_scalars <- list()
for (m in FAMILY) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", m$variant))
  if (!file.exists(path)) {
    cat(sprintf("MISSING: %s\n", path)); next
  }
  fit <- readRDS(path)
  d <- as_draws_df(fit)

  rows <- lapply(PARAMS, function(p) {
    if (!p %in% names(d)) return(NULL)
    x <- d[[p]]
    if (sd(x) < 1e-3) return(NULL)
    tibble(label = m$label, variant = m$variant, desc = m$desc,
           param = p, median = median(x),
           lo95 = quantile(x, 0.025, names = FALSE),
           hi95 = quantile(x, 0.975, names = FALSE))
  })
  # Also pull beta_c values per class if present
  if ("beta_c[1]" %in% names(d)) {
    bundle <- load_dataset_bundle("english")
    classes <- bundle$class_levels
    for (c in seq_along(classes)) {
      x <- d[[sprintf("beta_c[%d]", c)]]
      if (sd(x) < 1e-3) next
      rows <- c(rows, list(tibble(
        label = m$label, variant = m$variant, desc = m$desc,
        param = sprintf("beta_c[%s]", classes[c]),
        median = median(x),
        lo95 = quantile(x, 0.025, names = FALSE),
        hi95 = quantile(x, 0.975, names = FALSE))))
    }
  }
  all_scalars[[m$label]] <- bind_rows(rows)
}
scalars <- bind_rows(all_scalars) %>%
  mutate(label = factor(label, levels = sapply(FAMILY, function(m) m$label)))

cat("\nScalar posteriors across nested family:\n")
print(as.data.frame(scalars %>%
                    mutate(cell = sprintf("%.2f [%.2f, %.2f]",
                                           median, lo95, hi95)) %>%
                    select(label, param, cell) %>%
                    pivot_wider(names_from = param, values_from = cell)),
      row.names = FALSE, max = 200)

# ---- 2. Forest plot of key scalars across the chain ---- #
KEY_PARAMS <- c("sigma_alpha", "pi_alpha", "delta", "sigma_zeta",
                "sigma_lambda")
sf <- scalars %>% filter(param %in% KEY_PARAMS) %>%
  mutate(param = factor(param, levels = KEY_PARAMS))

p_forest <- ggplot(sf, aes(x = median, y = label,
                           xmin = lo95, xmax = hi95)) +
  geom_pointrange(color = "steelblue", size = 0.4) +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  scale_y_discrete(limits = rev) +
  labs(x = "posterior median, 95% CrI", y = NULL,
       title = "Nested family M0..M5: key scalar posteriors (English)",
       subtitle = "Forest by model. Empty cells = parameter not in that model.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_FIGS, "nested_family_scalars.png"),
       p_forest, width = 13, height = 4.5, dpi = 200)
cat(sprintf("\nWrote %s\n",
            file.path(OUT_FIGS, "nested_family_scalars.png")))

# ---- 3. LOO comparison (skip if any fit missing) ---- #
# This needs log_lik in the Stan model, which we don't currently emit;
# placeholder: report posterior log density at each draw for a rough
# fit-quality measure. Real LOO requires log_lik per observation.
cat("\nNote: LOO/WAIC requires per-obs log_lik in Stan generated quantities,\n")
cat("which the current models don't emit. Adding that is a small change\n")
cat("to log_irt_long.stan and worth doing for the formal comparison.\n")
cat("For now this script reports scalar posteriors only.\n")
