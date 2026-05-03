## Compare per-child ability trajectories for the input-observed
## datasets (BabyView, SEEDLingS), decomposing theta_it into its
## input-driven and efficiency-driven parts.
##
## In the io model, log_r_true_i and log_alpha_i are separate per-
## child latents (input is observed; the model can disambiguate them).
## So we can split:
##   theta_it = [ log_r_true_i + log H + log_age ]   ("input part")
##            + [ log_alpha_i + (delta + zeta_i) * log_age ]  ("efficiency part")
##
## Output: outputs/figs/io/io_ability_decomposition.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

IO_DATASETS <- list(
  babyview  = list(fit = "io_slopes.rds",           label = "BabyView (N=20)"),
  seedlings = list(fit = "io_slopes_seedlings.rds", label = "SEEDLingS (N=44)")
)

OUT_FIGS <- file.path(PATHS$figs_dir, "io")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)
AGE_GRID <- seq(8, 32, by = 0.25)

extract_decomp <- function(key) {
  ds  <- IO_DATASETS[[key]]
  bundle <- load_dataset_bundle(key)
  sd_    <- bundle$stan_data
  fit    <- readRDS(file.path(PATHS$fits_dir, ds$fit))
  draws  <- as_draws_df(fit)

  med <- function(p) median(draws[[p]])
  s     <- med("s")
  delta <- med("delta")
  log_H <- sd_$log_H
  a0    <- sd_$a0

  # Per-child latents: log_r_true and log_alpha are separate in io model
  log_r_true <- sapply(seq_len(sd_$I),
                        function(i) median(draws[[sprintf("log_r_true[%d]", i)]]))
  log_alpha  <- sapply(seq_len(sd_$I),
                        function(i) median(draws[[sprintf("log_alpha[%d]", i)]]))
  zeta       <- sapply(seq_len(sd_$I),
                        function(i) median(draws[[sprintf("zeta[%d]", i)]]))

  # Item difficulty (psi_j - log p_j) for the reference band
  psi  <- sapply(seq_len(sd_$J),
                 function(j) median(draws[[sprintf("psi[%d]", j)]]))
  log_p <- log(bundle$word_info$prob)
  beta  <- psi - log_p

  # theta on age grid for each kid, decomposed
  curves <- expand.grid(ii = seq_len(sd_$I), age = AGE_GRID,
                        KEEP.OUT.ATTRS = FALSE)
  curves <- curves %>% mutate(
    ae = pmax(age - s, 0.01),
    log_age = log(ae / a0),
    input_part      = log_r_true[ii] + log_H + log_age,
    efficiency_part = log_alpha[ii] + (delta + zeta[ii]) * log_age,
    theta           = input_part + efficiency_part
  )

  list(curves = curves %>% mutate(dataset = key, label = ds$label),
       beta   = tibble(beta = beta, dataset = key, label = ds$label))
}

cat("Extracting decompositions...\n")
out <- lapply(names(IO_DATASETS), extract_decomp)
curves <- bind_rows(lapply(out, `[[`, "curves")) %>%
  mutate(label = factor(label, levels = sapply(IO_DATASETS, `[[`, "label")))
beta_df <- bind_rows(lapply(out, `[[`, "beta")) %>%
  mutate(label = factor(label, levels = sapply(IO_DATASETS, `[[`, "label")))

# Pop mean per (dataset, age, component)
to_long <- function(df) {
  df %>%
    pivot_longer(c(input_part, efficiency_part, theta),
                 names_to = "component", values_to = "value") %>%
    mutate(component = recode(component,
                              input_part      = "Input contribution\n(log r + log H + log age)",
                              efficiency_part = "Efficiency contribution\n(log alpha + (delta+zeta) log age)",
                              theta           = "Total ability theta_it"),
           component = factor(component,
                              levels = c("Input contribution\n(log r + log H + log age)",
                                         "Efficiency contribution\n(log alpha + (delta+zeta) log age)",
                                         "Total ability theta_it")))
}

curves_long <- to_long(curves)

pop_mean <- curves_long %>%
  group_by(label, component, age) %>%
  summarise(value_mean = mean(value),
            lo10 = quantile(value, 0.10),
            hi90 = quantile(value, 0.90),
            .groups = "drop")

# Item difficulty band (for the total-theta panel only)
beta_summary <- beta_df %>%
  group_by(label) %>%
  summarise(beta_p10 = quantile(beta, 0.10),
            beta_p50 = quantile(beta, 0.50),
            beta_p90 = quantile(beta, 0.90),
            .groups = "drop") %>%
  mutate(component = factor("Total ability theta_it",
                            levels = levels(curves_long$component)))

p <- ggplot() +
  geom_rect(data = beta_summary,
            aes(ymin = beta_p10, ymax = beta_p90,
                xmin = -Inf, xmax = Inf),
            fill = "gray85", alpha = 0.4) +
  geom_hline(data = beta_summary,
             aes(yintercept = beta_p50),
             color = "gray50", linewidth = 0.4, linetype = "dashed") +
  geom_line(data = curves_long,
            aes(x = age, y = value, group = ii),
            color = "steelblue", alpha = 0.18, linewidth = 0.3) +
  geom_ribbon(data = pop_mean,
              aes(x = age, ymin = lo10, ymax = hi90),
              fill = "steelblue", alpha = 0.20) +
  geom_line(data = pop_mean,
            aes(x = age, y = value_mean),
            color = "steelblue", linewidth = 0.9) +
  facet_grid(component ~ label, scales = "free_y") +
  labs(x = "age (months)",
       y = "logit-scale contribution",
       title = "IO ability decomposition: input vs efficiency contributions to theta_it",
       subtitle = paste0("Faint lines: per-kid trajectories.  ",
                         "Steel ribbon: 80% population spread.  ",
                         "Gray band (total only): 10-90th percentile of item difficulty.")) +
  theme_minimal(base_size = 11) +
  theme(strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(size = 9))

ggsave(file.path(OUT_FIGS, "io_ability_decomposition.png"),
       p, width = 11, height = 9, dpi = 200)
cat("Wrote io_ability_decomposition.png\n")
