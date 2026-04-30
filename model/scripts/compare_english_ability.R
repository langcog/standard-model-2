## Compare per-child *ability trajectories* across English ablations.
##
## Re-factor the 2PL linear predictor as
##   eta_ijt = lambda_j * (theta_it - beta_j)
## with
##   theta_it = xi_i + log H + (1 + delta + zeta_i) * log((t - s) / a_0)
##   beta_j   = psi_j - log p_j
## So theta_it is the per-(child, age) ability on the logit scale, and
## beta_j is the per-item frequency-adjusted difficulty. This plot
## shows theta_it for each variant + the item-difficulty distribution
## as a reference; differences between variants reveal which model
## assumptions matter for "ability".
##
## Output: model/figs/longitudinal/english_ablations_ability.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(posterior); library(rstan)
})

VARIANTS <- c(
  "long_slopes"            = "Lean reference\n(slopes only)",
  "long_baseline"          = "Drop slopes\n(no zeta)",
  "long_fix_delta_slopes"  = "Pin delta = 0\n(no acceleration)",
  "long_free_s_slopes"     = "Free s\n(start time)",
  "long_2pl_slopes"        = "Add 2PL\n(item lambda)"
)

bundle <- load_dataset_bundle("english")
sd_    <- bundle$stan_data
log_p  <- log(bundle$word_info$prob)
log_H  <- sd_$log_H
a0     <- sd_$a0

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

AGE_GRID <- seq(11, 32, by = 0.25)

extract_theta_curves <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  draws <- as_draws_df(fit)

  med <- function(p) median(draws[[p]])
  s     <- med("s")
  delta <- med("delta")
  xi    <- sapply(seq_len(sd_$I),
                  function(i) median(draws[[sprintf("xi[%d]", i)]]))
  zeta  <- if ("zeta[1]" %in% names(draws))
            sapply(seq_len(sd_$I),
                   function(i) median(draws[[sprintf("zeta[%d]", i)]]))
          else rep(0, sd_$I)
  psi   <- sapply(seq_len(sd_$J),
                  function(j) median(draws[[sprintf("psi[%d]", j)]]))
  beta_j <- psi - log_p

  # theta_it on age grid for each kid
  theta_grid <- expand.grid(ii = seq_len(sd_$I),
                            age = AGE_GRID,
                            KEEP.OUT.ATTRS = FALSE)
  theta_grid$theta <- with(theta_grid, {
    ae <- pmax(age - s, 0.01)
    log_age <- log(ae / a0)
    xi[ii] + log_H + (1 + delta + zeta[ii]) * log_age
  })

  list(curves = theta_grid %>% mutate(variant = variant),
       beta = tibble(beta = beta_j, variant = variant))
}

cat("Extracting theta_it across variants...\n")
all_curves <- vector("list", length(VARIANTS))
all_beta   <- vector("list", length(VARIANTS))
for (k in seq_along(VARIANTS)) {
  v <- names(VARIANTS)[k]
  cat("  ", v, "\n")
  out <- extract_theta_curves(v)
  if (!is.null(out)) {
    all_curves[[k]] <- out$curves
    all_beta[[k]]   <- out$beta
  }
}
curves <- bind_rows(all_curves) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)))
beta_df <- bind_rows(all_beta) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)))

# Per-variant population mean theta
pop_mean <- curves %>%
  group_by(variant_label, age) %>%
  summarise(theta_mean = mean(theta),
            theta_lo10 = quantile(theta, 0.10),
            theta_hi90 = quantile(theta, 0.90),
            .groups = "drop")

# Item difficulty band: 10/50/90th percentiles of beta_j (per variant
# since psi differs slightly across variants)
beta_summary <- beta_df %>%
  group_by(variant_label) %>%
  summarise(beta_p10 = quantile(beta, 0.10),
            beta_p50 = quantile(beta, 0.50),
            beta_p90 = quantile(beta, 0.90),
            .groups = "drop")

# Plot: per-kid faint lines + pop mean + 10/90 ribbon, with item
# difficulty bands as horizontal lines.
p <- ggplot() +
  # Item-difficulty reference: shaded band from 10th to 90th
  # percentile of beta_j, with median as a darker line
  geom_rect(data = beta_summary,
            aes(ymin = beta_p10, ymax = beta_p90,
                xmin = -Inf, xmax = Inf),
            fill = "gray85", alpha = 0.4) +
  geom_hline(data = beta_summary,
             aes(yintercept = beta_p50),
             color = "gray50", linewidth = 0.4, linetype = "dashed") +
  # Per-kid theta trajectories (faint)
  geom_line(data = curves,
            aes(x = age, y = theta, group = ii),
            color = "steelblue", alpha = 0.10, linewidth = 0.3) +
  # Population mean + 80% spread
  geom_ribbon(data = pop_mean,
              aes(x = age, ymin = theta_lo10, ymax = theta_hi90),
              fill = "steelblue", alpha = 0.20) +
  geom_line(data = pop_mean,
            aes(x = age, y = theta_mean),
            color = "steelblue", linewidth = 0.9) +
  facet_wrap(~variant_label, ncol = 3) +
  coord_cartesian(ylim = c(-10, 25)) +
  labs(x = "age (months)",
       y = expression(theta[it] ~ " (logit scale)"),
       title = expression("Per-child ability trajectories: " ~
                          theta[it] ~ "= " ~ xi[i] ~ "+ log H + (1+" ~ delta ~ "+" ~ zeta[i] ~ ") log((t-s)/" ~ a[0] ~ ")"),
       subtitle = paste0("Faint lines: per-kid theta. Bold steel: pop. mean. ",
                         "Gray band: 10-90th percentile of item ",
                         "difficulty (",
                         "beta_j = psi_j - log p_j",
                         "); dashed: median item.")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(OUT_FIGS, "english_ablations_ability.png"),
       p, width = 12, height = 7.5, dpi = 200)
cat("Wrote english_ablations_ability.png\n")
