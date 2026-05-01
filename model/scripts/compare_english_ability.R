## Diagnostic: which ablations actually affect the *ability* term?
##
## The 2PL form factors as
##   eta_ijt = lambda_j * (theta_it - beta_j)
## where
##   theta_it = xi_i + log H + (1 + delta + zeta_i) * log((t - s) / a_0)
##   beta_j   = psi_j - log p_j
##
## Population-mean ability (sum-to-zero centering means mean(zeta_i)=0):
##   theta_pop(t) = mean(xi_i) + log H + (1 + delta) * log((t - s) / a_0)
## Only delta and s appear, so ablations that change neither -- drop
## slopes, add 2PL -- should give *identical* population-mean curves to
## the lean reference. pin delta=0 should flatten; free s should shift
## the kink.
##
## Two panels in one figure:
##   (1) Population-mean theta_pop(t) overlaid across 5 variants.
##   (2) SD of theta_it across kids vs age (per-child spread).
##       Reveals where slopes / pin-delta really differ.
##
## Output: model/figs/longitudinal/english_ablations_ability.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

VARIANTS <- c(
  "long_slopes"            = "lean ref",
  "long_baseline"          = "drop slopes",
  "long_fix_delta_slopes"  = "pin delta=0",
  "long_free_s_slopes"     = "free s",
  "long_2pl_slopes"        = "add 2PL"
)

# Color: emphasize that drop slopes and add 2PL should *not* differ
# from lean ref on the population-mean curve. Group them visually.
VARIANT_COLORS <- c(
  "lean ref"    = "#1f77b4",   # blue: ability-side reference
  "drop slopes" = "#aec7e8",   # light blue: same pop mean
  "add 2PL"     = "#7f7f7f",   # gray: discrimination, neither side
  "pin delta=0" = "#d62728",   # red: ability slope ablation
  "free s"      = "#ff7f0e"    # orange: ability start ablation
)

bundle <- load_dataset_bundle("english")
sd_    <- bundle$stan_data
log_p  <- log(bundle$word_info$prob)
log_H  <- sd_$log_H
a0     <- sd_$a0

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

# Plot from age 0 to 32 -- but each variant's curve only starts at its
# own s (the model's accumulation onset). Makes free-s shift visible.
AGE_GRID <- seq(0.0, 32, by = 0.25)

extract_theta <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  draws <- as_draws_df(fit)

  med <- function(p) median(draws[[p]])
  s     <- med("s")
  delta <- med("delta")
  xi    <- sapply(seq_len(sd_$I),
                  function(i) median(draws[[sprintf("xi[%d]", i)]]))
  mu_xi <- mean(xi)
  zeta  <- if ("zeta[1]" %in% names(draws))
            sapply(seq_len(sd_$I),
                   function(i) median(draws[[sprintf("zeta[%d]", i)]]))
          else rep(0, sd_$I)

  # Only define theta where t > s (the model's accumulation onset).
  # Below s, eta is set to a large negative value in code; conceptually
  # there is no ability trajectory pre-s.
  in_range <- AGE_GRID > s + 0.05
  age_v    <- AGE_GRID[in_range]
  log_age  <- log((age_v - s) / a0)

  pop_mean <- mu_xi + log_H + (1 + delta) * log_age

  # Per-kid spread: SD of theta_it across kids at each age.
  theta_kids <- outer(xi, rep(1, length(age_v))) +
                log_H +
                outer(1 + delta + zeta, log_age)
  spread <- apply(theta_kids, 2, sd)

  tibble(variant = variant,
         age = age_v,
         pop_mean = pop_mean,
         spread = spread,
         s_med = s,
         delta_med = delta)
}

cat("Extracting ability summaries...\n")
abilities <- bind_rows(lapply(names(VARIANTS), function(v) {
  cat("  ", v, "\n"); extract_theta(v)
})) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)))

# Item-difficulty band (from lean ref) for visual reference
ref_fit <- readRDS(file.path(PATHS$fits_dir, "long_slopes.rds"))
ref_draws <- as_draws_df(ref_fit)
psi_ref <- sapply(seq_len(sd_$J),
                  function(j) median(ref_draws[[sprintf("psi[%d]", j)]]))
beta_ref <- psi_ref - log_p
beta_q <- quantile(beta_ref, c(0.10, 0.50, 0.90))

# Annotate s for free_s variant (need to mark its kink)
s_free <- abilities %>% filter(variant == "long_free_s_slopes") %>%
  slice(1) %>% pull(s_med)
s_lean <- abilities %>% filter(variant == "long_slopes") %>%
  slice(1) %>% pull(s_med)

# ---- Panel 1: Population-mean theta_pop(t) overlay ---- #
p_mean <- ggplot(abilities,
                 aes(x = age, y = pop_mean, color = variant_label)) +
  # Item difficulty band as background reference
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = beta_q[1], ymax = beta_q[3],
           fill = "gray85", alpha = 0.5) +
  annotate("segment", x = -Inf, xend = Inf,
           y = beta_q[2], yend = beta_q[2],
           color = "gray55", linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = VARIANT_COLORS, name = "variant") +
  # Mark s for free-s and lean reference
  geom_vline(xintercept = s_free, linetype = "dotted",
             color = "#ff7f0e", alpha = 0.8) +
  geom_vline(xintercept = s_lean, linetype = "dotted",
             color = "#1f77b4", alpha = 0.6) +
  annotate("text", x = s_free, y = -10, label = sprintf("s=%.1f", s_free),
           color = "#ff7f0e", hjust = -0.1, vjust = 0, size = 3.2) +
  annotate("text", x = s_lean, y = -10, label = sprintf("s=%.1f", s_lean),
           color = "#1f77b4", hjust = -0.1, vjust = 0, size = 3.2) +
  coord_cartesian(xlim = c(0, 32), ylim = c(-12, 22)) +
  labs(x = "age (months)",
       y = expression("population-mean " ~ theta[pop]~"(t)  (logit scale)"),
       title = "(1) Population-mean ability trajectory",
       subtitle = expression("Only " ~ delta ~ " and s appear in " ~ theta[pop]*";  lean / drop-slopes / 2PL should overlay.  Gray band = 10-90% of item " ~ beta[j])) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# ---- Panel 2: spread of theta_it across kids vs age ---- #
SPREAD_CAP <- 30   # clip; pin-delta=0 goes off-scale
abilities_cap <- abilities %>%
  mutate(spread_capped = pmin(spread, SPREAD_CAP),
         off_scale     = spread > SPREAD_CAP)

# Off-scale annotation for pin delta=0
off_scale_var <- abilities_cap %>% filter(off_scale) %>%
  count(variant_label) %>% pull(variant_label)

p_spread <- ggplot(abilities_cap,
                   aes(x = age, y = spread_capped, color = variant_label)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = VARIANT_COLORS, name = "variant") +
  coord_cartesian(xlim = c(0, 32), ylim = c(0, SPREAD_CAP)) +
  annotate("text", x = 16, y = SPREAD_CAP * 0.95,
           label = "pin delta=0 saturates at top of scale\n(SD reaches ~1000 at young ages)",
           hjust = 0.5, vjust = 1, size = 3, color = "#d62728") +
  labs(x = "age (months)",
       y = expression("SD of " ~ theta[it] ~ " across kids  (capped at 30)"),
       title = "(2) Per-child ability spread vs. age",
       subtitle = expression("How wide is the population at age t?  drop-slopes => flat (no " ~ zeta[i] ~ ").  pin " ~ delta ~ "=0 => " ~ sigma[zeta] ~ " inflates wildly.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

p_combined <- p_mean / p_spread + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUT_FIGS, "english_ablations_ability.png"),
       p_combined, width = 10, height = 8, dpi = 200)
cat("Wrote english_ablations_ability.png\n")
