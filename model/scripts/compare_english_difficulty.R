## Diagnostic counterpart to compare_english_ability.R: how does the
## difficulty term beta_j = psi_j - log p_j vary across the ablation
## set? Prediction: nearly identical across all five variants, because
## *no current ablation directly touches beta_j*. This figure documents
## that gap and motivates planned beta-side ablations (free_beta,
## class_beta, no_class).
##
## Two panels in one figure:
##   (1) Density of beta_j across all items, one curve per variant.
##   (2) Per-item beta_j scatter vs. lean reference (each variant
##       plotted against the lean-ref baseline).
##
## Output: model/figs/longitudinal/english_ablations_difficulty.png

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

VARIANT_COLORS <- c(
  "lean ref"    = "#1f77b4",
  "drop slopes" = "#aec7e8",
  "add 2PL"     = "#7f7f7f",
  "pin delta=0" = "#d62728",
  "free s"      = "#ff7f0e"
)

bundle <- load_dataset_bundle("english")
sd_    <- bundle$stan_data
log_p  <- log(bundle$word_info$prob)
class_  <- bundle$word_info$cc

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

extract_beta <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  draws <- as_draws_df(fit)
  psi <- sapply(seq_len(sd_$J),
                function(j) median(draws[[sprintf("psi[%d]", j)]]))
  tibble(variant = variant, jj = seq_len(sd_$J),
         class = class_, psi = psi, log_p = log_p, beta = psi - log_p)
}

cat("Extracting beta_j across variants...\n")
betas <- bind_rows(lapply(names(VARIANTS), function(v) {
  cat("  ", v, "\n"); extract_beta(v)
})) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)))

# ---- Panel 1: density of beta_j per variant ---- #
p_dens <- ggplot(betas, aes(x = beta, color = variant_label,
                            group = variant_label)) +
  geom_density(linewidth = 0.9, alpha = 0.6) +
  scale_color_manual(values = VARIANT_COLORS, name = "variant") +
  labs(x = expression(beta[j] ~ "= " ~ psi[j] ~ "- log " ~ p[j]),
       y = "density (across J=200 items)",
       title = expression("(1) Distribution of item difficulty " ~ beta[j] ~ " by variant"),
       subtitle = "Densities should overlap nearly perfectly: no current ablation touches the difficulty term.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# ---- Panel 2: beta_j scatter against lean-ref baseline ---- #
ref <- betas %>% filter(variant == "long_slopes") %>%
  select(jj, beta_ref = beta, class)
non_ref <- betas %>%
  filter(variant != "long_slopes") %>%
  left_join(ref %>% select(jj, beta_ref), by = "jj")

cors <- non_ref %>% group_by(variant_label) %>%
  summarise(r = cor(beta, beta_ref),
            mae = mean(abs(beta - beta_ref)),
            .groups = "drop") %>%
  mutate(label = sprintf("r=%.3f\nMAE=%.2f", r, mae))

p_scatter <- ggplot(non_ref,
                    aes(x = beta_ref, y = beta, color = variant_label)) +
  geom_abline(slope = 1, intercept = 0, color = "gray50",
              linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.45, size = 0.9) +
  geom_text(data = cors, aes(label = label),
           x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
           inherit.aes = FALSE, size = 3, color = "gray20") +
  facet_wrap(~variant_label, nrow = 1) +
  scale_color_manual(values = VARIANT_COLORS, guide = "none") +
  labs(x = expression(beta[j] ~ "  (lean reference)"),
       y = expression(beta[j] ~ "  (variant)"),
       title = expression("(2) Per-item " ~ beta[j] ~ " against lean-reference baseline"),
       subtitle = "Off-diagonal scatter would mean the variant is shifting item difficulties.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

p_combined <- p_dens / p_scatter + plot_layout(heights = c(1, 1.1))

ggsave(file.path(OUT_FIGS, "english_ablations_difficulty.png"),
       p_combined, width = 11, height = 8, dpi = 200)
cat("Wrote english_ablations_difficulty.png\n")

# Print correlations for quick verification
cat("\nCorrelations of beta_j vs lean reference:\n")
print(as.data.frame(cors %>% select(variant_label, r, mae)),
      row.names = FALSE)
