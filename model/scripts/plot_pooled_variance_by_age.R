## Between-child variance composition by age for the pooled IO ladder.
##
## For each of baseline / gamma-additive / gamma-multiplicative,
## decompose Var_i(theta_i(t)) at log_age L = log(t / a_0) into:
##   * Input intercept            sigma_r^2
##   * Input slope (gamma)        L^2 * (gamma * sigma_r)^2  (+interaction in mult)
##   * Efficiency intercept       sigma_alpha^2
##   * Efficiency slope           L^2 * sigma_zeta^2
##   * Cross (input <-> slope)    2 L * Cov(xi, kappa)
## and plot 5-component stacked bars at AGES = {12, 16, 19, 24, 30}.
## Combined "Input" share = sum of input-int + input-slope + cross.
##
## NOTE: uses posterior medians of sigma_r, sigma_alpha, sigma_zeta,
## gamma, and A=mean(kappa_study). The MULTIPLICATIVE partition scales
## with A^2, so will rebase when the wide-delta refit lands.
##
## Output: figs/io/pooled_variance_by_age.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

fits <- list(
  baseline = readRDS(file.path(PATHS$fits_dir, "io_pooled_widedelta.rds")),
  add      = readRDS(file.path(PATHS$fits_dir, "io_pooled_gamma_widedelta_add.rds"))
  # mult excluded: wide-delta mult fit failed to mix (see experiments.md §30).
)

A0   <- 15                       # pivot age used in the IO bundle
AGES <- c(12, 16, 19, 24, 30)
L    <- log(AGES / A0)

get_p <- function(fit, model) {
  d <- fit$draws(format = "df")
  ks <- as.matrix(d[, grep("^kappa_study\\[", colnames(d), value = TRUE)])
  list(
    sigma_r     = median(d$sigma_r),
    sigma_alpha = median(d$sigma_alpha),
    sigma_zeta  = median(d$sigma_zeta),
    gamma       = if (model == "baseline") 0 else median(d$gamma),
    A           = median(rowMeans(ks))
  )
}

decompose <- function(p, model, L) {
  sr2 <- p$sigma_r^2; sa2 <- p$sigma_alpha^2; sz2 <- p$sigma_zeta^2
  g <- p$gamma; A <- p$A
  input_int <- rep(sr2, length(L))
  eff_int   <- rep(sa2, length(L))
  eff_slope <- L^2 * sz2
  if (model == "baseline") {
    input_slope <- rep(0, length(L))
    cross       <- rep(0, length(L))
  } else if (model == "add") {
    input_slope <- L^2 * g^2 * sr2
    cross       <- 2 * L * g * sr2
  } else { # mult
    # mult decomposition derived analytically (see header comment):
    input_slope <- g^2 * L^2 * sr2 * (A^2 + sz2)
    cross       <- 2 * L * g * A * sr2
  }
  total <- input_int + input_slope + eff_int + eff_slope + cross
  data.frame(model = model, age = AGES,
             `Input intercept`    = input_int,
             `Input slope`        = input_slope,
             `Efficiency intercept` = eff_int,
             `Efficiency slope`   = eff_slope,
             `Cross (input link)` = cross,
             Total = total,
             check.names = FALSE)
}

# Build
res <- list()
for (m in names(fits)) {
  p <- get_p(fits[[m]], m)
  res[[m]] <- decompose(p, m, L) |> mutate(p_input_share =
    (`Input intercept` + `Input slope` + `Cross (input link)`) / Total)
}
df <- bind_rows(res)
cat("=== variance partition by age and model ===\n")
print(df, digits = 3)

# Long form for stacking
comp_cols <- c("Input intercept","Input slope","Efficiency intercept",
               "Efficiency slope","Cross (input link)")
df_long <- df |>
  pivot_longer(all_of(comp_cols), names_to = "component", values_to = "var") |>
  mutate(component = factor(component, levels = comp_cols),
         model     = factor(model, levels = c("baseline","add"),
                            labels = c("baseline (no γ)", "γ additive")))
df$model <- factor(df$model, levels = c("baseline","add"),
                   labels = c("baseline (no γ)", "γ additive"))

# Palette: 2 input shades + 2 efficiency shades + cross
fill_pal <- c(
  "Input intercept"      = "#B4D7E8",   # light blue
  "Input slope"          = "#1F77B4",   # darker blue
  "Efficiency intercept" = "#003D7A",   # navy
  "Efficiency slope"     = "#E8A0A8",   # pink
  "Cross (input link)"   = "#F0AB58"    # orange
)

# Compute label positions for input-share at the top of bars
labels_df <- df |>
  mutate(label = sprintf("Var=%.1f\ninput %.0f%%",
                         Total, 100 * p_input_share))

p <- ggplot(df_long, aes(factor(age), var, fill = component)) +
  geom_col(width = 0.7, color = "white") +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
  geom_text(data = labels_df,
            aes(factor(age), Total + 0.4, label = label),
            inherit.aes = FALSE, size = 2.9, vjust = 0, lineheight = 0.95) +
  facet_grid(. ~ model) +
  scale_fill_manual(values = fill_pal, name = NULL) +
  labs(x = "Age (months)",
       y = expression(paste("Variance contribution to ", Var[i](theta[it]))),
       title = "Between-child variance decomposition by age — pooled IO ladder",
       subtitle = sprintf(
         "a_0 = %d mo (log_age = 0 at age %d). Two input shades reveal where σ_r enters each model.",
         A0, A0)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

out <- file.path("figs", "io", "pooled_variance_by_age.png")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 14, height = 6.2, dpi = 150)
cat(sprintf("\nWrote %s\n", out))
