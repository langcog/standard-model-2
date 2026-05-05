## Analytical σ_r sensitivity for M_best fits.
##
## In the model, σ_ξ² = σ_α² + σ_r². The data identifies σ_ξ²
## (per the §4b finding); σ_r is pinned externally. So
##   π_α(σ_r) = 1 − σ_r² / σ_ξ²
## given the σ_ξ posterior. This is "free" in the sense that we don't
## have to refit at different σ_r priors — we can compute the implied
## π_α from any single fit's σ_ξ posterior, evaluated at arbitrary σ_r.
##
## Validity check: §4b ran 4 actual refits at σ_r ∈ {0.3, 0.5, 0.8, 1.2}
## on the cross-sectional 2PL English fit. Their σ_α posteriors give
## σ_ξ² ≈ 4.8 essentially constant; π_α follows 1 − σ_r²/4.8 exactly.
## So the analytical extrapolation is well-supported there. This script
## extends it to M_best (longitudinal, slopes, no_freq) on English and
## Norwegian.
##
## Run:  Rscript model/scripts/sigma_r_analytical_sensitivity.R
## Output:
##   outputs/figs/longitudinal/sigma_r_analytical_sensitivity.png
##   outputs/figs/longitudinal/sigma_r_analytical_sensitivity.csv

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Fits to evaluate (M_best variants for the two production samples)
FITS <- list(
  english_no_freq = list(
    tag    = "long_no_freq_slopes",
    label  = "English, M_best (no_freq slopes)",
    n_kids = 200
  ),
  english_slopes = list(
    tag    = "long_slopes",
    label  = "English (slopes, +freq)",
    n_kids = 200
  ),
  norwegian_no_freq = list(
    tag    = "long_no_freq_slopes_norwegian",
    label  = "Norwegian, M_best (no_freq slopes)",
    n_kids = 200
  ),
  norwegian_slopes = list(
    tag    = "long_slopes_norwegian",
    label  = "Norwegian (slopes, +freq)",
    n_kids = 200
  ),
  # Cross-sectional 2PL §4b reference (validates that the analytical
  # extrapolation matches actual refits — σ_xi² ≈ 4.8, π_α curve hits
  # observed 0.98 / 0.94 / 0.87 / 0.70 at σ_r = 0.30 / 0.53 / 0.80 / 1.20).
  csec_2pl_ref = list(
    tag    = "wordbank_2pl",
    label  = "English §4b 2PL cross-sectional (validation)",
    n_kids = NA
  )
)

# σ_r grid: span the same range as §4b plus a couple of intermediate
SIGMA_R_GRID <- c(0.30, 0.40, 0.50, 0.534, 0.65, 0.80, 1.00, 1.20)

# §4b empirical anchors (cross-sectional 2PL, full data) for the validation
SECTION_4B <- tibble::tibble(
  fit_name = "csec_2pl_ref",
  sigma_r  = c(0.30, 0.534, 0.80, 1.20),
  sigma_alpha = c(2.18, 2.13, 2.05, 1.83),
  pi_alpha    = c(0.98, 0.94, 0.87, 0.70),
  source = "§4b refit"
)

# ---- Compute analytical π_α(σ_r) per fit, propagating σ_ξ posterior --
get_sigma_xi_draws <- function(tag) {
  pth <- file.path(PATHS$fits_dir, "summaries", paste0(tag, ".draws.rds"))
  if (!file.exists(pth)) {
    warning("Missing draws: ", pth); return(NULL)
  }
  d <- readRDS(pth)
  # σ_ξ may be present directly, or computable from σ_α and σ_r prior
  if ("sigma_xi" %in% names(d)) {
    return(as.numeric(unlist(d$sigma_xi)))
  }
  if ("sigma_alpha" %in% names(d)) {
    # If sigma_xi not in scalar pars, reconstruct from sigma_alpha + the
    # fit's σ_r prior. Need to load stan_data for that fit.
    # For simplicity require sigma_xi in the draws.
    warning("No sigma_xi in ", tag, "; using sigma_alpha alone (may misrepresent if σ_r differed).")
    return(NULL)
  }
  NULL
}

results <- list()
for (nm in names(FITS)) {
  cfg <- FITS[[nm]]
  sxi <- get_sigma_xi_draws(cfg$tag)
  if (is.null(sxi)) next
  rows <- lapply(SIGMA_R_GRID, function(sr) {
    pi_a <- 1 - (sr^2) / (sxi^2)
    pi_a <- pmax(pi_a, 0)  # negative just means σ_r prior is incompatible with data
    sigma_alpha_implied <- sqrt(pmax(sxi^2 - sr^2, 0))
    tibble::tibble(
      fit_name = nm, label = cfg$label,
      sigma_r = sr,
      sigma_xi_med = median(sxi),
      sigma_xi_q025 = quantile(sxi, 0.025),
      sigma_xi_q975 = quantile(sxi, 0.975),
      sigma_alpha_implied_med = median(sigma_alpha_implied),
      sigma_alpha_implied_q025 = quantile(sigma_alpha_implied, 0.025),
      sigma_alpha_implied_q975 = quantile(sigma_alpha_implied, 0.975),
      pi_alpha_med = median(pi_a),
      pi_alpha_q025 = quantile(pi_a, 0.025),
      pi_alpha_q975 = quantile(pi_a, 0.975)
    )
  })
  results[[nm]] <- bind_rows(rows)
}

df <- bind_rows(results)

# ---- Save table ------------------------------------------------------
out_csv <- file.path(OUT_DIR, "sigma_r_analytical_sensitivity.csv")
write.csv(df |> select(-label), out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")

# Compact summary table for the writeup
compact <- df |>
  mutate(pi_alpha = sprintf("%.2f [%.2f, %.2f]",
                            pi_alpha_med, pi_alpha_q025, pi_alpha_q975)) |>
  select(fit_name, sigma_r, pi_alpha) |>
  pivot_wider(names_from = sigma_r, values_from = pi_alpha,
              names_prefix = "σ_r=")
cat("\n=== π_α as a function of σ_r (analytical, M_best variants) ===\n")
print(as.data.frame(compact))

# Validation against §4b: for the csec_2pl_ref fit, compare analytical
# vs §4b empirical
val <- df |>
  filter(fit_name == "csec_2pl_ref") |>
  inner_join(SECTION_4B |> rename(pi_alpha_empirical = pi_alpha,
                                  sigma_alpha_empirical = sigma_alpha),
             by = c("fit_name", "sigma_r")) |>
  select(sigma_r, sigma_xi_med,
         pi_alpha_med, pi_alpha_empirical,
         sigma_alpha_implied_med, sigma_alpha_empirical) |>
  mutate(diff_pi_alpha = pi_alpha_med - pi_alpha_empirical,
         diff_sigma_alpha = sigma_alpha_implied_med - sigma_alpha_empirical)
cat("\n=== Validation: analytical (csec_2pl_ref) vs §4b refits ===\n")
print(as.data.frame(val))

# ---- Plot ------------------------------------------------------------
plot_df <- df |>
  filter(fit_name != "csec_2pl_ref") |>
  mutate(label = factor(label, levels = unique(df$label)))

ref_pts <- SECTION_4B |>
  rename(pi_alpha_med = pi_alpha) |>
  mutate(label = "English §4b 2PL cross-sectional (refit)",
         pi_alpha_q025 = pi_alpha_med, pi_alpha_q975 = pi_alpha_med)

p <- ggplot(plot_df, aes(sigma_r, pi_alpha_med, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = pi_alpha_q025, ymax = pi_alpha_q975),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.85) +
  geom_point(data = ref_pts, aes(sigma_r, pi_alpha_med),
             colour = "grey15", fill = "grey15", size = 2.5,
             shape = 21, inherit.aes = FALSE) +
  geom_vline(xintercept = 0.534, colour = "grey55",
             linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = 0.534, y = 0.62, hjust = -0.05, size = 3.0,
           colour = "grey40",
           label = "external σ_r = 0.534\n(Sperry / HR / WF)") +
  scale_x_continuous(breaks = c(0.3, 0.5, 0.534, 0.8, 1.0, 1.2),
                     labels = c("0.30", "0.50", "0.534", "0.80", "1.00", "1.20")) +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  labs(
    x = expression("Externally pinned " * sigma[r]),
    y = expression(pi[alpha]),
    colour = NULL, fill = NULL,
    title = expression("Analytical " * sigma[r] * " sensitivity for M_best"),
    subtitle = "Posterior median (line) and 95% CrI (band) of π_α = 1 − σ_r²/σ_ξ², propagating σ_ξ posterior. Black points: §4b refits at matched σ_r values (validation that analytical and refit match)."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, colour = "grey25"))

out_png <- file.path(OUT_DIR, "sigma_r_analytical_sensitivity.png")
ggsave(out_png, p, width = 9.5, height = 5.5, dpi = 150)
cat("\nWrote:", out_png, "\n")
