## Produce RQ-aligned plots and posterior summaries from the fit variants.
##
## Usage:   Rscript model/scripts/04_analyze.R [variant]
##          default variant = "baseline". Pass "all" to loop over all.
## Reads:   fits/wordbank_{variant}.rds, fits/subset_data.rds
## Outputs: model/figs/wordbank_{variant}_*.png
##          Prints scalar posteriors, class means, RQ1 correlations.

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/ppc.R")

args <- commandArgs(trailingOnly = TRUE)
variants_arg <- if (length(args) >= 1) args[1] else "baseline"
variants <- if (variants_arg == "all") {
  # The full set of registered variants in variant_hyperpriors().
  # Lean default first; legacy/full variants kept for re-loading old fits.
  c("baseline", "slopes", "2pl", "2pl_slopes",
    "free_s", "free_s_slopes",
    "fix_delta", "fix_s", "both_fixed", "2pl_fix_delta")
} else if (variants_arg == "ras-2pl") {
  c("baseline", "2pl")
} else if (variants_arg == "slopes-sweep") {
  c("baseline", "2pl", "slopes", "2pl_slopes")
} else {
  variants_arg
}

bundle <- readRDS(file.path(PATHS$fits_dir, "subset_data.rds"))

for (nm in variants) {
  fit_path <- file.path(PATHS$fits_dir, sprintf("wordbank_%s.rds", nm))
  if (!file.exists(fit_path)) {
    message(sprintf("[%s] fit not found at %s, skipping.", nm, fit_path))
    next
  }
  fit <- readRDS(fit_path)
  cat(sprintf("\n===== %s =====\n", nm))

  cat("\n-- Scalar posteriors --\n")
  print(summarize_fit(fit), digits = 3)

  cat("\n-- Class-level thresholds --\n")
  ctbl <- class_threshold_table(fit, bundle$class_levels)
  print(ctbl, n = Inf)

  psi_df <- extract_psi_df(fit, bundle$word_info, bundle$class_levels)
  r_glb <- cor(psi_df$psi_median, psi_df$log_p)
  cat(sprintf("\n-- RQ1 --\n  r(psi, log_p) = %.3f   R^2 = %.3f\n",
              r_glb, r_glb^2))
  print(psi_df %>% group_by(class) %>%
        summarise(n = n(), r = cor(psi_median, log_p), R2 = r^2), n = Inf)

  # plots
  plot_psi_vs_logp(psi_df,
    save_path = file.path(PATHS$figs_dir,
                         sprintf("wordbank_%s_psi_vs_logp.png", nm)),
    tag = nm)
  plot_class_means(ctbl,
    save_path = file.path(PATHS$figs_dir,
                         sprintf("wordbank_%s_class_means.png", nm)),
    tag = nm)
  plot_posterior_density(fit, "pi_alpha",
    save_path = file.path(PATHS$figs_dir,
                         sprintf("wordbank_%s_pi_alpha.png", nm)),
    tag = nm, vline = median(as_draws_df(fit)$pi_alpha), xlim = c(0, 1))
  plot_posterior_density(fit, "s",
    save_path = file.path(PATHS$figs_dir,
                         sprintf("wordbank_%s_s.png", nm)),
    tag = nm, xlab = "s (months)")
  plot_posterior_density(fit, "delta",
    save_path = file.path(PATHS$figs_dir,
                         sprintf("wordbank_%s_delta.png", nm)),
    tag = nm, vline = 0)

  # PPC suite: 4 panels saved individually + a combined overview.
  ppc_suite(fit, bundle,
    save_dir = PATHS$figs_dir,
    variant_tag = nm)
}

cat("\nDone. Figures in model/figs/\n")
