## Cross-model analysis of the M0..M5 nested family on English.
##
## Reads small summary / draws / loo artifacts from
## fits/summaries/ (synced from Sherlock).  Produces:
##
##   outputs/figs/longitudinal/nested_family_scalars.png  forest of headline
##                                                       scalars across models
##   outputs/figs/longitudinal/nested_family_loo.png      ELPD differences
##                                                       (cumulative, M_k - M_0)
##   outputs/figs/longitudinal/nested_family_summary.csv  scalar summary
##                                                       table for paper
##
## Robust to missing files: any model without a {summary,draws,loo}.rds
## is silently skipped, so we can re-run as more extractions land.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(loo)
})

SUM_DIR  <- file.path(PATHS$fits_dir, "summaries")
OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

FAMILY <- list(
  list(label = "M0", tag = "long_m0",                desc = "no time, no freq"),
  list(label = "M1", tag = "long_m1_time_only",      desc = "+ time only (no freq)"),
  list(label = "M2", tag = "long_m1",                desc = "+ frequency"),
  list(label = "M3", tag = "long_baseline",          desc = "+ free delta"),
  list(label = "M4", tag = "long_slopes",            desc = "+ per-child slopes"),
  list(label = "M5", tag = "long_class_beta_slopes", desc = "+ class-specific beta_c"),
  list(label = "M6", tag = "long_m5",                desc = "+ 2PL (lambda_j)")
)
EXTRA <- list(  # off-spine variants for later context
  list(label = "no_freq_M4", tag = "long_no_freq_slopes",
       desc = "M4 - log p_j (RQ2 robustness with full structure on top)"),
  list(label = "lmm_NO",  tag = "long_lmm_slopes_norwegian",
       desc = "LMM (Norwegian; off-spine)")
)
ALL <- c(FAMILY, EXTRA)

available <- function(tag, kind) {
  file.exists(file.path(SUM_DIR, sprintf("%s.%s.rds", tag, kind)))
}

# ---- 1. Scalar summary table ---- #
HEADLINE_PARS <- c("sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
                   "pi_alpha", "delta", "s", "beta_age", "sigma_lambda",
                   sprintf("beta_c[%d]", 1:4))

# Loader robust to both formats:
#   - new (posterior::summarise_draws tibble): columns variable, mean,
#     median, q025, q975, ess_bulk, rhat
#   - old (rstan summary matrix): rownames are variables, columns
#     include mean, 50%, 2.5%, 97.5%, n_eff, Rhat
load_summary <- function(path) {
  s <- readRDS(path)
  if (is.matrix(s)) {
    tibble(variable = rownames(s),
           mean   = s[, "mean"],
           median = s[, "50%"],
           q025   = s[, "2.5%"],
           q975   = s[, "97.5%"],
           ess_bulk = s[, "n_eff"],
           rhat   = s[, "Rhat"])
  } else {
    as_tibble(s) %>%
      transmute(variable, mean, median, q025, q975, ess_bulk, rhat)
  }
}

scalars <- bind_rows(lapply(ALL, function(m) {
  if (!available(m$tag, "summary")) return(NULL)
  load_summary(file.path(SUM_DIR, sprintf("%s.summary.rds", m$tag))) %>%
    filter(variable %in% HEADLINE_PARS) %>%
    mutate(label = m$label, tag = m$tag, desc = m$desc) %>%
    select(label, tag, variable, mean, median,
           lo95 = q025, hi95 = q975, ess_bulk, rhat)
}))

if (nrow(scalars) == 0) stop("No summary files found in ", SUM_DIR)

cat("\nLoaded scalars from", length(unique(scalars$label)), "models.\n")
cat("Headline scalar summary:\n")
print(as.data.frame(scalars %>%
                    mutate(cell = sprintf("%.2f [%.2f, %.2f]",
                                           median, lo95, hi95)) %>%
                    select(label, variable, cell) %>%
                    pivot_wider(names_from = variable, values_from = cell)),
      row.names = FALSE, max = 200)

# Write CSV for paper
write.csv(scalars,
          file.path(OUT_FIGS, "nested_family_summary.csv"),
          row.names = FALSE)

# ---- 2. Forest plot ---- #
KEY_FOREST <- c("sigma_alpha", "pi_alpha", "delta", "sigma_zeta",
                "sigma_lambda")
sf <- scalars %>%
  filter(variable %in% KEY_FOREST,
         label %in% sapply(FAMILY, function(m) m$label)) %>%
  mutate(label    = factor(label, levels = sapply(FAMILY, function(m) m$label)),
         variable = factor(variable, levels = KEY_FOREST))

p_forest <- ggplot(sf, aes(x = median, y = label,
                           xmin = lo95, xmax = hi95)) +
  geom_pointrange(color = "steelblue", size = 0.4) +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  scale_y_discrete(limits = rev) +
  labs(x = "posterior median, 95% CrI", y = NULL,
       title = "Nested family M0..M5 — headline scalars (English)",
       subtitle = "Each row = one model. Empty cells: parameter not present in that model.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_FIGS, "nested_family_scalars.png"),
       p_forest, width = 13, height = 4.5, dpi = 200)
cat(sprintf("\nWrote %s\n",
            file.path(OUT_FIGS, "nested_family_scalars.png")))

# ---- 3. LOO comparison ---- #
loo_objs <- list()
for (m in FAMILY) {
  p <- file.path(SUM_DIR, sprintf("%s.loo.rds", m$tag))
  if (file.exists(p)) loo_objs[[m$label]] <- readRDS(p)
}
# Off-spine variants (no_freq, etc.) for targeted RQ2/RQ3 tests
extra_loo <- list()
for (m in EXTRA) {
  p <- file.path(SUM_DIR, sprintf("%s.loo.rds", m$tag))
  if (file.exists(p)) extra_loo[[m$label]] <- readRDS(p)
}

if (length(loo_objs) >= 2) {
  cat(sprintf("\nLOO available for %d models: %s\n",
              length(loo_objs),
              paste(names(loo_objs), collapse = ", ")))

  # Pairwise comparison (loo_compare picks best as ref, reports
  # ELPD differences relative to it, with SE).
  comp <- loo_compare(loo_objs)
  cat("\nloo_compare (rows ordered best -> worst):\n")
  print(comp, simplify = FALSE)

  # Build a "M_k vs M_(k-1)" table for the spine
  spine_labels <- sapply(FAMILY, function(m) m$label)
  spine_present <- spine_labels[spine_labels %in% names(loo_objs)]
  if (length(spine_present) >= 2) {
    pair_diffs <- list()
    for (i in 2:length(spine_present)) {
      a <- spine_present[i - 1]; b <- spine_present[i]
      d <- loo_compare(list(a = loo_objs[[a]], b = loo_objs[[b]]))
      # Row 'b' minus row 'a' if a is ref (worse); we want b - a.
      ref <- rownames(d)[1]; other <- rownames(d)[2]
      diff <- if (ref == "a") d[other, "elpd_diff"] else -d[other, "elpd_diff"]
      se   <- d[other, "se_diff"]
      pair_diffs[[length(pair_diffs) + 1]] <- tibble(
        step = sprintf("%s vs %s", b, a),
        elpd_diff = diff,
        se_diff = se,
        z = diff / se
      )
    }
    pair <- bind_rows(pair_diffs)
    cat("\nStep-wise ELPD differences along M0 -> M5 spine (b minus a):\n")
    print(as.data.frame(pair), row.names = FALSE, digits = 3)

    # Persist both the full loo_compare ranking and the step-wise table
    # as CSVs alongside the figure (the header comment promises these).
    rank_df <- as.data.frame(comp)
    rank_df$model <- rownames(rank_df)
    rank_df <- rank_df[, c("model", setdiff(colnames(rank_df), "model"))]
    write.csv(rank_df,
              file.path(OUT_FIGS, "nested_family_loo_ranking.csv"),
              row.names = FALSE)
    write.csv(pair,
              file.path(OUT_FIGS, "nested_family_loo_steps.csv"),
              row.names = FALSE)

    # Off-spine pairwise tests. The spine itself now has a clean RQ2
    # test at M2 vs M1 (+ frequency on top of time-only); we keep the
    # M4-level no_freq comparison as a robustness check that the
    # frequency-redundancy result is stable under more model structure.
    if ("no_freq_M4" %in% names(extra_loo) && "M4" %in% names(loo_objs)) {
      d <- loo_compare(list(M4 = loo_objs[["M4"]],
                            no_freq = extra_loo[["no_freq_M4"]]))
      ref <- rownames(d)[1]; other <- rownames(d)[2]
      diff <- if (ref == "M4") -d[other, "elpd_diff"] else d[other, "elpd_diff"]
      se <- d[other, "se_diff"]
      cat(sprintf(
        "\nRQ2 robustness (no_freq_M4 vs M4): ELPD_diff = %.1f +/- %.1f (z=%.1f)\n",
        diff, se, diff / se))
      cat("   (Drops log p_j on top of the slopes-level model.\n")
      cat("    Companion to the on-spine M2 vs M1 test.)\n")
    }

    # Plot
    pair_p <- pair %>%
      mutate(step = factor(step, levels = step))
    p_loo <- ggplot(pair_p, aes(x = elpd_diff, y = step,
                                xmin = elpd_diff - 1.96 * se_diff,
                                xmax = elpd_diff + 1.96 * se_diff)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      geom_pointrange(color = "firebrick", size = 0.4) +
      geom_text(aes(label = sprintf("%.0f +/- %.0f  (z=%.1f)",
                                     elpd_diff, se_diff, z)),
                hjust = -0.1, vjust = -0.4, size = 3, color = "gray30") +
      scale_y_discrete(limits = rev) +
      labs(x = "ELPD_diff (later model minus earlier; >0 = later wins)",
           y = NULL,
           title = "Step-wise LOO ELPD differences along M0 -> M5",
           subtitle = "Each step adds one component. >0 means the new component improves out-of-sample fit.") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(OUT_FIGS, "nested_family_loo.png"),
           p_loo, width = 10, height = 4, dpi = 200)
    cat(sprintf("\nWrote %s\n",
                file.path(OUT_FIGS, "nested_family_loo.png")))
  }
} else {
  cat("\nLOO objects available for fewer than 2 models; skipping ELPD plot.\n")
  cat("(Loaded:", paste(names(loo_objs), collapse = ", "), ")\n")
}
