## Cross-fit results summary across the full ablation set, for resuming
## work after the prior session crashed mid-summary.
##
## Pulls scalar posteriors from every fit currently on disk:
##   - English longitudinal ablations (5)
##   - Norwegian longitudinal ablations (5)
##   - IO datasets (BabyView, SEEDLingS)
##   - Peekbank-linked proc fit
##
## Prints wide tables to stdout, saves a long-format CSV to
##   outputs/figs/results_summary.csv

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(posterior); library(rstan)
})

FITS <- list(
  # English ablations
  list(group = "English",   variant = "long_slopes",                label = "lean ref"),
  list(group = "English",   variant = "long_baseline",              label = "drop slopes"),
  list(group = "English",   variant = "long_fix_delta_slopes",      label = "pin delta=0"),
  list(group = "English",   variant = "long_free_s_slopes",         label = "free s"),
  list(group = "English",   variant = "long_2pl_slopes",            label = "add 2PL"),
  # Norwegian ablations
  list(group = "Norwegian", variant = "long_slopes_norwegian",            label = "lean ref"),
  list(group = "Norwegian", variant = "long_baseline_norwegian",          label = "drop slopes"),
  list(group = "Norwegian", variant = "long_fix_delta_slopes_norwegian",  label = "pin delta=0"),
  list(group = "Norwegian", variant = "long_free_s_slopes_norwegian",     label = "free s"),
  list(group = "Norwegian", variant = "long_2pl_slopes_norwegian",        label = "add 2PL"),
  # IO datasets
  list(group = "IO",        variant = "io_slopes",                  label = "BabyView"),
  list(group = "IO",        variant = "io_slopes_seedlings",        label = "SEEDLingS"),
  # Peekbank
  list(group = "Proc",      variant = "long_proc_slopes",           label = "Stanford+LWL"),
  # Difficulty-side ablations (English & Norwegian where available)
  list(group = "Difficulty", variant = "long_no_class_slopes",                label = "no_class (en)"),
  list(group = "Difficulty", variant = "long_no_class_slopes_norwegian",      label = "no_class (nor)"),
  list(group = "Difficulty", variant = "long_class_beta_slopes",              label = "class_beta (en)"),
  list(group = "Difficulty", variant = "long_class_beta_slopes_norwegian",    label = "class_beta (nor)"),
  # LMM (linear-in-age) comparison
  list(group = "LMM",       variant = "long_lmm_slopes",            label = "LMM (en)"),
  list(group = "LMM",       variant = "long_lmm_slopes_norwegian",  label = "LMM (nor)")
)

PARAMS <- c(
  "sigma_alpha", "sigma_xi", "sigma_zeta", "rho_xi_zeta",
  "pi_alpha", "delta", "s", "beta_age", "sigma_lambda",
  "beta_c[1]", "beta_c[2]", "beta_c[3]", "beta_c[4]",
  "beta_react", "sigma_within", "sigma_r",
  "gamma_rt", "mu_rtslope", "sigma_rtslope", "sigma_lwl"
)

summarize_fit <- function(entry) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", entry$variant))
  if (!file.exists(path)) {
    cat(sprintf("  [missing] %s\n", path)); return(NULL)
  }
  fit   <- readRDS(path)
  draws <- as_draws_df(fit)

  rows <- lapply(PARAMS, function(p) {
    if (!p %in% names(draws)) return(NULL)
    x <- draws[[p]]
    # Skip pinned-near-zero parameters (tight prior).
    if (sd(x) < 1e-3) return(NULL)
    tibble(group = entry$group, variant = entry$variant, label = entry$label,
           param = p,
           median = median(x),
           lo95   = quantile(x, 0.025, names = FALSE),
           hi95   = quantile(x, 0.975, names = FALSE),
           rhat   = tryCatch(rstan::Rhat(x), error = function(e) NA_real_))
  })
  bind_rows(rows)
}

cat("Loading fits...\n")
all_scalars <- bind_rows(lapply(FITS, summarize_fit))

# ---- Pretty-print each group as wide table ---- #
print_wide <- function(grp_name) {
  d <- all_scalars %>% filter(group == grp_name)
  if (nrow(d) == 0) { cat("  (no fits in group ", grp_name, ")\n"); return() }
  d <- d %>%
    mutate(cell = sprintf("%.2f [%.2f, %.2f]", median, lo95, hi95)) %>%
    select(label, param, cell) %>%
    pivot_wider(names_from = param, values_from = cell)
  cat("\n=========== ", grp_name, " ===========\n", sep = "")
  print(as.data.frame(d), row.names = FALSE)
}

print_wide("English")
print_wide("Norwegian")
print_wide("IO")
print_wide("Proc")
print_wide("Difficulty")
print_wide("LMM")

# ---- Save CSV for re-use ---- #
out_csv <- file.path(PATHS$figs_dir, "results_summary.csv")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(all_scalars, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
