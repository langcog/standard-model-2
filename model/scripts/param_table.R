## Build a parameter-estimate table across the 5 headline fits:
##   1. EN M_best (long_no_freq_slopes, I=500)
##   2. NO M_best (long_no_freq_slopes_norwegian, I=500)
##   3. BabyView IO  (io_no_freq_slopes)
##   4. SEEDLingS IO (io_no_freq_slopes_seedlings)
##   5. Peekbank/proc (long_proc_no_freq_slopes)
##
## Outputs:
##   outputs/param_table.csv  -- wide table for spreadsheet import
##   outputs/param_table.md   -- markdown table for slide pasting
##   console-printed pretty version

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(knitr)
})

SUMM_DIR <- file.path(PATHS$fits_dir, "summaries")
OUT_DIR  <- PATHS$outputs_dir
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FITS <- list(
  list(col = "EN M_best (I=500)",      tag = "long_no_freq_slopes"),
  list(col = "NO M_best (I=500)",      tag = "long_no_freq_slopes_norwegian"),
  list(col = "BabyView (IO)",          tag = "io_no_freq_slopes"),
  list(col = "SEEDLingS (IO)",         tag = "io_no_freq_slopes_seedlings"),
  list(col = "Peekbank (proc)",        tag = "long_proc_no_freq_slopes")
)

# Parameters to report, in display order. Some only appear in a subset
# of the fits; missing entries become "--".
PARAM_ORDER <- c(
  "I", "J", "N",
  "delta", "sigma_alpha", "sigma_zeta", "pi_alpha", "rho_xi_zeta",
  "mu_r", "sigma_r",
  "gamma_rt"
)

PARAM_LABELS <- c(
  I = "I (kids)", J = "J (items)", N = "N (obs)",
  delta        = "delta",
  sigma_alpha  = "sigma_alpha",
  sigma_zeta   = "sigma_zeta",
  pi_alpha     = "pi_alpha",
  rho_xi_zeta  = "rho(xi, zeta)",
  mu_r         = "mu_r",
  sigma_r      = "sigma_r",
  gamma_rt     = "gamma_rt"
)

# Bundle paths for sample-size lookup (I, J, N from stan_data).
BUNDLES <- list(
  long_no_freq_slopes               = "long_subset_data.rds",
  long_no_freq_slopes_norwegian     = "long_subset_data_nor.rds",
  io_no_freq_slopes                 = "babyview_subset_data.rds",
  io_no_freq_slopes_seedlings       = "seedlings_subset_data.rds",
  long_proc_no_freq_slopes          = "stanford_linked_subset_data.rds"
)

fmt_val <- function(mean, lo, hi, digits = 2) {
  if (is.na(mean)) return("--")
  sprintf("%.*f [%.*f, %.*f]", digits, mean, digits, lo, digits, hi)
}

extract_row <- function(fit) {
  path <- file.path(SUMM_DIR, paste0(fit$tag, ".summary.rds"))
  if (!file.exists(path)) stop("Missing summary: ", path)
  sm <- as.data.frame(readRDS(path))
  # Try to also load the extras from a draws_full.rds (we wrote one for
  # the proc fit with gamma_rt etc.) -- summary.rds doesn't include them.
  extras_path <- file.path(SUMM_DIR, paste0(fit$tag, ".draws_full.rds"))
  if (file.exists(extras_path)) {
    d <- as.data.frame(readRDS(extras_path))
    extra_pars <- setdiff(intersect(PARAM_ORDER, names(d)), sm$variable)
    if (length(extra_pars) > 0) {
      extras <- do.call(rbind, lapply(extra_pars, function(p) {
        v <- d[[p]]
        data.frame(variable = p, mean = mean(v),
                   q025 = quantile(v, 0.025, names = FALSE),
                   q975 = quantile(v, 0.975, names = FALSE))
      }))
      sm <- bind_rows(sm[, c("variable","mean","q025","q975")], extras)
    }
  }
  bundle_path <- file.path(PATHS$fits_dir, BUNDLES[[fit$tag]])
  if (!file.exists(bundle_path)) stop("Missing bundle: ", bundle_path)
  sd_b <- readRDS(bundle_path)$stan_data

  # Build per-param formatted strings.
  vals <- sapply(PARAM_ORDER, function(p) {
    if (p %in% c("I", "J", "N")) return(sprintf("%d", sd_b[[p]]))
    row <- sm[sm$variable == p, , drop = FALSE]
    if (nrow(row) == 0) return("--")
    fmt_val(row$mean[1], row$q025[1], row$q975[1])
  })
  data.frame(
    parameter = PARAM_ORDER,
    label     = PARAM_LABELS[PARAM_ORDER],
    value     = vals,
    column    = fit$col,
    stringsAsFactors = FALSE
  )
}

rows <- bind_rows(lapply(FITS, extract_row))

wide <- rows |>
  select(label, column, value) |>
  pivot_wider(names_from = column, values_from = value)

cat("=== Parameter table (rounded to 2 d.p., 95% CI in brackets) ===\n\n")
print(knitr::kable(wide, format = "pipe", align = "lrrrrr"))

write_csv(wide, file.path(OUT_DIR, "param_table.csv"))
cat(sprintf("\nWrote %s\n", file.path(OUT_DIR, "param_table.csv")))

md_lines <- capture.output(print(knitr::kable(wide, format = "pipe",
                                               align = "lrrrrr")))
writeLines(md_lines, file.path(OUT_DIR, "param_table.md"))
cat(sprintf("Wrote %s\n", file.path(OUT_DIR, "param_table.md")))
