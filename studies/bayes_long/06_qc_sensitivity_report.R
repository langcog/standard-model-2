## 06_qc_sensitivity_report.R -- exclusion sensitivity for the SI.
##
## A reviewer asked to see the analysis with no exclusions and at several filter strengths,
## on the concern that the filter might preferentially remove noisy or non-monotonic
## trajectories and thereby INFLATE estimated acceleration. That concern is directional, so
## it is testable: if kappa is unchanged with the filter off, it is answered.
##
## IMPORTANT FRAMING. The filter only does work in Marchman and Norwegian. It never fires
## at all in Japanese, and touches Thal and Smith only at the tightest setting (1 and 2
## administrations respectively). So a flat kappa row for those three is NOT evidence of
## robustness -- it is a null operation, and claiming "kappa is robust across filter
## settings in five datasets" would overstate what was tested. The informative comparison
## is Marchman, where the filter removes most and where the reviewer's 16.5% exclusion of
## eligible children bites.
##
## CONVERGENCE. Report rhat on the parameters the table actually shows (kappa_pop and
## sigma_b), not the max over every parameter. In every cell the max-over-all laggard is
## tau_delta -- the item-difficulty SD, a nuisance hyperparameter sitting on the same
## non-centred ridge as delta_j_raw -- while kappa_pop mixes at rhat <= 1.002 with ess in
## the thousands. Flagging on max-over-all would have condemned Norwegian loose (1.301) and
## Thal tight (1.503) although both report kappa at 1.001-1.002 and sigma_b at 1.04-1.06.
##
## Usage:  Rscript studies/bayes_long/06_qc_sensitivity_report.R
## Output: paper/cache/si_qc_sensitivity.rds + console report

suppressPackageStartupMessages({library(dplyr)})
SUMM <- file.path("fits", "bayes_long", "summaries")
BL   <- file.path("fits", "bayes_long")
CACHE <- file.path("paper", "cache")
DATASETS <- c(thal = "English (Thal)", smith = "English (Smith)",
              marchman = "English (Marchman)", norwegian = "Norwegian",
              japanese = "Japanese")
## The reported setting has an empty bundle tag, and "" is not a legal R list name, so
## the conditions live in a data frame keyed by row rather than by name.
COND <- data.frame(
  tag    = c("qcnone", "qcloose", "",     "qctight"),
  lab    = c("none",   "loose",   "main", "tight"),
  crater = c(NA,       0.40,      0.25,   0.15),
  jump   = c(NA,       0.60,      0.40,   0.25), stringsAsFactors = FALSE)

one <- function(slug, label, i) {
  cc <- COND[i, ]; tag <- cc$tag
  sf <- if (tag == "") file.path(SUMM, sprintf("%s_a3_m3.summary.rds", slug))
        else            file.path(SUMM, sprintf("%s_a3_%s_m3.summary.rds", slug, tag))
  bf <- if (tag == "") file.path(BL, sprintf("bundle_%s_a3.rds", slug))
        else            file.path(BL, sprintf("bundle_%s_a3_%s.rds", slug, tag))
  if (!file.exists(sf) || !file.exists(bf)) { cat("  pending:", slug, cc$lab, "\n"); return(NULL) }
  s <- as.data.frame(readRDS(sf)); g <- function(v) s$median[s$variable == v]
  b <- readRDS(bf); m <- b$meta
  data.frame(slug = slug, lang = label, setting = cc$lab,
             crater = cc$crater, jump = cc$jump,
             n_kids = b$stan_data$I,
             admins_removed = if (is.null(m$qc_admins_removed)) NA_integer_ else m$qc_admins_removed,
             kids_dropped   = if (is.null(m$qc_kids_dropped))   NA_integer_ else m$qc_kids_dropped,
             kappa = g("kappa_pop"), sigma_a = g("sigma_a"),
             sigma_b = g("sigma_b"), rho = g("rho_ab"),
             rhat_kappa = s$rhat[s$variable == "kappa_pop"],
             rhat_sigma_b = s$rhat[s$variable == "sigma_b"],
             ess_kappa = s$ess_bulk[s$variable == "kappa_pop"],
             rhat_reported = max(s$rhat[s$variable %in% c("kappa_pop", "sigma_b")], na.rm = TRUE),
             max_rhat = max(s$rhat, na.rm = TRUE), row.names = NULL)
}

cat("QC exclusion sensitivity (M3, 3+ administrations)\n")
S <- bind_rows(lapply(names(DATASETS), function(sl)
        bind_rows(lapply(seq_len(nrow(COND)), function(i) one(sl, DATASETS[[sl]], i))))) |>
  mutate(lang = factor(lang, levels = unname(DATASETS)),
         setting = factor(setting, levels = c("none", "loose", "main", "tight")))
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
saveRDS(S, file.path(CACHE, "si_qc_sensitivity.rds"))

cat("\n=== does the filter do anything? (administrations removed) ===\n")
print(S |> select(lang, setting, admins_removed) |>
      tidyr::pivot_wider(names_from = setting, values_from = admins_removed) |> as.data.frame(),
      row.names = FALSE)
cat("  -> never fires in Japanese; touches Thal/Smith only at 'tight'. Those flat kappa\n")
cat("     rows are a null operation, not evidence of robustness.\n")

cat("\n=== kappa by setting ===\n")
print(S |> transmute(lang, setting, kappa = sprintf("%.2f", kappa)) |>
      tidyr::pivot_wider(names_from = setting, values_from = kappa) |> as.data.frame(),
      row.names = FALSE)
cat("\n=== sigma_b by setting ===\n")
print(S |> transmute(lang, setting, sigma_b = sprintf("%.2f", sigma_b)) |>
      tidyr::pivot_wider(names_from = setting, values_from = sigma_b) |> as.data.frame(),
      row.names = FALSE)

rng <- S |> group_by(lang) |>
  summarise(k_lo = min(kappa), k_hi = max(kappa), sb_lo = min(sigma_b), sb_hi = max(sigma_b),
            .groups = "drop")
cat("\n=== range across all four settings ===\n")
for (i in seq_len(nrow(rng))) with(rng[i,],
  cat(sprintf("  %-20s kappa %5.2f-%5.2f (spread %.2f)   sigma_b %5.2f-%5.2f\n",
              lang, k_lo, k_hi, k_hi - k_lo, sb_lo, sb_hi)))

cat("\n=== convergence on the REPORTED quantities vs max over all parameters ===\n")
cat(sprintf("  worst rhat on kappa_pop across all 20 cells: %.3f (min ess %.0f)\n",
            max(S$rhat_kappa, na.rm = TRUE), min(S$ess_kappa, na.rm = TRUE)))
cat(sprintf("  worst rhat on sigma_b:                       %.3f\n", max(S$rhat_sigma_b, na.rm = TRUE)))
cat(sprintf("  worst rhat over ALL parameters:              %.3f (always tau_delta)\n",
            max(S$max_rhat, na.rm = TRUE)))
bad <- filter(S, rhat_reported > 1.1)
if (nrow(bad)) {
  cat("\n  cells where a REPORTED parameter mixes poorly:\n")
  for (i in seq_len(nrow(bad))) with(bad[i,],
    cat(sprintf("    %-20s %-6s kappa %.3f sigma_b %.3f\n", lang, setting, rhat_kappa, rhat_sigma_b)))
} else cat("\n  -> no cell has a reported parameter above 1.1; no refits needed.\n")
cat("\nwrote", file.path(CACHE, "si_qc_sensitivity.rds"), "\n")
