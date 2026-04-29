## Centralized paths and model configuration.
## Every driver script starts with: source("model/R/config.R")

# Project root: prefer env-var override (useful on Sherlock); otherwise
# auto-detect from the working directory; fall back to the original dev
# path as a last resort.
.default_root <- "/Users/mcfrank/Projects/standard_model_2"
PROJECT_ROOT <- Sys.getenv("STANDARD_MODEL_ROOT", unset = "")
if (!nzchar(PROJECT_ROOT)) {
  # Look for the Makefile as a landmark in the cwd and parents
  cwd <- getwd()
  while (cwd != "/" && !file.exists(file.path(cwd, "Makefile"))) {
    cwd <- dirname(cwd)
  }
  PROJECT_ROOT <- if (file.exists(file.path(cwd, "Makefile"))) cwd
                  else .default_root
}

PATHS <- list(
  # Inputs (from the original standard_model codebase)
  wordbank = file.path(PROJECT_ROOT,
                      "standard_model/scripts/data/engWS_preprocessed.Rdata"),
  input_rate = file.path(PROJECT_ROOT,
                      "standard_model/scripts/data/hourly_tokens_Sperry_HartRisley.csv"),

  # Stan model
  stan_model = file.path(PROJECT_ROOT, "model/stan/log_irt.stan"),

  # Outputs — may be redirected via STANDARD_MODEL_FITS_DIR / _FIGS_DIR
  fits_dir = Sys.getenv("STANDARD_MODEL_FITS_DIR",
                        unset = file.path(PROJECT_ROOT, "model/fits")),
  figs_dir = Sys.getenv("STANDARD_MODEL_FIGS_DIR",
                        unset = file.path(PROJECT_ROOT, "model/figs")),
  notes_dir = file.path(PROJECT_ROOT, "notes")
)

for (d in c(PATHS$fits_dir, PATHS$figs_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## Model constants
MODEL_CONSTANTS <- list(
  log_H = log(365),   # 12 waking hr/day * 30.44 days/mo
  a0    = 20          # reference age in months
)

## Default priors define the LEAN BASELINE model:
##   * Rasch (no 2PL discrimination)
##   * No per-child slopes (added explicitly via 'slopes' variant for
##     longitudinal data, where they're identifiable)
##   * Start time s pinned at 0 (rarely identified well; opt in via
##     'free_s' for the RQ2 robustness check)
##   * delta (population acceleration) free — it's load-bearing
##   * Frequency enters with unit coefficient on log p_j (separate from psi_j)
##
## Variants OPT IN to extra parameters via variant_hyperpriors() in
## helpers.R. This inverts the older default (which started with the
## full model and let variants pin parameters off).
DEFAULT_PRIORS <- list(
  mu_mu_c          = 8,
  sigma_mu_c       = 3,
  s_prior_mean     = 0,
  s_prior_sd       = 0.001,        # s pinned at 0
  delta_prior_mean = 0,
  delta_prior_sd   = 0.5,
  sigma_lambda_prior_sd = 0.001,   # no 2PL by default
  sigma_zeta_prior_sd   = 0.001    # no slopes by default; opt in for longitudinal
)

## Defaults for fitting.
## (iter=1500, warmup=750, adapt_delta=0.9 matches the 37-min baseline fit;
##  adapt_delta=0.95 with iter=2000 was ~30x slower, so reverted.)
DEFAULT_FIT_CONFIG <- list(
  chains      = 4,
  iter        = 1500,
  warmup      = 750,
  adapt_delta = 0.9,
  max_treedepth = 10,
  seed        = 20250420
)
