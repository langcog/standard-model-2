## Core helpers for the log-linear IRT accumulator model.
## Loaded by every driver script via source("model/R/helpers.R").
##
## Provides:
##   load_wordbank_data()
##   load_input_rate_prior()
##   subsample_wordbank()
##   build_stan_data()
##   simulate_data()
##   fit_variant()
##   summarize_fit()
##   variant_hyperpriors()
##   plot_delta_j_vs_logp(), plot_class_means(), plot_posterior_density()

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(ggplot2)
  library(tibble)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ===========================================================================
# Data loading
# ===========================================================================

## Load English cross-sectional Wordbank data, drawn from the
## WG-and-WS-combined long_items.rds. One admin per (child, form):
## the most recent admin per child in each form. Children may appear
## twice if they have data on both WG and WS — they're treated as
## separate "persons" for the cross-sectional fit (admin_key as id).
##
## This replaces the older WS-only engWS_preprocessed.Rdata.
load_wordbank_data <- function(
  long_items_path = file.path(PATHS$fits_dir, "long_items.rds"),
  language = "English (American)"
) {
  long <- readRDS(long_items_path)
  long %>%
    filter(language == !!language,
           !is.na(prob), prob > 0,
           !is.na(produces)) %>%
    rename(lexical_class = lexical_category) %>%
    # Treat each (child, form, age) admin as its own "person" for
    # cross-sectional purposes; collapses into a single row per child
    # if they only have one admin.
    mutate(person = paste(child_id, form, age, sep = "_")) %>%
    select(person, age, item, lexical_class, prob, produces, child_id, form)
}

load_input_rate_prior <- function(path = PATHS$input_rate,
                                  column = "adult_child_tokens_hr") {
  hr <- read.csv(path)
  x <- hr[[column]]
  x <- x[!is.na(x)]
  list(mu_r = mean(log(x)), sigma_r = sd(log(x)), n = length(x))
}

# ===========================================================================
# Subsampling + Stan data construction
# ===========================================================================

subsample_wordbank <- function(df, n_children, n_items,
                               seed = 20250420, age_bins = 5,
                               class_col = "lexical_class") {
  set.seed(seed)

  persons <- df %>% distinct(person, age) %>%
    mutate(age_bin = cut(age, breaks = age_bins)) %>%
    group_by(age_bin) %>%
    slice_sample(n = max(1, floor(n_children / age_bins))) %>%
    ungroup() %>%
    pull(person)
  if (length(persons) > n_children) persons <- sample(persons, n_children)

  items <- df %>% distinct(item, !!sym(class_col)) %>%
    group_by(!!sym(class_col)) %>%
    slice_sample(n = max(1, floor(n_items / 5))) %>%
    ungroup() %>%
    pull(item)
  if (length(items) > n_items) items <- sample(items, n_items)

  df %>% filter(person %in% persons, item %in% items)
}

build_stan_data <- function(df, prior_r,
                            priors = DEFAULT_PRIORS,
                            constants = MODEL_CONSTANTS) {
  df <- df %>%
    mutate(ii = as.integer(factor(person)),
           jj = as.integer(factor(item)),
           cc = as.integer(factor(lexical_class)))

  class_levels <- levels(factor(df$lexical_class))
  child_age <- df %>% distinct(ii, age) %>% arrange(ii) %>% pull(age)
  word_info <- df %>% distinct(jj, item, prob, cc) %>% arrange(jj)

  stan_data <- c(
    list(
      N = nrow(df),
      I = max(df$ii), J = max(df$jj), C = max(df$cc),
      ii = df$ii, jj = df$jj, cc = word_info$cc,
      y  = df$produces,
      age = child_age,
      log_p = log(word_info$prob),
      log_H = constants$log_H,
      a0    = constants$a0,
      mu_r = prior_r$mu_r,
      sigma_r = prior_r$sigma_r
    ),
    priors
  )

  list(stan_data = stan_data,
       word_info = word_info,
       class_levels = class_levels,
       child_age = child_age,
       df = df)
}

# Variants opt in to model components ABOVE the lean baseline defined
# in DEFAULT_PRIORS (Rasch + frequency + per-class delta_j + free delta,
# with s pinned at 0 and no per-child slopes).
#
# Naming convention: variants are named for what they ADD.
#   baseline      - lean defaults; the cross-sectional default
#   slopes        - + per-child slopes zeta (default for longitudinal)
#   2pl           - + 2PL item discrimination
#   2pl_slopes    - + both
#   free_s        - frees the start time s (RQ2 robustness)
#   fix_delta     - pins delta = 0 (RQ3 ablation; "no acceleration")
#   no_freq       - drops frequency by setting beta=0 in prep (handled
#                   in data prep, not here; left as a sentinel)
#
# Long-form variants (long_*) prefix is preserved for clarity in the
# longitudinal pipeline but resolves to the same overrides.
variant_hyperpriors <- function(name) {
  # Strip pipeline prefixes ("long_proc_", "long_", "io_") so all
  # pipelines share the same variant grammar. Order matters: longer
  # prefixes first so they match before their substrings.
  base <- sub("^(long_proc_|long_|io_)", "", name)

  # === ACTIVE / HEADLINE VARIANTS ===
  # The five "5-panel build" variants for the main figure are:
  #   demo_pure, demo_alpha, demo_kappa_pop, no_freq_slopes,
  #   no_freq_slopes_si_signed
  # The headline (full model) is no_freq_slopes_si_signed.
  # All other variants below are legacy / robustness / ablation; some are
  # retained for cross-referencing fits on disk and for the §14 nested
  # family analysis (m0..m5).
  switch(base,
    # --- Legacy: spine ablations (no s_i; primarily used in §14 LOO) ---
    baseline      = list(),
    slopes        = list(sigma_zeta_prior_sd = 1),
    `2pl`         = list(sigma_lambda_prior_sd = 1),
    `2pl_slopes`  = list(sigma_lambda_prior_sd = 1,
                         sigma_zeta_prior_sd = 1),
    free_s        = list(s_prior_mean = 4.5, s_prior_sd = 2),
    free_s_slopes = list(s_prior_mean = 4.5, s_prior_sd = 2,
                         sigma_zeta_prior_sd = 1),
    fix_delta     = list(delta_prior_mean = 0, delta_prior_sd = 0.001),
    fix_delta_slopes = list(delta_prior_mean = 0, delta_prior_sd = 0.001,
                            sigma_zeta_prior_sd = 1),
    # --- Legacy: difficulty-side ablations ---
    class_beta        = list(beta_c_prior_sd = 0.5),
    class_beta_slopes = list(beta_c_prior_sd = 0.5,
                             sigma_zeta_prior_sd = 1),
    no_class          = list(),  # data-side override; see variant_data_overrides
    no_class_slopes   = list(sigma_zeta_prior_sd = 1),
    # LMM (linear-in-age) variants -- distinct Stan file log_irt_long_lmm.stan
    lmm               = list(),
    lmm_slopes        = list(sigma_zeta_prior_sd = 1),
    # Nested family canonical labels.  Spine has 7 stages (M0..M6) so
    # that "+ time" and "+ frequency" are factored as separate steps.
    # The fit-tag column shows which fit on disk each label maps to
    # (so M3..M6 reuse existing baseline/slopes/class_beta/2pl fits).
    #
    #   label | tag                           | what it adds
    #   M0    | long_m0                       | minimal IRT (no time, no freq)
    #   M1    | long_m1_time_only             | + time only (delta=0, freq=0)
    #   M2    | long_m1                       | + frequency (delta=0)
    #   M3    | long_baseline                 | + free delta
    #   M4    | long_slopes                   | + per-child slope zeta_i
    #   M5    | long_class_beta_slopes        | + class-specific beta_c
    #   M6    | long_m5                       | + 2PL (lambda_j)
    m0 = list(time_baseline = 0,
              beta_c_prior_mean = 0,
              delta_prior_mean = 0, delta_prior_sd = 0.001),
    m1_time_only = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                        delta_prior_mean = 0, delta_prior_sd = 0.001),
    m1 = list(delta_prior_mean = 0, delta_prior_sd = 0.001),
    m2 = list(),  # = baseline (drops slopes, frees delta)
    m3 = list(sigma_zeta_prior_sd = 1),  # = slopes (lean ref)
    m4 = list(sigma_zeta_prior_sd = 1, beta_c_prior_sd = 0.5),  # + class_beta
    m5 = list(sigma_zeta_prior_sd = 1, beta_c_prior_sd = 0.5,
              sigma_lambda_prior_sd = 1),  # + 2PL
    # === ACTIVE: 5-panel "model spine" variants for the quantile-demo
    # figure. Each progressively frees one parameter on top of the pure
    # accumulator. All are no_freq (beta_c pinned at 0) so delta_j
    # absorbs item difficulty without a separate frequency channel.
    #
    #   variant         | sigma_alpha | sigma_zeta | delta  | what it represents
    #   demo_pure       | pinned ~0   | pinned ~0  | pinned | McMurray pure accumulator
    #   demo_alpha      | free        | pinned ~0  | pinned | + per-kid efficiency (Kachergis 2021)
    #   demo_kappa_pop  | free        | pinned ~0  | free   | + population acceleration (Hidaka rate-change)
    #   no_freq_slopes  | free        | free       | free   | + per-kid acceleration (was M_best)
    #   no_freq_slopes_si_signed | + sigma_s free  |        | + per-kid trajectory phase (this work)
    #
    # Legacy demo_kappa (sigma_alpha pinned, sigma_zeta + delta free) is
    # kept for backward compatibility but superseded by demo_kappa_pop
    # in the 5-panel build.
    demo_pure  = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                      sigma_alpha_prior_sd = 0.001,
                      delta_prior_mean = 0, delta_prior_sd = 0.001),
    demo_alpha = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                      delta_prior_mean = 0, delta_prior_sd = 0.001),
    # demo_kappa_pop: + population scaling kappa_pop free (was implicit in
    # demo_kappa, but demo_kappa also turns on sigma_zeta). This variant
    # frees kappa_pop alone, with sigma_zeta still pinned -- gives the
    # "all kids share the same kappa > 1" panel between demo_alpha and
    # slopes in the 5-panel additive build.
    demo_kappa_pop = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                          sigma_zeta_prior_sd = 0.001),
    demo_kappa = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                      sigma_alpha_prior_sd = 0.001,
                      sigma_zeta_prior_sd = 1),
    demo_full  = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                      sigma_zeta_prior_sd = 1),
    # no_freq variant: drops log p_j entirely (beta_c pinned at 0).
    # With delta and slopes free this is the off-spine "M4 without freq"
    # robustness check; m1_time_only above is the equivalent test on the spine.
    no_freq        = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001),
    # M_best (no_freq + slopes): the headline model after we excised
    # s and s_i (2026-05-22). Rationale: at I=500 the (s, delta) ridge
    # and the fmax(age - s - s_i, 0.01) floor created multi-modal
    # posteriors that did not appear at I=50. The May 3 fit with
    # s_prior_mean=0.5 converged cleanly; we go one step further and
    # pin s at 0 (s_prior_sd=0.001) so it's effectively a constant.
    # s_i remains pinned off via the default sigma_s_prior_sd=0.001.
    # Linear predictor reduces to:
    #   eta = xi + log H + (1 + delta + zeta) * log(age/a0) - delta_j
    no_freq_slopes = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                          sigma_zeta_prior_sd = 1,
                          s_prior_mean = 0, s_prior_sd = 0.001),
    # M_best (no_freq + slopes) with global start time s freed.
    # Used to investigate whether onset-time effects (per-child or
    # population) help with the 4-panel demo "compressed at top, wide
    # at bottom" pattern at very young ages.
    free_s_no_freq_slopes = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                                 sigma_zeta_prior_sd = 1,
                                 s_prior_mean = 4.5, s_prior_sd = 2),
    # Per-child onset variants. Population s pinned tight near 0 (via
    # DEFAULT_PRIORS s_prior_mean=0.5, s_prior_sd=0.05); per-child onset
    # offset s_i ~ Normal_+(0, sigma_s) with sigma_s ~ HN(0, 3) lets
    # individual kids start late by up to ~9 mo. Two variants:
    #   no_freq_si_only -- sigma_zeta pinned at 0 (no kappa variation),
    #                      sigma_s free; tests whether per-child onset
    #                      ALONE can match the panel-4 fan asymmetry.
    #   no_freq_slopes_si -- sigma_zeta free AND sigma_s free; full
    #                      M_best plus per-child onset.
    # The two together let us see whether sigma_s competes with
    # sigma_zeta (the identifiability worry) or adds genuinely new
    # variance (closing the q10-at-floor gap at age 16).
    # sigma_s_prior_sd = 2 (was 3 in the first round): the first si_only
    # fit pushed sigma_s to 6.5 (~ 2.2 prior SDs above mean), which over-
    # predicted the q90 fan at older ages. HN(0, 2) keeps individual
    # onsets in a more empirically plausible range (~0 to 6 mo).
    no_freq_si_only = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                           sigma_zeta_prior_sd = 0.001,
                           sigma_s_prior_sd = 2),
    no_freq_slopes_si = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                             sigma_zeta_prior_sd = 1,
                             sigma_s_prior_sd = 2),
    # Same priors as no_freq_slopes_si, but uses log_irt_long_si_corr.stan
    # via the _si_corr dispatch in fit_longitudinal.R: trivariate
    # Cholesky on (xi, zeta, s_i_lat) with LKJ(2) prior on the 3x3
    # correlation matrix, plus Tobit clipping s_i = fmax(s_i_lat, 0).
    # Tests whether kid-level (xi, zeta, s_i) co-vary.
    no_freq_slopes_si_corr = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                                  sigma_zeta_prior_sd = 1,
                                  sigma_s_prior_sd = 2),
    # Signed-normal s_i variant: combines the (sigma_total, p_zeta)
    # reparam with signed s_i (sum-to-zero centered). Two-pronged fix:
    # the reparam handles the (sigma_zeta, sigma_s) ridge; the signed
    # s_i handles the (sigma_s, delta) population-mean compensation
    # ridge that remained in the half-normal variant. With signed s_i,
    # E[s_i] = 0 by construction regardless of sigma_s, so delta no
    # longer trades off with sigma_s. Interpretation: s_i is a
    # developmental offset (deviation from population reference), not a
    # literal delay -- kids with s_i < 0 are "ahead" of the population
    # reference at any given calendar age. Routes via _si_signed
    # dispatch in fit_longitudinal.R.
    no_freq_slopes_si_signed = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                                    sigma_zeta_prior_sd = 1,
                                    sigma_s_prior_sd = 2),
    # Comprehension-channel variants. Used only by log_irt_long_io_comp.stan;
    # fit_io.R selects that Stan file when the variant matches `comp_*`.
    # The bundle already carries N_comp, aa_comp, jj_comp, y_comp; this
    # override unlocks the gamma_0 / gamma_1 priors so the comp shift
    # is freely estimated. Pair with whatever production-side variant
    # you want (slopes, no_freq_slopes, etc.).
    comp_slopes = list(sigma_zeta_prior_sd = 1,
                       gamma_0_prior_sd = 2, gamma_1_prior_sd = 2),
    comp_no_freq_slopes = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                               sigma_zeta_prior_sd = 1,
                               gamma_0_prior_sd = 2, gamma_1_prior_sd = 2),
    # Standardized-test channel variants. Same Stan file; unlocks
    # gamma_std prior (default tight at 0). The standardized-score
    # observations sigma-anchor sigma_alpha jointly with the CDI side,
    # so the per-child latent log_alpha is identified from two readouts
    # rather than one. Combined `comp_std_*` variants enable both.
    std_slopes = list(sigma_zeta_prior_sd = 1,
                      gamma_std_prior_sd = 2),
    std_no_freq_slopes = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                              sigma_zeta_prior_sd = 1,
                              gamma_std_prior_sd = 2),
    comp_std_no_freq_slopes = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                                   sigma_zeta_prior_sd = 1,
                                   gamma_0_prior_sd = 2, gamma_1_prior_sd = 2,
                                   gamma_std_prior_sd = 2),
    # comp + si_signed: same as comp_no_freq_slopes plus the per-child
    # trajectory phase s_i (signed-normal, sum-to-zero). Used by the
    # log_irt_long_io_comp_si_signed.stan dispatch in fit_io.R for the
    # SEEDLingS s=6 refit.
    comp_no_freq_slopes_si_signed = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                                         sigma_zeta_prior_sd = 1,
                                         gamma_0_prior_sd = 2, gamma_1_prior_sd = 2,
                                         sigma_s_prior_sd = 2),
    comp_std_no_freq_slopes_si_signed = list(beta_c_prior_mean = 0, beta_c_prior_sd = 0.001,
                                             sigma_zeta_prior_sd = 1,
                                             gamma_0_prior_sd = 2, gamma_1_prior_sd = 2,
                                             gamma_std_prior_sd = 2,
                                             sigma_s_prior_sd = 2),
    # Legacy variants for re-loading old fits / explicit comparison
    fix_s         = list(s_prior_mean = 2, s_prior_sd = 0.001),
    both_fixed    = list(delta_prior_mean = 0, delta_prior_sd = 0.001,
                         s_prior_mean = 2,    s_prior_sd = 0.001),
    `2pl_fix_delta` = list(sigma_lambda_prior_sd = 1,
                           delta_prior_mean = 0, delta_prior_sd = 0.001),
    stop(sprintf("Unknown variant: %s", name))
  )
}

# Some variants need to mutate stan_data structure, not just priors.
# Currently `no_class` is the only such variant: it collapses the
# class-hierarchical prior on delta_j into a single global N(mu, tau)
# by overriding cc and C in stan_data. Equivalent to "C = 1 class".
variant_data_overrides <- function(stan_data, variant) {
  base <- sub("^(long_proc_|long_|io_)", "", variant)
  if (base %in% c("no_class", "no_class_slopes")) {
    stan_data$cc <- rep(1L, length(stan_data$cc))
    stan_data$C  <- 1L
    message(sprintf("[variant %s] data override: cc <- 1, C <- 1 (class hierarchy disabled)",
                    variant))
  }
  stan_data
}

# ===========================================================================
# Simulation (for recovery test)
# ===========================================================================

simulate_data <- function(I = 250, J = 150, C = 3,
                          mu_r = log(1198), sigma_r = 0.4,
                          sigma_alpha_true = 0.5,
                          mu_c_true  = c(6.5, 8.0, 9.5),
                          tau_c_true = c(0.5, 0.7, 0.7),
                          s_true = 4.5, delta_true = 0.1,
                          log_H = log(365), a0 = 20,
                          age_range = c(12, 30),
                          log_p_range = c(log(1e-5), log(1e-3)),
                          seed = 42) {
  set.seed(seed)
  stopifnot(length(mu_c_true) == C, length(tau_c_true) == C)

  cc <- sort(rep(seq_len(C), length.out = J))
  cc <- sample(cc)
  log_p <- runif(J, log_p_range[1], log_p_range[2])
  delta_j   <- rnorm(J, mu_c_true[cc], tau_c_true[cc])

  log_r     <- rnorm(I, mu_r, sigma_r)
  log_alpha <- rnorm(I, 0, sigma_alpha_true)
  ages      <- runif(I, age_range[1], age_range[2])

  idx <- expand.grid(ii = seq_len(I), jj = seq_len(J),
                     KEEP.OUT.ATTRS = FALSE)
  ae <- pmax(ages[idx$ii] - s_true, 0.01)
  eta <- log_r[idx$ii] + log_alpha[idx$ii] +
         log_p[idx$jj] + log_H +
         (1 + delta_true) * log(ae / a0) - delta_j[idx$jj]
  y <- rbinom(length(eta), 1, plogis(eta))

  list(
    obs = data.frame(ii = idx$ii, jj = idx$jj,
                     age = ages[idx$ii],
                     log_p = log_p[idx$jj], y = y),
    true = list(log_r = log_r, log_alpha = log_alpha, ages = ages,
                delta_j = delta_j, cc = cc, log_p = log_p,
                sigma_alpha = sigma_alpha_true,
                mu_c = mu_c_true, tau_c = tau_c_true,
                s = s_true, delta = delta_true,
                mu_r = mu_r, sigma_r = sigma_r),
    constants = list(log_H = log_H, a0 = a0, I = I, J = J, C = C)
  )
}

build_stan_data_from_sim <- function(sim, priors = DEFAULT_PRIORS) {
  obs <- sim$obs
  c(
    list(N = nrow(obs), I = sim$constants$I, J = sim$constants$J,
         C = sim$constants$C,
         ii = obs$ii, jj = obs$jj, cc = sim$true$cc,
         y = obs$y, age = sim$true$ages, log_p = sim$true$log_p,
         log_H = sim$constants$log_H, a0 = sim$constants$a0,
         mu_r = sim$true$mu_r, sigma_r = sim$true$sigma_r),
    priors
  )
}

# ===========================================================================
# Fitting
# ===========================================================================

fit_variant <- function(stan_data, tag,
                        cfg = DEFAULT_FIT_CONFIG,
                        model_path = PATHS$stan_model,
                        fits_dir = PATHS$fits_dir,
                        force = FALSE) {
  fit_file <- file.path(fits_dir, sprintf("%s.rds", tag))
  if (!force && file.exists(fit_file)) {
    message(sprintf("[%s] already fit at %s, loading.", tag, fit_file))
    return(readRDS(fit_file))
  }

  message(sprintf("[%s] fitting (chains=%d, iter=%d, warmup=%d)...",
                  tag, cfg$chains, cfg$iter, cfg$warmup))
  t0 <- Sys.time()
  fit <- stan(
    file    = model_path,
    data    = stan_data,
    chains  = cfg$chains,
    iter    = cfg$iter,
    warmup  = cfg$warmup,
    seed    = cfg$seed,
    control = list(adapt_delta = cfg$adapt_delta,
                   max_treedepth = cfg$max_treedepth)
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("[%s] sampling time: %.1f min", tag, dt))
  saveRDS(fit, fit_file)
  fit
}

# cmdstanr-backed version with reduce_sum threading. Compatible Stan
# files have a `partial_sum_lpmf` function in their `functions` block
# called via `reduce_sum`; non-threaded files compile and run fine
# too (threads_per_chain just goes unused).
#
# Output is converted to a `stanfit` via rstan::read_stan_csv(), so
# downstream code (summarize_fit, posterior::as_draws_df, loo, etc.)
# works unchanged.
# Build a list of per-chain init values from a previous fit's
# summary.rds (the small one written by sherlock/extract_summary_table_only.R).
#
# Stan can only init RAW (parameters-block) variables, not derived
# (transformed-parameters or generated-quantities) ones. So this
# function filters to known-safe raw param names. The set varies by
# Stan file; common-safe defaults are sigma_alpha, delta, s,
# sigma_lambda. For the signed-s_i variant (log_irt_long_si_signed.stan)
# we additionally compute the reparameterized raw params (sigma_total,
# p_zeta) from sigma_zeta and sigma_s in the source summary.
#
# Note: things NOT in raw-params and intentionally excluded:
#   - sigma_xi (derived from sigma_alpha + sigma_r)
#   - sigma_zeta (derived from sigma_child[2] or from sigma_total*sqrt(p_zeta))
#   - sigma_s in the signed variant (derived from sigma_total*sqrt(1-p_zeta))
#   - rho_xi_zeta, rho_xi_s, rho_zeta_s (derived from L_child Cholesky)
#   - pi_alpha (generated quantity from sigma_alpha)
#
# Args:
#   source_tag        - tag of a prior fit (without .summary.rds extension).
#   n_chains          - number of chains to build inits for.
#   jitter_frac       - per-param jitter scale, fraction of source posterior SD.
#   target_variant    - optional variant name to choose appropriate raw params.
#                       Currently supports "signed" (adds sigma_total + p_zeta)
#                       or "default" (sigma_alpha, delta, s, sigma_lambda,
#                       sigma_s if present). Default auto-detects from tag.
build_init_from_summary <- function(source_tag,
                                    n_chains = 4,
                                    jitter_frac = 0.05,
                                    target_variant = "auto",
                                    summaries_dir = NULL) {
  if (is.null(summaries_dir)) {
    candidates <- c(
      file.path(PATHS$fits_dir, "summaries", paste0(source_tag, ".summary.rds")),
      file.path("fits", "summaries", paste0(source_tag, ".summary.rds"))
    )
    summary_path <- candidates[file.exists(candidates)][1]
  } else {
    summary_path <- file.path(summaries_dir, paste0(source_tag, ".summary.rds"))
    candidates <- summary_path
  }
  if (is.na(summary_path) || !file.exists(summary_path)) {
    stop(sprintf("No summary file found for source_tag '%s' (checked: %s)",
                 source_tag, paste(candidates, collapse = ", ")))
  }
  summ <- readRDS(summary_path)

  # Auto-detect target variant from source_tag name.
  if (target_variant == "auto") {
    target_variant <- if (grepl("_si_signed$", source_tag)) "signed" else "default"
  }

  # Raw params Stan can accept inits for, by variant type.
  # In log_irt_long.stan, sigma_zeta is DERIVED from sigma_child[2] (a
  # vector raw param) -- so for the default variant we init sigma_child
  # instead, with element [2] = source sigma_zeta. In signed variant
  # the Cholesky structure is different so we use sigma_total/p_zeta.
  raw_params <- switch(target_variant,
    signed  = c("sigma_alpha", "delta", "s", "sigma_lambda"),
    default = c("sigma_alpha", "delta", "s", "sigma_lambda", "sigma_s"),
    stop("Unknown target_variant: ", target_variant)
  )

  # Cap SD floor to avoid zero-jitter on tight posteriors.
  med_lookup <- setNames(summ$median, summ$variable)
  sd_lookup  <- setNames(summ$sd, summ$variable)

  # For signed variant: reconstruct sigma_total and p_zeta from
  # sigma_zeta and sigma_s in the source. Both quantities must exist
  # in the summary.
  signed_extras <- NULL
  if (target_variant == "signed") {
    if (all(c("sigma_zeta", "sigma_s") %in% summ$variable)) {
      sz <- as.numeric(med_lookup["sigma_zeta"])
      ss <- as.numeric(med_lookup["sigma_s"])
      sigma_total_med <- sqrt(sz^2 + ss^2)
      p_zeta_med      <- sz^2 / sigma_total_med^2
      # Derive SDs via first-order Taylor (rough, just for jitter scale)
      sz_sd <- max(as.numeric(sd_lookup["sigma_zeta"]), 1e-3)
      ss_sd <- max(as.numeric(sd_lookup["sigma_s"]),    1e-3)
      sigma_total_sd <- sqrt((sz * sz_sd)^2 + (ss * ss_sd)^2) / sigma_total_med
      p_zeta_sd      <- 2 * sz * sz_sd / sigma_total_med^2  # rough
      signed_extras <- list(
        sigma_total = list(med = sigma_total_med, sd = max(sigma_total_sd, 1e-3)),
        p_zeta      = list(med = p_zeta_med,      sd = max(p_zeta_sd, 1e-3))
      )
      message(sprintf("  signed-variant init: derived sigma_total=%.3f, p_zeta=%.3f",
                      sigma_total_med, p_zeta_med))
    } else {
      message("  signed-variant requested but source summary lacks sigma_zeta/sigma_s; ",
              "skipping reparam-raw inits")
    }
  }

  # Filter raw_params to those actually present in source.
  vars <- intersect(raw_params, summ$variable)

  # Known >=0 params to clip after jitter (avoid Stan boundary errors).
  positive_params <- c("sigma_alpha", "sigma_zeta", "sigma_s",
                       "sigma_total", "sigma_lambda", "p_zeta", "s",
                       "tau_c")
  upper_one_params <- c("p_zeta")
  s_upper <- 14.99  # leave headroom from the upper=15 bound

  set.seed(20260515)
  inits <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    chain_init <- list()
    for (p in vars) {
      v <- as.numeric(med_lookup[p]) +
           rnorm(1, 0, jitter_frac * max(as.numeric(sd_lookup[p]), 1e-3))
      if (p %in% positive_params) v <- max(v, 1e-4)
      if (p %in% upper_one_params) v <- min(v, 0.999)
      if (p == "s") v <- min(v, s_upper)
      chain_init[[p]] <- v
    }
    # Inject signed-variant reparam params if applicable
    if (!is.null(signed_extras)) {
      for (p in names(signed_extras)) {
        v <- signed_extras[[p]]$med +
             rnorm(1, 0, jitter_frac * signed_extras[[p]]$sd)
        if (p %in% positive_params) v <- max(v, 1e-4)
        if (p %in% upper_one_params) v <- min(v, 0.999)
        chain_init[[p]] <- v
      }
    }
    # For default variant: log_irt_long.stan declares sigma_child as a
    # 2-vector raw param. sigma_zeta is derived as sigma_child[2].
    # If the source has sigma_zeta, build a sigma_child init.
    if (target_variant == "default" && "sigma_zeta" %in% summ$variable) {
      sz <- as.numeric(med_lookup["sigma_zeta"])
      sz_sd <- max(as.numeric(sd_lookup["sigma_zeta"]), 1e-3)
      sz_init <- max(sz + rnorm(1, 0, jitter_frac * sz_sd), 1e-4)
      # sigma_child[1] gets internally replaced by sigma_xi; init at 1
      # for the prior, harmless.
      chain_init$sigma_child <- c(1.0, sz_init)
    }
    inits[[ch]] <- chain_init
  }
  message(sprintf("Built init from %s: %d raw params, %d chains, jitter=%.2f*SD",
                  source_tag, length(inits[[1]]), n_chains, jitter_frac))
  inits
}

fit_variant_cmdstanr <- function(stan_data, tag,
                                  cfg = DEFAULT_FIT_CONFIG,
                                  model_path = PATHS$stan_model,
                                  fits_dir = PATHS$fits_dir,
                                  threads_per_chain = NULL,
                                  init = NULL,
                                  force = FALSE) {
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("cmdstanr not installed; install via install.packages(\"cmdstanr\", repos = \"https://stan-dev.r-universe.dev\") and cmdstanr::install_cmdstan().")

  fit_file <- file.path(fits_dir, sprintf("%s.rds", tag))
  if (!force && file.exists(fit_file)) {
    message(sprintf("[%s] already fit at %s, loading.", tag, fit_file))
    return(readRDS(fit_file))
  }

  if (is.null(threads_per_chain)) {
    n_cores <- as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN",
                                      unset = max(1, parallel::detectCores() %/% cfg$chains)))
    threads_per_chain <- n_cores
  }

  message(sprintf("[%s] cmdstanr fitting (chains=%d, iter=%d, warmup=%d, threads_per_chain=%d)...",
                  tag, cfg$chains, cfg$iter, cfg$warmup, threads_per_chain))

  # Stan's TBB build needs TBB_CXX_TYPE set explicitly when the
  # compiler can't be auto-detected (Sherlock case). cmdstanr will
  # rebuild the binary whenever the .stan source is newer than the
  # cached .rds, and the compile env doesn't inherit TBB_CXX_TYPE
  # unless we set it here. Match what setup_R.R does on install.
  if (!nzchar(Sys.getenv("TBB_CXX_TYPE"))) {
    Sys.setenv(TBB_CXX_TYPE = "gcc")
  }

  m <- cmdstanr::cmdstan_model(model_path,
                                cpp_options = list(stan_threads = TRUE))

  # Persistent CSV output directory. By default cmdstanr writes CSVs
  # to a per-session tempdir that gets auto-cleaned when R exits, so if
  # save_object() crashes (e.g., disk-full during fread of the log_lik
  # column) the draws are lost. Pointing output_dir at a real directory
  # under fits/ means the CSVs survive any post-sampling failure and
  # scalar params can be recovered via cmdstanr::as_cmdstan_fit().
  csv_dir <- file.path(fits_dir, sprintf("csvs_%s", tag))
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  t0 <- Sys.time()
  sample_args <- list(
    data              = stan_data,
    chains            = cfg$chains,
    parallel_chains   = cfg$chains,
    threads_per_chain = threads_per_chain,
    iter_sampling     = cfg$iter - cfg$warmup,
    iter_warmup       = cfg$warmup,
    seed              = cfg$seed,
    adapt_delta       = cfg$adapt_delta,
    max_treedepth     = cfg$max_treedepth,
    refresh           = max(1L, (cfg$iter %/% 10L)),
    output_dir        = csv_dir
  )
  # Optional init values: if supplied (typically a list of length
  # n_chains, each a named list of parameter starting values), pass to
  # cmdstanr. Missing params are drawn from the prior, which is the
  # standard cmdstanr fallback.
  if (!is.null(init)) sample_args$init <- init
  csm <- do.call(m$sample, sample_args)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("[%s] cmdstanr sampling time: %.1f min", tag, dt))

  # Save the CmdStanMCMC object directly. cmdstanr's save_object()
  # reads all the CSV draws into memory before serializing, so the
  # resulting .rds is portable across machines and self-contained.
  # Downstream code uses posterior::as_draws_df(fit), which works on
  # both stanfit and CmdStanMCMC.
  # Wrapped in try() so a save_object failure (typically OOM/disk-full
  # during the internal read_cmdstan_csv) doesn't kill the script
  # before we've at least preserved the CSVs in csv_dir.
  saved <- try(csm$save_object(file = fit_file), silent = FALSE)
  if (inherits(saved, "try-error")) {
    message(sprintf("[%s] save_object FAILED -- CSVs preserved at %s",
                    tag, csv_dir))
    message("Recover via: cmdstanr::as_cmdstan_fit(list.files('",
            csv_dir, "', pattern = 'csv$', full.names = TRUE))")
  }
  csm
}

# ===========================================================================
# Summaries
# ===========================================================================

# Backend-agnostic summary. Works on both rstan stanfit and cmdstanr
# CmdStanMCMC via posterior::summarise_draws(); each row is one
# parameter, columns are mean / lo95 / median / hi95 / n_eff / Rhat.
summarize_fit <- function(fit, pars = c("sigma_alpha", "s", "delta",
                                        "pi_alpha", "sigma_xi")) {
  d <- posterior::as_draws_df(fit)
  pars_present <- intersect(pars, names(d))
  if (length(pars_present) == 0)
    return(tibble(param = character(), mean = double(),
                  lo95 = double(), median = double(), hi95 = double(),
                  n_eff = double(), Rhat = double()))
  s <- posterior::summarise_draws(
    posterior::subset_draws(d, variable = pars_present),
    "mean", q025 = ~ stats::quantile(.x, 0.025, names = FALSE),
    "median",
    q975 = ~ stats::quantile(.x, 0.975, names = FALSE),
    "ess_bulk", "rhat"
  )
  tibble(
    param  = s$variable,
    mean   = s$mean,
    lo95   = s$q025,
    median = s$median,
    hi95   = s$q975,
    n_eff  = s$ess_bulk,
    Rhat   = s$rhat
  )
}

class_threshold_table <- function(fit, class_levels) {
  C <- length(class_levels)
  d <- posterior::as_draws_df(fit)
  q <- function(x, p) stats::quantile(x, p, names = FALSE)
  tibble(
    class       = class_levels,
    mu_median   = sapply(seq_len(C), function(c) median(d[[sprintf("mu_c[%d]",  c)]])),
    mu_lo       = sapply(seq_len(C), function(c) q(d[[sprintf("mu_c[%d]",  c)]], 0.025)),
    mu_hi       = sapply(seq_len(C), function(c) q(d[[sprintf("mu_c[%d]",  c)]], 0.975)),
    tau_median  = sapply(seq_len(C), function(c) median(d[[sprintf("tau_c[%d]", c)]])),
    tau_lo      = sapply(seq_len(C), function(c) q(d[[sprintf("tau_c[%d]", c)]], 0.025)),
    tau_hi      = sapply(seq_len(C), function(c) q(d[[sprintf("tau_c[%d]", c)]], 0.975))
  )
}

extract_delta_j_df <- function(fit, word_info, class_levels) {
  draws <- as_draws_df(fit)
  delta_j_cols <- grep("^delta_j\\[", names(draws), value = TRUE)
  word_info %>%
    mutate(delta_j_median = sapply(delta_j_cols, function(p) median(draws[[p]])),
           delta_j_lo = sapply(delta_j_cols, function(p) quantile(draws[[p]], .025)),
           delta_j_hi = sapply(delta_j_cols, function(p) quantile(draws[[p]], .975)),
           class = class_levels[cc],
           log_p = log(prob))
}

# ===========================================================================
# Plots
# ===========================================================================

plot_delta_j_vs_logp <- function(delta_j_df, save_path = NULL, tag = "") {
  r <- cor(delta_j_df$delta_j_median, delta_j_df$log_p)
  p <- ggplot(delta_j_df, aes(log_p, delta_j_median, colour = class)) +
    geom_errorbar(aes(ymin = delta_j_lo, ymax = delta_j_hi), width = 0, alpha = 0.3) +
    geom_point(size = 1.6, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey40",
                aes(group = 1), linewidth = 0.5, linetype = "dashed") +
    labs(title = sprintf("RQ1: delta_j vs log p%s",
                         if (nzchar(tag)) sprintf(" (%s)", tag) else ""),
         subtitle = sprintf("r = %.2f, R^2 = %.2f", r, r^2),
         x = "log p_j", y = "Posterior median delta_j") +
    theme_minimal(base_size = 12)
  if (!is.null(save_path)) ggsave(save_path, p, width = 7, height = 5, dpi = 150)
  invisible(p)
}

plot_class_means <- function(class_tbl, save_path = NULL, tag = "") {
  p <- ggplot(class_tbl, aes(reorder(class, mu_median), mu_median)) +
    geom_pointrange(aes(ymin = mu_lo, ymax = mu_hi)) +
    coord_flip() +
    labs(title = sprintf("Class thresholds%s",
                         if (nzchar(tag)) sprintf(" (%s)", tag) else ""),
         x = NULL, y = "Log-threshold (log-tokens)") +
    theme_minimal(base_size = 12)
  if (!is.null(save_path)) ggsave(save_path, p, width = 6, height = 3.5, dpi = 150)
  invisible(p)
}

plot_posterior_density <- function(fit, param, save_path = NULL,
                                   xlab = param, tag = "",
                                   vline = NULL, xlim = NULL) {
  draws <- as_draws_df(fit)
  p <- ggplot(tibble(x = draws[[param]]), aes(x)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    labs(title = sprintf("%s%s", param,
                         if (nzchar(tag)) sprintf(" (%s)", tag) else ""),
         x = xlab, y = "density") +
    theme_minimal(base_size = 12)
  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed")
  if (!is.null(xlim)) p <- p + xlim(xlim)
  if (!is.null(save_path)) ggsave(save_path, p, width = 6, height = 4, dpi = 150)
  invisible(p)
}
