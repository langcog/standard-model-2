// Longitudinal log-linear IRT accumulator model.
//
// Extends log_irt.stan to properly handle longitudinal data:
// - Observations are indexed by (admin, item); each admin has its own age.
// - Children (I) are indexed by admin (A), via admin_to_child.
// - Per-child latents xi_i (intercept) and zeta_i (slope deviation) are
//   modeled with a BIVARIATE prior (LKJ on correlation), since
//   longitudinal LMM shows corr(xi, zeta) ~ 0.7.
// - 2PL discrimination lambda_j per word, hierarchical.
// - Global start time s, age rate-change exponent delta.
// - Optional per-child onset offset s_i (>= 0) drawn from a half-normal
//   N_+(0, sigma_s); enabled by sigma_s_prior_sd > 0, pinned off otherwise.
//   Effective onset for child i is (s + s_i[i]); when sigma_s_prior_sd is
//   tiny (default 0.001), s_i collapses to ~0 and the model reduces to
//   the standard global-s formulation.
//
// Linear predictor:
//   eta = lambda_j * [xi_i + beta_c[cc[j]] * log p_j + log H
//                     + (1 + delta + zeta_i) * log((age_a - s - s_i)/a0)
//                     - delta_j]
//
// beta_c is a per-class log-p slope that is pinned at 1 by default
// (tight prior beta_c_prior_sd ~ 0.001 in DEFAULT_PRIORS). The
// `class_beta` variant frees it (e.g., 0.5) to test whether frequency
// enters with class-specific weight.
//
// reduce_sum: the per-observation likelihood is parallelized across
// threads via reduce_sum. Each thread computes eta locally for its
// slice of y, then evaluates bernoulli_logit_lpmf. Linear speedup with
// `STAN_NUM_THREADS` (set in the cmdstanr sample call). The
// per-observation Stan code is duplicated inside `partial_sum_lpmf`;
// the unused `eta` local in the model block is retained for
// generated_quantities (log_lik) which is single-threaded.

functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        // observation indices
                        array[] int aa, array[] int jj,
                        array[] int admin_to_child,
                        // global / scalar
                        vector admin_age, real s, real a0,
                        real time_baseline, real delta,
                        real log_H,
                        // per-item / per-child
                        vector delta_j,
                        vector xi, vector zeta, vector s_i) {
    int n_slice = end - start + 1;
    vector[n_slice] eta_slice;
    for (i in 1:n_slice) {
      int n = start + i - 1;
      int ch = admin_to_child[aa[n]];
      real ae = fmax(admin_age[aa[n]] - s - s_i[ch], 0.01);
      real log_age_n = log(ae / a0);
      real slope_n = time_baseline + delta + zeta[ch];
      eta_slice[i] = xi[ch] + log_H + slope_n * log_age_n - delta_j[jj[n]];
    }
    return bernoulli_logit_lpmf(y_slice | eta_slice);
  }
}

data {
  int<lower=1> N;                       // observations
  int<lower=1> A;                       // admins
  int<lower=1> I;                       // children
  int<lower=1> J;                       // words
  int<lower=1> C;                       // lexical classes (unused in this Stan file
                                        // after the 2026-05-21 cleanup; kept so
                                        // bundles built for earlier versions
                                        // remain accepted)

  array[N] int<lower=1, upper=A> aa;    // admin index per obs
  array[N] int<lower=1, upper=J> jj;    // word index per obs
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;    // unused (see above)
  array[N] int<lower=0, upper=1> y;

  vector[A] admin_age;                  // age per admin (months)
  vector[J] log_p;                      // unused after cleanup; kept for bundle compat

  real log_H;
  real<lower=0> a0;

  real mu_r;
  real<lower=0> sigma_r;

  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  // Toggles (near-zero SD disables feature)
  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;  // unused; kept for bundle compat
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> sigma_alpha_prior_sd;   // tight => sigma_alpha pinned at 0 (pure-accumulator demo variants)
  real<lower=0> sigma_s_prior_sd;       // per-child onset s_i; tight => off
  real<lower=0> beta_c_prior_sd;        // unused; kept for bundle compat
  real beta_c_prior_mean;               // unused; kept for bundle compat
  real time_baseline;                   // 1 = unit-rate accumulator; 0 = no time term (M0)
}

parameters {
  // Child-level bivariate (xi, zeta) using non-centered Cholesky
  matrix[2, I] z_child;                 // standardized
  vector<lower=0>[2] sigma_child;       // (sigma_xi_effective, sigma_zeta)
  cholesky_factor_corr[2] L_child;      // Cholesky of correlation

  real<lower=0> sigma_alpha;

  vector[J] delta_j_raw;
  // FLAT (not class-stratified) hyperprior on delta_j. Replaces the
  // earlier vector[C] mu_c + vector<lower=0>[C] tau_c hierarchy on
  // 2026-05-21 along with the s-prior tightening: we don't use the
  // per-class delta_j shrinkage anywhere downstream (the headline
  // params sigma_alpha, kappa_pop, sigma_zeta, sigma_s don't depend
  // on which class shrinkage delta_j gets), and a single global hyperprior
  // is fewer parameters and simpler to explain. cc[] and C remain in
  // the data block since beta_c uses them.
  real mu_delta;
  real<lower=0> tau_delta;

  real<lower=0, upper=15> s;
  real delta;

  // 2PL discrimination (log_lambda, sigma_lambda) and per-class log-p
  // frequency slope (beta_c) removed 2026-05-21. Both were pinned off
  // in every active variant; removing them simplifies the linear
  // predictor to:  eta = xi + log_H + (1 + delta + zeta) * log_age - delta_j .

  // Per-child onset offset, drawn from N_+(0, sigma_s). Lower-bounded
  // at 0 so kid effective onsets s + s_i[ch] >= s (no negative onsets).
  // CENTERED parameterization: we briefly tried non-centered
  // (s_i = sigma_s * s_i_raw) but that exposed/created a bimodal
  // posterior where one mode had sigma_s ~ 36 with sigma_alpha ~ 0.9
  // -- the multiplicative factorization lets the model trade kid-onset
  // variation for kid-efficiency variation in geometrically distant
  // configurations. Centered ties each s_i directly to the data, which
  // empirically gives single-mode posteriors. When sigma_s_prior_sd is
  // tight (~0.001), sigma_s collapses to ~0 and s_i ~ 0.
  vector<lower=0>[I] s_i;
  real<lower=0> sigma_s;
}

transformed parameters {
  // Enforce sigma_child[1]^2 = sigma_r^2 + sigma_alpha^2 so mu_r and sigma_r
  // pinning remain meaningful. We do this via a deterministic constraint:
  // don't declare sigma_xi as a parameter; compute it.
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  // Build child-level effects with the required marginal SDs.
  // Note: sigma_child[1] is NOT used as a free parameter for the xi SD
  // (we replace it with sigma_xi below), but we keep sigma_child[2] as
  // sigma_zeta freely estimated.
  matrix[I, 2] child_effs;
  {
    matrix[2, 2] L_scaled;
    L_scaled[1, 1] = sigma_xi * L_child[1, 1];
    L_scaled[1, 2] = 0;
    L_scaled[2, 1] = sigma_child[2] * L_child[2, 1];
    L_scaled[2, 2] = sigma_child[2] * L_child[2, 2];
    child_effs = (L_scaled * z_child)';
  }
  // Sum-to-zero centering on each random-effect column. Without this
  // the (delta, mean(zeta)) split is partially unidentified: the
  // likelihood only sees (1 + delta + zeta_i) so the random-effect
  // mean can absorb part of the population slope. Centering forces
  // delta to carry the full population slope and zeta_i to be
  // deviations centered at 0. Same for xi: its mean should equal mu_r.
  vector[I] xi   = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta = child_effs[, 2] - mean(child_effs[, 2]);
  real<lower=0> sigma_zeta = sigma_child[2];

  // Flat (non-class-stratified) delta_j hierarchy.
  // (We briefly tried sum-to-zero centering here on 2026-05-21 to
  // match the kid-level random effects; the 50-kid pilot showed worse
  // chain agreement on rho_xi_zeta and sigma_alpha, so the centering
  // was reverted. The mean(delta_j_unc) subtraction introduces a
  // J-wide coupling between every delta_j_raw[k] that NUTS finds harder
  // to traverse than the small mu_delta-vs-globals soft correlation
  // it was meant to fix.)
  vector[J] delta_j = mu_delta + tau_delta * delta_j_raw;
}

model {
  // Priors
  to_vector(z_child) ~ std_normal();
  // sigma_child[1] is replaced by sigma_xi (deterministic; see
  // transformed parameters), so its raw N(0,1) prior is harmless.
  // sigma_child[2] = sigma_zeta is the actual slopes-toggle param;
  // route the variant-grammar prior here so `baseline` (tight prior
  // ~0.001) actually pins zeta off, and `slopes` (~1) frees it.
  sigma_child[1] ~ normal(0, 1);
  sigma_child[2] ~ normal(0, sigma_zeta_prior_sd);
  L_child       ~ lkj_corr_cholesky(2);  // mild prior toward 0 corr
  sigma_alpha   ~ normal(0, sigma_alpha_prior_sd);

  delta_j_raw ~ std_normal();
  mu_delta  ~ normal(mu_mu_c, sigma_mu_c);   // reuse the old hyperprior numerics
  tau_delta ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  // Per-child onset offset (half-normal via <lower=0> on s_i).
  // sigma_s_prior_sd ~ 0.001 effectively pins s_i at 0 (legacy global-s
  // behavior); sd = 2 frees per-child onset variation. Centered
  // parameterization (data directly informs each s_i) -- non-centered
  // was tried and produced multimodal posteriors; centered is single-mode.
  s_i     ~ normal(0, sigma_s);
  sigma_s ~ normal(0, sigma_s_prior_sd);

  // Likelihood: parallelize per-observation lpmf via reduce_sum.
  // grainsize = 1 lets Stan's TBB scheduler auto-tune the slice size.
  target += reduce_sum(partial_sum_lpmf, y, 1,
                       aa, jj, admin_to_child,
                       admin_age, s, a0,
                       time_baseline, delta, log_H,
                       delta_j,
                       xi, zeta, s_i);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real rho_xi_zeta = L_child[2, 1];  // since L_child is 2x2 Cholesky

  // Per-observation log-likelihood for LOO/WAIC across model variants.
  // Recomputes the linear predictor; doesn't slow the sampling phase.
  vector[N] log_lik;
  {
    vector[N] log_age;
    for (n in 1:N) {
      int ch_n = admin_to_child[aa[n]];
      log_age[n] = log(fmax(admin_age[aa[n]] - s - s_i[ch_n], 0.01) / a0);
    }
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      real slope_n = time_baseline + delta + zeta[ch];
      real eta_n = xi[ch] + log_H + slope_n * log_age[n] - delta_j[jj[n]];
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
    }
  }

  // posterior-expected log_alpha given xi
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
