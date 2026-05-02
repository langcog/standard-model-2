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
//
// Linear predictor:
//   eta = lambda_j * [xi_i + beta_c[cc[j]] * log p_j + log H
//                     + (1 + delta + zeta_i) * log((age_a - s)/a0) - psi_j]
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
                        array[] int admin_to_child, array[] int cc,
                        // global / scalar
                        vector admin_age, real s, real a0,
                        real time_baseline, real delta,
                        real log_H,
                        // per-item / per-child
                        vector log_p,
                        vector psi, vector lambda,
                        vector beta_c,
                        vector xi, vector zeta) {
    int n_slice = end - start + 1;
    vector[n_slice] eta_slice;
    for (i in 1:n_slice) {
      int n = start + i - 1;
      real ae = fmax(admin_age[aa[n]] - s, 0.01);
      real log_age_n = log(ae / a0);
      int ch = admin_to_child[aa[n]];
      real slope_n = time_baseline + delta + zeta[ch];
      real beta_n  = beta_c[cc[jj[n]]];
      real base = xi[ch] + beta_n * log_p[jj[n]] + log_H
                + slope_n * log_age_n - psi[jj[n]];
      eta_slice[i] = lambda[jj[n]] * base;
    }
    return bernoulli_logit_lpmf(y_slice | eta_slice);
  }
}

data {
  int<lower=1> N;                       // observations
  int<lower=1> A;                       // admins
  int<lower=1> I;                       // children
  int<lower=1> J;                       // words
  int<lower=1> C;                       // lexical classes

  array[N] int<lower=1, upper=A> aa;    // admin index per obs
  array[N] int<lower=1, upper=J> jj;    // word index per obs
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;
  array[N] int<lower=0, upper=1> y;

  vector[A] admin_age;                  // age per admin (months)
  vector[J] log_p;

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
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> beta_c_prior_sd;        // tight => beta_c pinned at beta_c_prior_mean
  real beta_c_prior_mean;               // 1 = unit accumulator; 0 = no_freq
  real time_baseline;                   // 1 = unit-rate accumulator; 0 = no time term (M0)
}

parameters {
  // Child-level bivariate (xi, zeta) using non-centered Cholesky
  matrix[2, I] z_child;                 // standardized
  vector<lower=0>[2] sigma_child;       // (sigma_xi_effective, sigma_zeta)
  cholesky_factor_corr[2] L_child;      // Cholesky of correlation

  real<lower=0> sigma_alpha;

  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real<lower=0, upper=15> s;
  real delta;

  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;

  vector[C] beta_c;                     // per-class log-p slope (pinned 1 by default)
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

  vector[J] psi;
  for (j in 1:J) {
    psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];
  }
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
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
  sigma_alpha   ~ normal(0, 1);

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);

  beta_c ~ normal(beta_c_prior_mean, beta_c_prior_sd);

  // Likelihood: parallelize per-observation lpmf via reduce_sum.
  // grainsize = 1 lets Stan's TBB scheduler auto-tune the slice size.
  target += reduce_sum(partial_sum_lpmf, y, 1,
                       aa, jj, admin_to_child, cc,
                       admin_age, s, a0,
                       time_baseline, delta, log_H,
                       log_p, psi, lambda, beta_c,
                       xi, zeta);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real rho_xi_zeta = L_child[2, 1];  // since L_child is 2x2 Cholesky

  // Per-observation log-likelihood for LOO/WAIC across model variants.
  // Recomputes the linear predictor; doesn't slow the sampling phase.
  vector[N] log_lik;
  {
    vector[N] ae;
    for (n in 1:N) ae[n] = fmax(admin_age[aa[n]] - s, 0.01);
    vector[N] log_age = log(ae / a0);
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      real slope_n = time_baseline + delta + zeta[ch];
      real base_n = xi[ch] + beta_c[cc[jj[n]]] * log_p[jj[n]] + log_H
                  + slope_n * log_age[n] - psi[jj[n]];
      real eta_n = lambda[jj[n]] * base_n;
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
    }
  }

  // posterior-expected log_alpha given xi
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
