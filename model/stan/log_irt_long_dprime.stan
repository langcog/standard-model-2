// Longitudinal log-linear IRT accumulator model -- D' (input-on-slope).
//
// Variant of log_irt_long.stan (M_best / no_freq_slopes) that adds an
// INPUT-ON-SLOPE channel to test whether the per-child acceleration is
// (partly) driven by the imputed input rate.
//
// Two changes relative to log_irt_long.stan:
//
//  (1) rho_xi_zeta is PINNED to 0. In M_best the intercept-slope
//      coupling is a free correlation rho on the bivariate (xi, zeta).
//      In D' that channel is removed (xi and zeta are independent), so
//      the ONLY intercept-slope coupling runs through the input term
//      below. This is what makes gamma interpretable: with rho free,
//      gamma and rho would be confounded.
//
//  (2) The slope gains a term gamma * log_r_dev_i, where log_r_dev_i is
//      the IMPUTED input deviation of child i -- the sigma_r share of
//      the intercept deviation:
//          log_r_dev_i = (sigma_r^2 / sigma_xi^2) * (xi_i - mu_r),
//      with sigma_xi^2 = sigma_r^2 + sigma_alpha^2 (the same external
//      sigma_r pin that identifies pi_alpha). Because zeta _|_ xi,
//          Cov(xi_i, kappa_i) = gamma * Var(log_r_dev_i) = gamma * sigma_r^2,
//      i.e. gamma = Cov(xi, kappa) / sigma_r^2. gamma is therefore
//      sensitive to the external sigma_r pin -- run the sigma_r sweep
//      (STAN_SIGMA_R_OVERRIDE) to map that dependence.
//
// Linear predictor:
//   eta = xi_i + log H
//         + (1 + delta + zeta_i + gamma * log_r_dev_i) * log((age_a - s - s_i)/a0)
//         - delta_j
//
// s and s_i remain in the model but are pinned off by the no_freq_slopes
// priors (s_prior_sd ~ 0.001, sigma_s_prior_sd ~ 0.001), exactly as in
// log_irt_long.stan. reduce_sum threading is unchanged.

functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        // observation indices
                        array[] int aa, array[] int jj,
                        array[] int admin_to_child,
                        // global / scalar
                        vector admin_age, real s, real a0,
                        real time_baseline, real delta,
                        real log_H, real gamma_in,
                        // per-item / per-child
                        vector delta_j,
                        vector xi, vector zeta, vector s_i,
                        vector log_r_dev) {
    int n_slice = end - start + 1;
    vector[n_slice] eta_slice;
    for (i in 1:n_slice) {
      int n = start + i - 1;
      int ch = admin_to_child[aa[n]];
      real ae = fmax(admin_age[aa[n]] - s - s_i[ch], 0.01);
      real log_age_n = log(ae / a0);
      real slope_n = time_baseline + delta + zeta[ch] + gamma_in * log_r_dev[ch];
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
  int<lower=1> C;                       // lexical classes (unused; bundle compat)

  array[N] int<lower=1, upper=A> aa;    // admin index per obs
  array[N] int<lower=1, upper=J> jj;    // word index per obs
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;    // unused
  array[N] int<lower=0, upper=1> y;

  vector[A] admin_age;                  // age per admin (months)
  vector[J] log_p;                      // unused after cleanup; bundle compat

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
  real<lower=0> sigma_lambda_prior_sd;  // unused; bundle compat
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> sigma_alpha_prior_sd;
  real<lower=0> sigma_s_prior_sd;
  real<lower=0> beta_c_prior_sd;        // unused; bundle compat
  real beta_c_prior_mean;               // unused; bundle compat
  real time_baseline;

  // Input-on-slope channel (D'). gamma_in pinned off (tight prior)
  // unless the `_dprime` variant frees it.
  real gamma_in_prior_mean;
  real<lower=0> gamma_in_prior_sd;
}

parameters {
  // Child-level (xi, zeta) -- INDEPENDENT (rho pinned 0 in D').
  matrix[2, I] z_child;                 // standardized
  vector<lower=0>[2] sigma_child;       // [1] vestigial; [2] = sigma_zeta

  real<lower=0> sigma_alpha;

  vector[J] delta_j_raw;
  real mu_delta;
  real<lower=0> tau_delta;

  real<lower=0, upper=15> s;
  real delta;

  // Input-on-slope coefficient (D').
  real gamma_in;

  // Per-child onset offset, pinned off by default (sigma_s_prior_sd ~ 0).
  vector<lower=0>[I] s_i;
  real<lower=0> sigma_s;
}

transformed parameters {
  // sigma_xi^2 = sigma_r^2 + sigma_alpha^2 (external sigma_r pin).
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  // Build child-level effects with DIAGONAL (independent) scaling:
  // rho_xi_zeta is pinned to 0 in D', so no Cholesky correlation.
  matrix[I, 2] child_effs;
  child_effs[, 1] = (sigma_xi       * z_child[1, ])';
  child_effs[, 2] = (sigma_child[2] * z_child[2, ])';

  // Sum-to-zero centering on each random-effect column (see
  // log_irt_long.stan for rationale): delta carries the full population
  // slope; xi's mean equals mu_r.
  vector[I] xi   = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta = child_effs[, 2] - mean(child_effs[, 2]);
  real<lower=0> sigma_zeta = sigma_child[2];

  // Flat (non-class-stratified) delta_j hierarchy.
  vector[J] delta_j = mu_delta + tau_delta * delta_j_raw;

  // Imputed per-child input deviation: the sigma_r share of (xi - mu_r).
  vector[I] log_r_dev;
  for (i in 1:I)
    log_r_dev[i] = (square(sigma_r) / square(sigma_xi)) * (xi[i] - mu_r);
}

model {
  // Priors
  to_vector(z_child) ~ std_normal();
  sigma_child[1] ~ normal(0, 1);                       // vestigial
  sigma_child[2] ~ normal(0, sigma_zeta_prior_sd);     // sigma_zeta
  sigma_alpha    ~ normal(0, sigma_alpha_prior_sd);

  delta_j_raw ~ std_normal();
  mu_delta  ~ normal(mu_mu_c, sigma_mu_c);
  tau_delta ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  gamma_in ~ normal(gamma_in_prior_mean, gamma_in_prior_sd);

  s_i     ~ normal(0, sigma_s);
  sigma_s ~ normal(0, sigma_s_prior_sd);

  // Likelihood: parallelize per-observation lpmf via reduce_sum.
  target += reduce_sum(partial_sum_lpmf, y, 1,
                       aa, jj, admin_to_child,
                       admin_age, s, a0,
                       time_baseline, delta, log_H, gamma_in,
                       delta_j,
                       xi, zeta, s_i, log_r_dev);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real rho_xi_zeta = 0;  // pinned to 0 in D'

  // Per-observation log-likelihood for LOO/WAIC.
  vector[N] log_lik;
  {
    vector[N] log_age;
    for (n in 1:N) {
      int ch_n = admin_to_child[aa[n]];
      log_age[n] = log(fmax(admin_age[aa[n]] - s - s_i[ch_n], 0.01) / a0);
    }
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      real slope_n = time_baseline + delta + zeta[ch] + gamma_in * log_r_dev[ch];
      real eta_n = xi[ch] + log_H + slope_n * log_age[n] - delta_j[jj[n]];
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
    }
  }

  // posterior-expected log_alpha given xi (complement of log_r_dev)
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
