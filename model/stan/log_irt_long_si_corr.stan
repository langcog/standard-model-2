// Longitudinal log-linear IRT accumulator model WITH correlated per-child
// (xi, zeta, s_i) random effects.
//
// Extends log_irt_long.stan to model correlation between per-child onset
// s_i and (xi, zeta). Replaces the bivariate (xi, zeta) Cholesky in the
// base model with a trivariate Cholesky for (xi, zeta, s_i_lat), with
// LKJ(2) prior on the 3x3 correlation matrix. Positivity of s_i is
// enforced via Tobit-style clipping: s_i = fmax(s_i_lat, 0) when used in
// the linear predictor. This allows s_i to share a linear correlation
// structure with zeta (and xi) without breaking the LKJ machinery.
//
// Motivation: hypothesis is that kids with higher kappa (faster log-time
// scaling) may correlate with onset s_i -- e.g., late-onset kids may
// later accelerate faster ("catch-up" learners), or early-onset kids may
// also be faster. The bivariate (xi, zeta, indep s_i) model can't
// represent this; this file's trivariate LKJ can.
//
// Tobit clipping consequence: marginally, P(s_i = 0) ~ 0.5 (latent half
// is negative). The "delayed" subgroup is the half with s_i_lat > 0;
// sigma_s parametrizes the SD of s_i_lat. Effective mean delay across
// all kids ~ sigma_s / sqrt(2*pi) ~ 0.4 * sigma_s.
//
// Linear predictor (identical to log_irt_long.stan):
//   eta = lambda_j * [xi_i + beta_c[cc[j]] * log p_j + log H
//                     + (1 + delta + zeta_i) * log((age_a - s - s_i_eff)/a0)
//                     - delta_j]
// where s_i_eff = fmax(s_i_lat[i], 0).

functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        array[] int aa, array[] int jj,
                        array[] int admin_to_child, array[] int cc,
                        vector admin_age, real s, real a0,
                        real time_baseline, real delta,
                        real log_H,
                        vector log_p,
                        vector delta_j, vector lambda,
                        vector beta_c,
                        vector xi, vector zeta, vector s_i) {
    int n_slice = end - start + 1;
    vector[n_slice] eta_slice;
    for (i in 1:n_slice) {
      int n = start + i - 1;
      int ch = admin_to_child[aa[n]];
      real ae = fmax(admin_age[aa[n]] - s - s_i[ch], 0.01);
      real log_age_n = log(ae / a0);
      real slope_n = time_baseline + delta + zeta[ch];
      real beta_n  = beta_c[cc[jj[n]]];
      real base = xi[ch] + beta_n * log_p[jj[n]] + log_H
                + slope_n * log_age_n - delta_j[jj[n]];
      eta_slice[i] = lambda[jj[n]] * base;
    }
    return bernoulli_logit_lpmf(y_slice | eta_slice);
  }
}

data {
  int<lower=1> N;
  int<lower=1> A;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> C;

  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;
  array[N] int<lower=0, upper=1> y;

  vector[A] admin_age;
  vector[J] log_p;

  real log_H;
  real<lower=0> a0;

  real mu_r;
  real<lower=0> sigma_r;

  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> sigma_alpha_prior_sd;
  real<lower=0> sigma_s_prior_sd;
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;
  real time_baseline;
}

parameters {
  // Child-level trivariate (xi_raw, zeta_raw, s_lat_raw) with LKJ on
  // the 3x3 correlation matrix. Non-centered: z_child are standardized.
  matrix[3, I] z_child;
  vector<lower=0>[3] sigma_child;   // (sigma_xi_placeholder, sigma_zeta, sigma_s)
  cholesky_factor_corr[3] L_child;

  real<lower=0> sigma_alpha;

  vector[J] delta_j_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real<lower=0, upper=15> s;
  real delta;

  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;

  vector[C] beta_c;
}

transformed parameters {
  // Enforce sigma_child[1]^2 = sigma_r^2 + sigma_alpha^2 (xi SD constraint).
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  matrix[I, 3] child_effs;
  {
    matrix[3, 3] L_scaled;
    // Row 1 (xi): use sigma_xi (not sigma_child[1])
    L_scaled[1, 1] = sigma_xi * L_child[1, 1];
    L_scaled[1, 2] = 0;
    L_scaled[1, 3] = 0;
    // Row 2 (zeta): use sigma_child[2]
    L_scaled[2, 1] = sigma_child[2] * L_child[2, 1];
    L_scaled[2, 2] = sigma_child[2] * L_child[2, 2];
    L_scaled[2, 3] = 0;
    // Row 3 (s_lat): use sigma_child[3] (= sigma_s)
    L_scaled[3, 1] = sigma_child[3] * L_child[3, 1];
    L_scaled[3, 2] = sigma_child[3] * L_child[3, 2];
    L_scaled[3, 3] = sigma_child[3] * L_child[3, 3];
    child_effs = (L_scaled * z_child)';
  }
  // Sum-to-zero centering on xi and zeta (same as bivariate model).
  // We do NOT center s_lat -- its population mean is structurally 0 by
  // the MVN draw, and centering would shift fmax(s_lat, 0) marginal in
  // an interpretation-breaking way.
  vector[I] xi    = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta  = child_effs[, 2] - mean(child_effs[, 2]);
  vector[I] s_lat = child_effs[, 3];
  // Tobit-style positivity for the time-shift: kids with s_lat <= 0 get
  // s_i = 0 (effective onset = population mean s). Kids with s_lat > 0
  // get delayed onset s + s_lat.
  vector[I] s_i;
  for (i in 1:I) s_i[i] = fmax(s_lat[i], 0);

  real<lower=0> sigma_zeta = sigma_child[2];
  real<lower=0> sigma_s    = sigma_child[3];

  vector[J] delta_j;
  for (j in 1:J) {
    delta_j[j] = mu_c[cc[j]] + tau_c[cc[j]] * delta_j_raw[j];
  }
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
}

model {
  // Priors
  to_vector(z_child) ~ std_normal();
  sigma_child[1] ~ normal(0, 1);                    // harmless (replaced)
  sigma_child[2] ~ normal(0, sigma_zeta_prior_sd);  // sigma_zeta
  sigma_child[3] ~ normal(0, sigma_s_prior_sd);     // sigma_s
  L_child        ~ lkj_corr_cholesky(2);            // mild toward 0 corrs
  sigma_alpha    ~ normal(0, sigma_alpha_prior_sd);

  delta_j_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);

  beta_c ~ normal(beta_c_prior_mean, beta_c_prior_sd);

  target += reduce_sum(partial_sum_lpmf, y, 1,
                       aa, jj, admin_to_child, cc,
                       admin_age, s, a0,
                       time_baseline, delta, log_H,
                       log_p, delta_j, lambda, beta_c,
                       xi, zeta, s_i);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));

  // Extract the three pairwise correlations from L_child.
  // L_child is the Cholesky factor of the correlation matrix, so
  // Corr = L_child * L_child' is the correlation matrix.
  matrix[3, 3] corr_child = L_child * L_child';
  real rho_xi_zeta = corr_child[1, 2];
  real rho_xi_s    = corr_child[1, 3];
  real rho_zeta_s  = corr_child[2, 3];

  // Per-observation log-likelihood
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
      real base_n = xi[ch] + beta_c[cc[jj[n]]] * log_p[jj[n]] + log_H
                  + slope_n * log_age[n] - delta_j[jj[n]];
      real eta_n = lambda[jj[n]] * base_n;
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
    }
  }

  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
