// Longitudinal log-linear IRT accumulator with REPARAMETERIZED variance
// channels. Drop-in replacement for the slopes_si variant in
// log_irt_long.stan, designed to break the (sigma_zeta, sigma_s)
// trade-off that produced Rhat 1.5-1.9 mixing problems.
//
// Reparameterization: instead of sampling sigma_zeta and sigma_s
// directly (with diagonal-ridge geometry in their joint posterior),
// sample axis-aligned coordinates:
//   sigma_total = sqrt(sigma_zeta^2 + sigma_s^2)
//   p_zeta      = sigma_zeta^2 / sigma_total^2  ∈ (0, 1)
// with priors
//   sigma_total ~ HN(0, 5)   (broad; data dominates with N=145k)
//   p_zeta      ~ Beta(2, 2) (weakly uniform on the partition)
// Back-transform in `transformed parameters`:
//   sigma_zeta = sigma_total * sqrt(p_zeta)
//   sigma_s    = sigma_total * sqrt(1 - p_zeta)
// Both back-transformed parameters retain their original interpretation
// and are reported in summaries.
//
// Note on dimensional handling: sigma_zeta scales log-time (multiplicative
// on log_age) while sigma_s is in age-months (additive inside log()).
// Treating sigma_total^2 as their sum-of-squares is a sampling choice,
// not a claim about dimensional commensurability -- it gives Stan an
// axis-aligned coordinate system, and what the data identifies is what
// the marginal posterior reflects.
//
// Linear predictor (identical to log_irt_long.stan with independent s_i):
//   eta = lambda_j * [xi_i + beta_c[cc[j]] * log p_j + log H
//                     + (1 + delta + zeta_i) * log((age_a - s - s_i)/a0)
//                     - psi_j]
// (xi, zeta) drawn from a bivariate Normal via 2x2 Cholesky with LKJ(2)
// prior on the correlation; s_i drawn independently from Normal_+(0, sigma_s).

functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        array[] int aa, array[] int jj,
                        array[] int admin_to_child, array[] int cc,
                        vector admin_age, real s, real a0,
                        real time_baseline, real delta,
                        real log_H,
                        vector log_p,
                        vector psi, vector lambda,
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
                + slope_n * log_age_n - psi[jj[n]];
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
  // sigma_zeta_prior_sd and sigma_s_prior_sd are NOT used here -- the
  // reparameterization replaces them with sigma_total + p_zeta. Kept
  // in the data block for compatibility with the standard bundle.
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> sigma_alpha_prior_sd;
  real<lower=0> sigma_s_prior_sd;
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;
  real time_baseline;
}

parameters {
  // (xi, zeta) drawn jointly via 2x2 Cholesky with LKJ correlation.
  // sigma_zeta is DERIVED below from sigma_total + p_zeta (not sampled
  // directly).
  matrix[2, I] z_child;
  cholesky_factor_corr[2] L_child;

  real<lower=0> sigma_alpha;

  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real<lower=0, upper=15> s;
  real delta;

  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;

  vector[C] beta_c;

  // Per-child onset (centered half-normal; sigma_s derived below).
  vector<lower=0>[I] s_i;

  // REPARAMETERIZED variance channels.
  real<lower=0> sigma_total;          // total between-kid SD
  real<lower=0, upper=1> p_zeta;      // fraction of variance in zeta vs s
}

transformed parameters {
  // sigma_xi tied to sigma_alpha + sigma_r (same as base model).
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  // Back-transform: sigma_zeta and sigma_s recovered from sigma_total + p_zeta.
  // These have the same interpretation as in the base model; only the
  // sampling coordinates differ.
  real<lower=0> sigma_zeta = sigma_total * sqrt(p_zeta);
  real<lower=0> sigma_s    = sigma_total * sqrt(1 - p_zeta);

  // Build child effects for (xi, zeta) with the recovered SDs.
  matrix[I, 2] child_effs;
  {
    matrix[2, 2] L_scaled;
    L_scaled[1, 1] = sigma_xi * L_child[1, 1];
    L_scaled[1, 2] = 0;
    L_scaled[2, 1] = sigma_zeta * L_child[2, 1];
    L_scaled[2, 2] = sigma_zeta * L_child[2, 2];
    child_effs = (L_scaled * z_child)';
  }
  // Sum-to-zero center on xi and zeta (same as base model).
  vector[I] xi   = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta = child_effs[, 2] - mean(child_effs[, 2]);

  vector[J] psi;
  for (j in 1:J) {
    psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];
  }
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
}

model {
  // Priors on reparameterized variance channels.
  sigma_total ~ normal(0, 5);    // broad; data with N=145k dominates
  p_zeta      ~ beta(2, 2);       // weakly uniform on (0, 1)

  // Bivariate (xi, zeta) prior: standard non-centered with LKJ correlation.
  to_vector(z_child) ~ std_normal();
  L_child           ~ lkj_corr_cholesky(2);
  sigma_alpha       ~ normal(0, sigma_alpha_prior_sd);

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);

  beta_c ~ normal(beta_c_prior_mean, beta_c_prior_sd);

  // Per-child onset (centered half-normal with derived sigma_s).
  s_i ~ normal(0, sigma_s);

  target += reduce_sum(partial_sum_lpmf, y, 1,
                       aa, jj, admin_to_child, cc,
                       admin_age, s, a0,
                       time_baseline, delta, log_H,
                       log_p, psi, lambda, beta_c,
                       xi, zeta, s_i);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real rho_xi_zeta = L_child[2, 1];

  // Per-observation log-likelihood for PSIS-LOO.
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
                  + slope_n * log_age[n] - psi[jj[n]];
      real eta_n = lambda[jj[n]] * base_n;
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
    }
  }

  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
