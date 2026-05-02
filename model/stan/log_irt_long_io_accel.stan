// Input-uptake + acceleration model: tests the Fernald (2006) virtuous
// cycle hypothesis that input drives acceleration, not just level.
//
// Same CDI / per-recording observation structure as log_irt_io.stan.
// The key difference: instead of an independent prior on per-child
// slope deviation zeta_i, we regress zeta_i on the inferred per-child
// input rate:
//
//   zeta_i = gamma_input * log_r_true_dev_i + sigma_zeta_resid * eps_i
//
// where log_r_true_dev_i = log_r_true_i - mu_r is the per-child
// deviation from the population mean input rate. Then:
//
//   gamma_input > 0  (CrI excludes 0)  ->  input drives acceleration
//                                          (virtuous cycle confirmed)
//   gamma_input ~ 0                    ->  input drives level only
//   gamma_input < 0                    ->  compensation / anti-virtuous
//
// Designed to support a joint fit across multiple input-uptake datasets
// (BabyView, SEEDLingS, ...) with per-dataset reactivity inflation
// beta_react[d]. Each child belongs to exactly one dataset.

data {
  // ---- CDI side (same shape as log_irt_io.stan) ----
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

  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;

  // beta_c stays for class-specific frequency slopes (default unit)
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;
  real time_baseline;

  // ---- Input-uptake side (multi-dataset capable) ----
  int<lower=1> V;                                  // total recordings across datasets
  int<lower=1> D;                                  // number of datasets (e.g. 2 for BV+SD)
  array[V] int<lower=1, upper=I> video_to_child;
  array[V] int<lower=1, upper=D> dataset_per_video;
  array[I] int<lower=1, upper=D> dataset_per_child;
  vector[V] log_r_obs;
  vector[V] log_r_obs_weight;                      // per-recording weight; 0 = unused

  // Population priors
  real mu_r_prior_mean;
  real<lower=0> mu_r_prior_sd;
  real<lower=0> sigma_r_prior_sd;

  // Per-dataset measurement priors
  vector[D] beta_react_prior_mean;
  vector<lower=0>[D] beta_react_prior_sd;
  real<lower=0> sigma_within_prior_sd;

  // ---- Acceleration coupling priors ----
  real gamma_input_prior_mean;        // 0 = no virtuous cycle
  real<lower=0> gamma_input_prior_sd; // weakly informative, e.g. 0.5
  real<lower=0> sigma_zeta_resid_prior_sd;
}

parameters {
  // Population
  real mu_r;
  real<lower=0> sigma_r;
  real<lower=0> sigma_alpha;

  // Per-child latents (non-centered)
  vector[I] log_r_true_raw;
  vector[I] log_alpha_raw;
  vector[I] eps_zeta_raw;             // residual zeta (after coupling to log_r)

  // Word-level
  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;
  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;
  vector[C] beta_c;

  // Global
  real<lower=0, upper=15> s;
  real delta;

  // Measurement
  vector[D] beta_react;
  real<lower=0> sigma_within;

  // Coupling
  real gamma_input;
  real<lower=0> sigma_zeta_resid;
}

transformed parameters {
  // Per-child input rate (sum-to-zero centered on log_r_true_dev so mu_r
  // carries the population mean)
  vector[I] log_r_true_dev_uncentered = sigma_r * log_r_true_raw;
  vector[I] log_r_true_dev = log_r_true_dev_uncentered
                              - mean(log_r_true_dev_uncentered);
  vector[I] log_r_true     = mu_r + log_r_true_dev;

  // Per-child efficiency (sum-to-zero)
  vector[I] log_alpha_uncentered = sigma_alpha * log_alpha_raw;
  vector[I] log_alpha = log_alpha_uncentered - mean(log_alpha_uncentered);

  // Per-child slope deviation: regression on log_r_true_dev plus residual.
  // Sum-to-zero centering on the residual ensures mean(zeta) = 0 since
  // mean(log_r_true_dev) = 0 by construction.
  vector[I] eps_zeta_uncentered = sigma_zeta_resid * eps_zeta_raw;
  vector[I] eps_zeta = eps_zeta_uncentered - mean(eps_zeta_uncentered);
  vector[I] zeta = gamma_input * log_r_true_dev + eps_zeta;

  vector[I] xi = log_r_true + log_alpha;

  vector[J] psi;
  for (j in 1:J) psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];

  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
}

model {
  // Population priors
  mu_r          ~ normal(mu_r_prior_mean, mu_r_prior_sd);
  sigma_r       ~ normal(0, sigma_r_prior_sd);
  sigma_alpha   ~ normal(0, 1);
  log_r_true_raw ~ std_normal();
  log_alpha_raw  ~ std_normal();

  // Coupling priors
  gamma_input        ~ normal(gamma_input_prior_mean, gamma_input_prior_sd);
  sigma_zeta_resid   ~ normal(0, sigma_zeta_resid_prior_sd);
  eps_zeta_raw       ~ std_normal();

  // Word-level
  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);
  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);
  beta_c         ~ normal(beta_c_prior_mean, beta_c_prior_sd);

  // Global
  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  // Measurement (per-dataset reactivity inflation)
  for (d in 1:D) {
    beta_react[d] ~ normal(beta_react_prior_mean[d], beta_react_prior_sd[d]);
  }
  sigma_within ~ normal(0, sigma_within_prior_sd);

  // Per-recording observations
  {
    vector[V] mu_obs = log_r_true[video_to_child] + beta_react[dataset_per_video];
    log_r_obs ~ normal(mu_obs, sigma_within);
  }

  // CDI likelihood
  vector[N] eta;
  {
    vector[N] ae;
    for (n in 1:N) ae[n] = fmax(admin_age[aa[n]] - s, 0.01);
    vector[N] log_age = log(ae / a0);
    vector[N] xi_per_obs;
    vector[N] zeta_per_obs;
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      xi_per_obs[n]   = xi[ch];
      zeta_per_obs[n] = zeta[ch];
    }
    vector[N] slope_per_obs = time_baseline + delta + zeta_per_obs;
    vector[N] beta_per_obs  = beta_c[cc[jj]];
    vector[N] base = xi_per_obs + beta_per_obs .* log_p[jj] + log_H
                   + slope_per_obs .* log_age - psi[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  vector[D] reactivity_multiplier;
  for (d in 1:D) reactivity_multiplier[d] = exp(beta_react[d]);

  // Marginal sigma_zeta implied by the coupling
  real sigma_zeta_marginal = sqrt(square(gamma_input) * square(sigma_r)
                                   + square(sigma_zeta_resid));

  // Fraction of zeta variance explained by input
  real input_share_zeta = square(gamma_input) * square(sigma_r)
                          / (square(gamma_input) * square(sigma_r)
                             + square(sigma_zeta_resid));

  // Per-observation log-likelihood for LOO comparison
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
}
