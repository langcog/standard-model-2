// Strict superset of log_irt_io.stan: adds an optional comprehension
// observation channel.
//
// For comp observations on the same (admin, word) grid:
//   eta_comp = lambda * base + gamma_0 + gamma_1 * log_age
// where `base` is the same per-(admin, word) linear predictor that
// production uses. gamma_0 captures the constant logit-scale shift
// from production to comprehension; gamma_1 captures how that shift
// changes with log-age. Equivalent reparameterization on the
// difficulty side: psi_j_comp = psi_j - gamma_0 - gamma_1 * log_age.
//
// Modularity contract: when N_comp = 0 the comp block contributes
// nothing to the likelihood, and with gamma_0_prior_sd / gamma_1_prior_sd
// pinned tight near 0 the gamma parameters are pinned at their prior
// mean (0 by default). In that configuration the posterior on every
// other parameter is identical to log_irt_io.stan. The fit_io.R
// dispatcher selects this file only for variants matching `comp_*`.

data {
  // ---- CDI side (production; identical to log_irt_io.stan) ----
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
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;

  // ---- Input-observed side ----
  int<lower=1> V;
  array[V] int<lower=1, upper=I> video_to_child;
  vector[V] log_r_obs;
  vector[V] log_r_obs_weight;

  real mu_r_prior_mean;
  real<lower=0> mu_r_prior_sd;
  real beta_react_prior_mean;
  real<lower=0> beta_react_prior_sd;
  real<lower=0> sigma_r_prior_sd;
  real<lower=0> sigma_within_prior_sd;

  // ---- Comprehension side (new; pass N_comp=0 to disable) ----
  int<lower=0> N_comp;
  array[N_comp] int<lower=1, upper=A> aa_comp;
  array[N_comp] int<lower=1, upper=J> jj_comp;
  array[N_comp] int<lower=0, upper=1> y_comp;

  // Comp shift priors. Default tight at 0 (variant-overridable):
  //   gamma_0 = 0  ->  comp uses same linear predictor as production
  //   gamma_1 = 0  ->  no age-varying gap
  real gamma_0_prior_mean;
  real<lower=0> gamma_0_prior_sd;
  real gamma_1_prior_mean;
  real<lower=0> gamma_1_prior_sd;
}

parameters {
  // Population
  real mu_r;
  real<lower=0> sigma_r;
  real<lower=0> sigma_alpha;

  // Per-child latents
  vector[I] log_r_true_raw;
  vector[I] log_alpha_raw;

  // Per-child slopes
  vector[I] zeta_raw;
  real<lower=0> sigma_zeta;

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
  real beta_react;
  real<lower=0> sigma_within;

  // Comp shift
  real gamma_0;
  real gamma_1;
}

transformed parameters {
  vector[I] log_r_true_dev_uncentered = sigma_r * log_r_true_raw;
  vector[I] log_r_true_dev = log_r_true_dev_uncentered
                              - mean(log_r_true_dev_uncentered);
  vector[I] log_r_true     = mu_r + log_r_true_dev;

  vector[I] log_alpha_uncentered = sigma_alpha * log_alpha_raw;
  vector[I] log_alpha = log_alpha_uncentered - mean(log_alpha_uncentered);

  vector[I] zeta_uncentered = sigma_zeta * zeta_raw;
  vector[I] zeta = zeta_uncentered - mean(zeta_uncentered);

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

  zeta_raw   ~ std_normal();
  sigma_zeta ~ normal(0, sigma_zeta_prior_sd);

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);
  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);
  beta_c         ~ normal(beta_c_prior_mean, beta_c_prior_sd);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  // Measurement (videos)
  beta_react   ~ normal(beta_react_prior_mean, beta_react_prior_sd);
  sigma_within ~ normal(0, sigma_within_prior_sd);
  log_r_obs ~ normal(log_r_true[video_to_child] + beta_react, sigma_within);

  // Comp shift priors
  gamma_0 ~ normal(gamma_0_prior_mean, gamma_0_prior_sd);
  gamma_1 ~ normal(gamma_1_prior_mean, gamma_1_prior_sd);

  // ---- Production likelihood (same as log_irt_io.stan) ----
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
    vector[N] slope_per_obs = 1 + delta + zeta_per_obs;
    vector[N] beta_per_obs  = beta_c[cc[jj]];
    vector[N] base = xi_per_obs + beta_per_obs .* log_p[jj] + log_H
                   + slope_per_obs .* log_age - psi[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);

  // ---- Comprehension likelihood (only if N_comp > 0) ----
  if (N_comp > 0) {
    vector[N_comp] eta_c;
    {
      vector[N_comp] ae_c;
      for (n in 1:N_comp) ae_c[n] = fmax(admin_age[aa_comp[n]] - s, 0.01);
      vector[N_comp] log_age_c = log(ae_c / a0);
      vector[N_comp] xi_per_obs_c;
      vector[N_comp] zeta_per_obs_c;
      for (n in 1:N_comp) {
        int ch = admin_to_child[aa_comp[n]];
        xi_per_obs_c[n]   = xi[ch];
        zeta_per_obs_c[n] = zeta[ch];
      }
      vector[N_comp] slope_per_obs_c = 1 + delta + zeta_per_obs_c;
      vector[N_comp] beta_per_obs_c  = beta_c[cc[jj_comp]];
      vector[N_comp] base_c = xi_per_obs_c + beta_per_obs_c .* log_p[jj_comp]
                            + log_H + slope_per_obs_c .* log_age_c
                            - psi[jj_comp];
      // Comp shift sits OUTSIDE the lambda multiplier so that gamma_0,
      // gamma_1 are interpretable as a logit-scale shift independent
      // of word discrimination. With lambda == 1 (Rasch / 1PL) this
      // distinction is moot.
      eta_c = lambda[jj_comp] .* base_c + gamma_0 + gamma_1 * log_age_c;
    }
    y_comp ~ bernoulli_logit(eta_c);
  }
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real reactivity_multiplier = exp(beta_react);
  // Comp shift in production-equivalent age units: at a_0, comprehension
  // is gamma_0 logits "easier" than production; at age t the gap is
  // gamma_0 + gamma_1 * log(t/a_0).
}
