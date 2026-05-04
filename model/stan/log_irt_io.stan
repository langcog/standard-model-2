// Input-observed extension of log_irt_long.stan.
//
// Same accumulator + 2PL + per-child slopes structure as the longitudinal
// model, but instead of pinning the population input distribution
// from external data, we observe per-video log token rates and infer
// per-child true input rate, the population mean / SD of input rate, a
// reactivity bias, and within-child measurement noise jointly.
//
// Per-video:    log_r_obs[v] ~ N(log_r_true[child[v]] + beta_react, sigma_within)
// Per-child:    log_r_true[i] ~ N(mu_r, sigma_r)
//               log_alpha[i] ~ N(0, sigma_alpha)
//               zeta[i] ~ N(0, sigma_zeta)
//               xi[i] = log_r_true[i] + log_alpha[i]
//
// pi_alpha = sigma_alpha^2 / (sigma_alpha^2 + sigma_r^2) — same definition,
// but now both numerator and denominator are estimated from data.

data {
  // ---- CDI side (same shape as log_irt_long.stan) ----
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
  // beta_c is a per-class slope on log p_j. Default behaviour is the
  // unit accumulator (beta_c == 1) via DEFAULT_PRIORS{mean=1, sd=0.001}.
  // The `no_freq` variant pins beta_c at 0 (drops frequency entirely);
  // `class_beta` frees it (sd=0.5). Mirrors log_irt_long.stan.
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;

  // ---- Input-observed side (new) ----
  int<lower=1> V;                        // total videos
  array[V] int<lower=1, upper=I> video_to_child;
  vector[V] log_r_obs;                   // log tokens/hr per video
  vector[V] log_r_obs_weight;            // per-video weight (e.g. log duration); use 0s if unused

  // Priors on the input-side scalars
  real mu_r_prior_mean;       // weakly anchored at Sperry mean ≈ 7.34
  real<lower=0> mu_r_prior_sd;
  real beta_react_prior_mean; // ~0.4 = videos ~1.4× real life
  real<lower=0> beta_react_prior_sd;
  real<lower=0> sigma_r_prior_sd;
  real<lower=0> sigma_within_prior_sd;
}

parameters {
  // Population
  real mu_r;
  real<lower=0> sigma_r;
  real<lower=0> sigma_alpha;

  // Per-child latents
  vector[I] log_r_true_raw;
  vector[I] log_alpha_raw;

  // Per-child slope deviations (toggled by sigma_zeta_prior_sd)
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
}

transformed parameters {
  // Sum-to-zero centering on each per-child latent. Without this the
  // (mu_r, mean(log_r_true)) split (and (delta, mean(zeta)) split)
  // are partially unidentified -- the random-effect means can absorb
  // part of the corresponding population fixed effect. Centering
  // pins those means and lets mu_r / delta carry their intended role.
  vector[I] log_r_true_dev_uncentered = sigma_r * log_r_true_raw;
  vector[I] log_r_true_dev = log_r_true_dev_uncentered
                              - mean(log_r_true_dev_uncentered);
  vector[I] log_r_true     = mu_r + log_r_true_dev;

  vector[I] log_alpha_uncentered = sigma_alpha * log_alpha_raw;
  vector[I] log_alpha = log_alpha_uncentered - mean(log_alpha_uncentered);

  vector[I] zeta_uncentered = sigma_zeta * zeta_raw;
  vector[I] zeta = zeta_uncentered - mean(zeta_uncentered);

  vector[I] xi         = log_r_true + log_alpha;

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

  // Slopes (toggled)
  zeta_raw   ~ std_normal();
  sigma_zeta ~ normal(0, sigma_zeta_prior_sd);

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

  // Measurement (videos)
  beta_react   ~ normal(beta_react_prior_mean, beta_react_prior_sd);
  sigma_within ~ normal(0, sigma_within_prior_sd);
  log_r_obs ~ normal(log_r_true[video_to_child] + beta_react, sigma_within);

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
    vector[N] slope_per_obs = 1 + delta + zeta_per_obs;
    vector[N] beta_per_obs  = beta_c[cc[jj]];
    vector[N] base = xi_per_obs + beta_per_obs .* log_p[jj] + log_H
                   + slope_per_obs .* log_age - psi[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);
}

generated quantities {
  // Variance decomposition: now BOTH sigma_alpha and sigma_r are
  // estimated from data, not assumed.
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  // Reactivity multiplier: videos / real-life ratio
  real reactivity_multiplier = exp(beta_react);
}
