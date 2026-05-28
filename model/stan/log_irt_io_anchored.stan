// Input-observed model with delta_j ANCHORED to an external fit.
//
// Identical to log_irt_io.stan except item difficulties delta_j are no
// longer built from a class-hierarchical prior (mu_c + tau_c * raw).
// Instead each delta_j gets an informative prior centered on the
// posterior-median delta_j from the big English longitudinal fit:
//
//     delta_j[j] ~ normal(delta_j_prior_mean[j], delta_j_prior_sd[j])
//
// Anchored items use a tight sd (~0.10); the handful of items with no
// anchor get a weak prior at the population mean (sd ~5). This lets the
// small IO samples (N = 20-66 kids) spend their data on the per-child
// input / efficiency / slope params instead of re-estimating ~600 word
// difficulties from scratch.
//
// beta_react is retained but is pinned at 0 via its prior in the IO
// bundles (passive LENA; no body-cam reactivity). Everything else --
// accumulator + 2PL + per-child slopes + input channel -- matches
// log_irt_io.stan.

data {
  // ---- CDI side ----
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

  real mu_mu_c;              // (unused here; kept for bundle compatibility)
  real<lower=0> sigma_mu_c;  // (unused here)

  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;

  // ---- delta_j anchor (NEW) ----
  vector[J] delta_j_prior_mean;
  vector<lower=0>[J] delta_j_prior_sd;

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
}

parameters {
  // Population
  real mu_r;
  real<lower=0> sigma_r;
  real<lower=0> sigma_alpha;

  // Per-child latents
  vector[I] log_r_true_raw;
  vector[I] log_alpha_raw;

  // Per-child slope deviations
  vector[I] zeta_raw;
  real<lower=0> sigma_zeta;

  // Word-level: delta_j is now free (anchored via prior), no class hierarchy
  vector[J] delta_j;
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
  vector[I] log_r_true_dev_uncentered = sigma_r * log_r_true_raw;
  vector[I] log_r_true_dev = log_r_true_dev_uncentered
                              - mean(log_r_true_dev_uncentered);
  vector[I] log_r_true     = mu_r + log_r_true_dev;

  vector[I] log_alpha_uncentered = sigma_alpha * log_alpha_raw;
  vector[I] log_alpha = log_alpha_uncentered - mean(log_alpha_uncentered);

  vector[I] zeta_uncentered = sigma_zeta * zeta_raw;
  vector[I] zeta = zeta_uncentered - mean(zeta_uncentered);

  vector[I] xi = log_r_true + log_alpha;

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

  // Slopes
  zeta_raw   ~ std_normal();
  sigma_zeta ~ normal(0, sigma_zeta_prior_sd);

  // Word-level: anchored item difficulties
  delta_j        ~ normal(delta_j_prior_mean, delta_j_prior_sd);
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
                   + slope_per_obs .* log_age - delta_j[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real reactivity_multiplier = exp(beta_react);
}
