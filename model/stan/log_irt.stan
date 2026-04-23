// Log-linear IRT accumulator model for early word learning.
//
// - Child-level latent collapsed into xi_i = log r_i + log alpha_i
//   with xi_i ~ N(mu_r, sigma_r^2 + sigma_alpha^2); sigma_r pinned externally.
// - Word threshold psi_j hierarchical by lexical class.
// - Global start time s, age rate-change exponent delta.
// - 2PL discrimination lambda_j per word (toggled via data flag):
//     eta_ij = lambda_j * [xi_i + log p_j + log H + (1+delta) log((a-s)/a0) - psi_j]
//   In Rasch mode (use_2pl=0), sigma_lambda_prior_sd is passed as ~0 so
//   lambda_j is pinned at 1 (log_lambda = 0). In 2PL mode, sigma_lambda
//   is freely estimated.

data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> C;

  array[N] int<lower=1, upper=I> ii;
  array[N] int<lower=1, upper=J> jj;
  array[J] int<lower=1, upper=C> cc;
  array[N] int<lower=0, upper=1> y;

  vector[I] age;
  vector[J] log_p;

  real log_H;
  real<lower=0> a0;

  real mu_r;
  real<lower=0> sigma_r;

  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  // Diagnostic hyperpriors: near-zero SD effectively fixes parameter.
  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;

  // 2PL toggle: sigma_lambda_prior_sd near 0 -> Rasch (lambda ≈ 1).
  //   use 1.0 for 2PL, 0.001 for Rasch.
  real<lower=0> sigma_lambda_prior_sd;

  // Per-child slope toggle: sigma_zeta_prior_sd near 0 disables per-child
  // random slopes (zeta_i ≈ 0). Use 1.0 to estimate freely, 0.001 to disable.
  // zeta_i adds to delta: effective slope = (1 + delta + zeta_i).
  real<lower=0> sigma_zeta_prior_sd;
}

parameters {
  vector[I] xi_raw;
  real<lower=0> sigma_alpha;

  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real<lower=0, upper=15> s;
  real delta;

  // 2PL word discrimination, non-centered
  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;

  // Per-child slope deviations, non-centered
  vector[I] zeta_raw;
  real<lower=0> sigma_zeta;
}

transformed parameters {
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));
  vector[I] xi = mu_r + sigma_xi * xi_raw;
  vector[J] psi;
  for (j in 1:J) {
    psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];
  }
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
  vector[I] zeta = sigma_zeta * zeta_raw;
}

model {
  xi_raw      ~ std_normal();
  sigma_alpha ~ normal(0, 1);

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);

  zeta_raw   ~ std_normal();
  sigma_zeta ~ normal(0, sigma_zeta_prior_sd);

  vector[N] eta;
  {
    vector[N] ae;
    for (n in 1:N) ae[n] = fmax(age[ii[n]] - s, 0.01);
    vector[N] log_age = log(ae / a0);
    // effective slope per child = (1 + delta + zeta_i)
    vector[N] slope_per_obs = 1 + delta + zeta[ii];
    vector[N] base = xi[ii] + log_p[jj] + log_H
                   + slope_per_obs .* log_age - psi[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  // Posterior-expected log_alpha given xi (closed form).
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
