// Linear-in-time LMM variant of log_irt_long.stan.
//
// The headline log-linear model has eta = lambda_j * (theta_it - beta_j)
// with theta_it = xi_i + (1+delta+zeta_i) * log((age - s)/a0). With
// delta ~ 9 the time term is effectively (t/a_0)^10 -- super-linear in
// real tokens. This file replaces the time term with a plain linear
// random-intercept-and-slope LMM:
//
//   theta_it = xi_i + (beta_age + zeta_i) * (age - a0)
//
// where beta_age is the population growth rate (logits per month) and
// zeta_i is the per-child slope deviation in the same units. xi_i and
// the difficulty term beta_j are unchanged. No s, no delta, no
// (t-s)/a0 -- just linear age. Used for direct comparison against the
// log-linear baseline (LOO-CV, posterior predictive checks).

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

  // Variant toggles -- mostly the same as log_irt_long.stan, but
  // beta_age replaces (1+delta) and there is no s.
  real beta_age_prior_mean;            // pop. growth rate (logits/mo)
  real<lower=0> beta_age_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> beta_c_prior_sd;       // tight => beta_c pinned at 1
}

parameters {
  matrix[2, I] z_child;
  vector<lower=0>[2] sigma_child;
  cholesky_factor_corr[2] L_child;

  real<lower=0> sigma_alpha;

  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real beta_age;                       // pop. age slope (logits/mo)

  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;

  vector[C] beta_c;
}

transformed parameters {
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  matrix[I, 2] child_effs;
  {
    matrix[2, 2] L_scaled;
    L_scaled[1, 1] = sigma_xi * L_child[1, 1];
    L_scaled[1, 2] = 0;
    L_scaled[2, 1] = sigma_child[2] * L_child[2, 1];
    L_scaled[2, 2] = sigma_child[2] * L_child[2, 2];
    child_effs = (L_scaled * z_child)';
  }
  vector[I] xi   = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta = child_effs[, 2] - mean(child_effs[, 2]);
  real<lower=0> sigma_zeta = sigma_child[2];

  vector[J] psi;
  for (j in 1:J) psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];

  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
}

model {
  to_vector(z_child) ~ std_normal();
  sigma_child[1] ~ normal(0, 1);
  sigma_child[2] ~ normal(0, sigma_zeta_prior_sd);
  L_child       ~ lkj_corr_cholesky(2);
  sigma_alpha   ~ normal(0, 1);

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  beta_age ~ normal(beta_age_prior_mean, beta_age_prior_sd);

  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);

  beta_c ~ normal(1, beta_c_prior_sd);

  vector[N] eta;
  {
    vector[N] dt;
    for (n in 1:N) dt[n] = admin_age[aa[n]] - a0;       // linear age, NOT log
    vector[N] xi_per_obs;
    vector[N] zeta_per_obs;
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      xi_per_obs[n]   = xi[ch];
      zeta_per_obs[n] = zeta[ch];
    }
    vector[N] slope_per_obs = beta_age + zeta_per_obs;
    vector[N] beta_per_obs  = beta_c[cc[jj]];
    vector[N] base = xi_per_obs + beta_per_obs .* log_p[jj] + log_H
                   + slope_per_obs .* dt - psi[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real rho_xi_zeta = L_child[2, 1];
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
