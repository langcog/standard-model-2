// Longitudinal log-linear IRT + LWL processing channel.
//
// Extends log_irt_long.stan with a second observation channel:
//   log(rt_{i,t}) ~ Normal(mu_rt + mu_rtslope * lz
//                          - gamma_rt * log_alpha_i
//                          + lz * rtslope_i, sigma_lwl)
// where lz = log(lwl_age / a0).
//
// Key change vs. log_irt_long.stan: log_alpha_i is now a free per-
// child latent (not derived by shrinkage from xi). Together with an
// independent log_r_dev_i ~ N(0, sigma_r), we have
//   xi_i = mu_r + log_r_dev_i + log_alpha_i.
// Adding the LWL channel breaks the log_r vs log_alpha exchangeability
// that holds in CDI-only fits, so the decomposition is identified by
// data, not just by the externally-pinned sigma_r.
//
// Per-child triple (log_alpha_i, zeta_i, rtslope_i) is drawn from a
// 3-D MVN with LKJ on the correlation matrix. log_r_dev_i is drawn
// independently (input rate is conceptually exogenous to processing
// efficiency and growth-rate deviation).

data {
  // ---- CDI / vocab side ---- //
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

  // Variant toggles (tight prior near 0 disables the component)
  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  // Per-class slope on log p_j. Default behavior (DEFAULT_PRIORS) pins
  // beta_c at 1 (unit-accumulator). The `no_freq*` variants pin it at 0.
  // Mirrors log_irt_long.stan / log_irt_io.stan.
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;

  // ---- LWL processing side ---- //
  int<lower=0> N_lwl;                          // LWL admins
  array[N_lwl] int<lower=1, upper=I> lwl_to_child;
  vector[N_lwl] lwl_log_age;                   // log(lwl_age / a0)
  vector[N_lwl] lwl_log_rt;                    // observed log(rt_ms)

  real mu_rt_prior_mean;
  real<lower=0> mu_rt_prior_sd;
  real<lower=0> mu_rtslope_prior_sd;
  real<lower=0> gamma_rt_prior_sd;
  real<lower=0> sigma_rtslope_prior_sd;
  real<lower=0> sigma_lwl_prior_sd;
}

parameters {
  // 3-D MVN per child: (log_alpha, zeta, rtslope)
  matrix[3, I] z_child;
  vector<lower=0>[3] sigma_child;     // (sigma_alpha, sigma_zeta, sigma_rtslope)
  cholesky_factor_corr[3] L_child;

  // log_r_dev: independent N(0, sigma_r) per child (not in the MVN)
  vector[I] log_r_dev_raw;

  // Word-level
  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  // Population structural
  real<lower=0, upper=15> s;
  real delta;

  // 2PL discrimination
  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;
  vector[C] beta_c;

  // LWL channel
  real mu_rt;
  real mu_rtslope;
  real gamma_rt;
  real<lower=0> sigma_lwl;
}

transformed parameters {
  real<lower=0> sigma_alpha   = sigma_child[1];
  real<lower=0> sigma_zeta    = sigma_child[2];
  real<lower=0> sigma_rtslope = sigma_child[3];
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  matrix[I, 3] child_effs;
  {
    matrix[3, 3] L_scaled = diag_pre_multiply(sigma_child, L_child);
    child_effs = (L_scaled * z_child)';
  }
  // Sum-to-zero centering on every random-effect column. Without this
  // the (delta, mean(zeta)) and (mu_rtslope, mean(rtslope)) splits
  // are partially unidentified -- each random-effect mean can absorb
  // part of its corresponding population fixed effect. Centering pins
  // the means at 0 so the fixed effects carry their full intended
  // role.
  vector[I] log_alpha = child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta      = child_effs[, 2] - mean(child_effs[, 2]);
  vector[I] rtslope   = child_effs[, 3] - mean(child_effs[, 3]);

  // Center log_r_dev too so mu_r captures the full population mean
  // input rate; otherwise log_r_dev's mean and log_alpha's mean
  // both contribute to E[xi] - mu_r.
  vector[I] log_r_dev_uncentered = sigma_r * log_r_dev_raw;
  vector[I] log_r_dev = log_r_dev_uncentered - mean(log_r_dev_uncentered);
  vector[I] xi = mu_r + log_r_dev + log_alpha;

  vector[J] psi;
  for (j in 1:J) psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
}

model {
  // Latent priors
  to_vector(z_child) ~ std_normal();
  log_r_dev_raw     ~ std_normal();
  sigma_child[1]    ~ normal(0, 1);                              // sigma_alpha
  sigma_child[2]    ~ normal(0, sigma_zeta_prior_sd);            // sigma_zeta
  sigma_child[3]    ~ normal(0, sigma_rtslope_prior_sd);         // sigma_rtslope
  L_child           ~ lkj_corr_cholesky(2);

  // Word-level priors
  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  // Population structural
  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);
  beta_c         ~ normal(beta_c_prior_mean, beta_c_prior_sd);

  // LWL priors
  mu_rt      ~ normal(mu_rt_prior_mean, mu_rt_prior_sd);
  mu_rtslope ~ normal(0, mu_rtslope_prior_sd);
  gamma_rt   ~ normal(0, gamma_rt_prior_sd);
  sigma_lwl  ~ normal(0, sigma_lwl_prior_sd);

  // ---- CDI likelihood (same structure as log_irt_long.stan) ---- //
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

  // ---- LWL likelihood ---- //
  if (N_lwl > 0) {
    vector[N_lwl] lwl_mean;
    for (n in 1:N_lwl) {
      int ch = lwl_to_child[n];
      lwl_mean[n] = mu_rt
                  + mu_rtslope * lwl_log_age[n]
                  - gamma_rt * log_alpha[ch]
                  + lwl_log_age[n] * rtslope[ch];
    }
    lwl_log_rt ~ normal(lwl_mean, sigma_lwl);
  }
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  // Reactivity-style summary: how many e-fold reductions in RT per
  // unit of log_alpha. (gamma_rt is in log-rt units per log_alpha.)
  real reactivity_multiplier = exp(gamma_rt);
  // Correlation matrix Sigma_corr = L_child * L_child'
  matrix[3, 3] Sigma_corr = multiply_lower_tri_self_transpose(L_child);
  real rho_alpha_zeta    = Sigma_corr[2, 1];
  real rho_alpha_rtslope = Sigma_corr[3, 1];
  real rho_zeta_rtslope  = Sigma_corr[3, 2];
}
