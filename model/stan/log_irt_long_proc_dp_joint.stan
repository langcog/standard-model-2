// Longitudinal log-linear IRT + LWL processing REGRESSION ladder (D'0-D'3).
//
// Re-spec of log_irt_long_proc.stan. Two key changes:
//   (1) RT enters as a PREDICTOR of xi/kappa (regression), not as an
//       indicator of log_alpha. Per-child RT level/slope (rt0_i, rt1_i)
//       are latent traits MEASURED by the LWL channel; their deviations
//       predict vocab efficiency (xi) and acceleration (kappa).
//   (2) An observed-input channel (LENA): a standardized within-study log
//       input rate is a noisy read of the child's input position z_r.
//       sigma_r (input variance) and sigma_lena (LENA noise) are PINNED;
//       the LENA anchors WHO is high/low, not the variance.
//
//   xi_i    = mu_r + sigma_r*z_r_i + beta_xi*rt0_i + log_alpha_i
//   kappa_i = (1+delta) + gamma_in*sigma_r*z_r_i + beta_k0*rt0_i
//             + beta_k1*rt1_i + zeta_i
//   theta_i(t) = xi_i + kappa_i*log((t-s)/a0) + item terms
//   lwl_log_rt[n] ~ N(tau_s + rt0_i + (psi_s + rt1_i)*lz_n, sigma_lwl)
//   z_lena[v]     ~ N(z_r_i, sigma_lena)                       # pinned
//
// Residuals (log_alpha, zeta) are INDEPENDENT of input and of RT, so
// Var(xi) = sigma_r^2 (input) + beta_xi^2 sigma_rt0^2 (processing)
//           + sigma_alpha^2 (residual) is a clean orthogonal partition.
//
// Ladder via prior SDs (tiny -> pinned at 0):
//   D'0 = {gamma_in on, betas off};  D'1 +beta_xi;  D'2 +beta_k0;  D'3 +beta_k1.

data {
  // ---- CDI / vocab side ---- //
  int<lower=1> N;
  int<lower=1> A;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> C;
  int<lower=1> S;                              // studies (datasets)

  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;
  array[N] int<lower=0, upper=1> y;
  array[I] int<lower=1, upper=S> study_of_child;

  vector[A] admin_age;
  vector[J] log_p;

  real log_H;
  real<lower=0> a0;

  real mu_r;
  real<lower=0> sigma_r;                        // PINNED input SD
  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> sigma_alpha_prior_sd;
  real beta_c_prior_mean;
  real<lower=0> beta_c_prior_sd;

  // ---- LWL processing side ---- //
  int<lower=0> N_lwl;
  array[N_lwl] int<lower=1, upper=I> lwl_to_child;
  vector[N_lwl] lwl_log_age;                    // log(lwl_age / a0)
  vector[N_lwl] lwl_log_rt;                     // observed log(rt_ms)

  real mu_rt_prior_mean;
  real<lower=0> mu_rt_prior_sd;
  real<lower=0> mu_rtslope_prior_sd;
  real<lower=0> sigma_rt0_prior_sd;
  real<lower=0> sigma_rt1_prior_sd;
  real<lower=0> sigma_lwl_prior_sd;

  // ---- Observed-input (LENA) channel ---- //
  int<lower=0> V;                               // input observations
  array[V] int<lower=1, upper=I> rec_to_child;
  vector[V] z_lena;                             // within-study standardized log input
  vector<lower=0>[S] sigma_lena;                // PINNED per-study input measurement noise (sd units)

  // ---- Ladder regression coefficients (prior sds toggle rungs) ---- //
  real<lower=0> gamma_in_prior_sd;              // D'0 input-on-slope
  real<lower=0> beta_xi_prior_sd;               // D'1 rt0 -> xi
  real<lower=0> beta_k0_prior_sd;               // D'2 rt0 -> kappa
  real<lower=0> beta_k1_prior_sd;               // D'3 rt1 -> kappa
}

parameters {
  // Per-child RT latents (bivariate level/slope)
  matrix[2, I] z_rt;
  vector<lower=0>[2] sigma_rt;                  // (sigma_rt0, sigma_rt1)
  cholesky_factor_corr[2] L_rt;

  // Per-child residuals (independent of input & RT)
  vector[I] log_alpha_raw;
  vector[I] zeta_raw;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_zeta;

  // Input position (standardized); log_r_dev = sigma_r * z_r
  vector[I] z_r;

  // Per-study RT level / slope means
  vector[S] tau_s;
  vector[S] psi_s;

  // Word-level
  vector[J] delta_j_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;
  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;
  vector[C] beta_c;

  // Population structural
  real<lower=0, upper=15> s;
  real delta;
  real<lower=0> sigma_lwl;

  // Regression coefficients (ladder)
  real gamma_in;
  real beta_xi;
  real beta_k0;
  real beta_k1;
}

transformed parameters {
  real<lower=0> sigma_rt0 = sigma_rt[1];
  real<lower=0> sigma_rt1 = sigma_rt[2];

  // RT latents (sum-to-zero centered within sample; per-study level/slope
  // live in tau_s/psi_s, so these are pure deviations)
  matrix[I, 2] rt_effs;
  {
    matrix[2, 2] L_scaled = diag_pre_multiply(sigma_rt, L_rt);
    rt_effs = (L_scaled * z_rt)';
  }
  vector[I] rt0 = rt_effs[, 1] - mean(rt_effs[, 1]);
  vector[I] rt1 = rt_effs[, 2] - mean(rt_effs[, 2]);

  vector[I] log_alpha = sigma_alpha * log_alpha_raw;
  log_alpha = log_alpha - mean(log_alpha);
  vector[I] zeta = sigma_zeta * zeta_raw;
  zeta = zeta - mean(zeta);

  // Input deviation (z_r centered so mu_r carries the population mean)
  vector[I] z_r_c = z_r - mean(z_r);
  vector[I] log_r_dev = sigma_r * z_r_c;

  // Ability intercept and slope per child
  vector[I] xi    = mu_r + log_r_dev + beta_xi * rt0 + log_alpha;
  vector[I] kappa = 1 + delta + gamma_in * log_r_dev
                    + beta_k0 * rt0 + beta_k1 * rt1 + zeta;

  vector[J] delta_j;
  for (j in 1:J) delta_j[j] = mu_c[cc[j]] + tau_c[cc[j]] * delta_j_raw[j];
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
}

model {
  // Latent priors
  to_vector(z_rt) ~ std_normal();
  sigma_rt[1] ~ normal(0, sigma_rt0_prior_sd);
  sigma_rt[2] ~ normal(0, sigma_rt1_prior_sd);
  L_rt ~ lkj_corr_cholesky(2);

  log_alpha_raw ~ std_normal();
  zeta_raw      ~ std_normal();
  sigma_alpha   ~ normal(0, sigma_alpha_prior_sd);
  sigma_zeta    ~ normal(0, sigma_zeta_prior_sd);

  z_r ~ std_normal();

  tau_s ~ normal(mu_rt_prior_mean, mu_rt_prior_sd);
  psi_s ~ normal(0, mu_rtslope_prior_sd);

  // Word-level priors
  delta_j_raw ~ std_normal();
  mu_c  ~ normal(mu_mu_c, sigma_mu_c);
  tau_c ~ normal(0, 1);
  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);
  beta_c         ~ normal(beta_c_prior_mean, beta_c_prior_sd);

  // Population structural
  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);
  sigma_lwl ~ normal(0, sigma_lwl_prior_sd);

  // Regression coefficients (tiny prior sd -> pinned at 0)
  gamma_in ~ normal(0, gamma_in_prior_sd);
  beta_xi  ~ normal(0, beta_xi_prior_sd);
  beta_k0  ~ normal(0, beta_k0_prior_sd);
  beta_k1  ~ normal(0, beta_k1_prior_sd);

  // ---- CDI likelihood ---- //
  vector[N] eta;
  {
    vector[N] ae;
    for (n in 1:N) ae[n] = fmax(admin_age[aa[n]] - s, 0.01);
    vector[N] log_age = log(ae / a0);
    vector[N] xi_per_obs;
    vector[N] kappa_per_obs;
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      xi_per_obs[n]    = xi[ch];
      kappa_per_obs[n] = kappa[ch];
    }
    vector[N] beta_per_obs = beta_c[cc[jj]];
    vector[N] base = xi_per_obs + beta_per_obs .* log_p[jj] + log_H
                   + kappa_per_obs .* log_age - delta_j[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);

  // ---- LWL measurement likelihood ---- //
  if (N_lwl > 0) {
    vector[N_lwl] lwl_mean;
    for (n in 1:N_lwl) {
      int ch = lwl_to_child[n];
      int st = study_of_child[ch];
      lwl_mean[n] = tau_s[st] + rt0[ch]
                  + (psi_s[st] + rt1[ch]) * lwl_log_age[n];
    }
    lwl_log_rt ~ normal(lwl_mean, sigma_lwl);
  }

  // ---- Observed-input likelihood ---- //
  if (V > 0) {
    vector[V] z_mean;
    vector[V] sd_v;
    for (v in 1:V) {
      z_mean[v] = z_r_c[rec_to_child[v]];
      sd_v[v]   = sigma_lena[study_of_child[rec_to_child[v]]];  // per-study input noise
    }
    z_lena ~ normal(z_mean, sd_v);
  }
}

generated quantities {
  // Intercept-side variance partition (the payoff)
  real var_input_xi = square(sigma_r);
  real var_proc_xi  = square(beta_xi) * square(sigma_rt0);
  real var_resid_xi = square(sigma_alpha);
  real sigma_xi = sqrt(var_input_xi + var_proc_xi + var_resid_xi);
  real share_input_xi = var_input_xi / square(sigma_xi);
  real share_proc_xi  = var_proc_xi  / square(sigma_xi);
  real share_resid_xi = var_resid_xi / square(sigma_xi);

  // Slope-side variance partition
  real var_input_k = square(gamma_in) * square(sigma_r);
  real var_proc_k  = square(beta_k0) * square(sigma_rt0)
                   + square(beta_k1) * square(sigma_rt1);
  real var_resid_k = square(sigma_zeta);
  real sigma_kappa = sqrt(var_input_k + var_proc_k + var_resid_k);

  // RT level/slope correlation
  matrix[2, 2] Rt_corr = multiply_lower_tri_self_transpose(L_rt);
  real rho_rt = Rt_corr[2, 1];

  // Pointwise CDI log-likelihood for LOO/WAIC across the ladder (the rungs
  // differ in how well xi/kappa -- hence vocab -- are explained).
  vector[N] log_lik;
  {
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      real ae = fmax(admin_age[aa[n]] - s, 0.01);
      real base = xi[ch] + beta_c[cc[jj[n]]] * log_p[jj[n]] + log_H
                + kappa[ch] * log(ae / a0) - delta_j[jj[n]];
      log_lik[n] = bernoulli_logit_lpmf(y[n] | lambda[jj[n]] * base);
    }
  }
}
