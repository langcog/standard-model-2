// Longitudinal log-linear IRT + LWL processing REGRESSION ladder (D'0-D'3).
//
// Re-spec of log_irt_long_proc.stan. Two key changes:
//   (1) RT enters as a PREDICTOR of xi/kappa (regression), not as an
//       indicator of log_alpha. Per-child RT level/slope (rt0_i, rt1_i)
//       are latent traits MEASURED by the LWL channel; their deviations
//       predict vocab efficiency (xi) and acceleration (kappa).
//   (2) MEASUREMENT-MODEL re-spec of the input channel (this _mm variant).
//       RAW per-recording log input is a noisy read of the child's latent
//       within-study input deviation d_i. sigma_r (between-child input SD) is
//       now ESTIMATED with a literature prior (symmetric with RT's sigma_rt0);
//       sigma_meas is per-INSTRUMENT noise (head-cam vs LENA), mu_r_s per-study
//       level. d_i is centered WITHIN study (between-study/instrument level
//       lives in mu_r_s, not in d), so xi's input component is within-study only.
//
//   xi_i    = mu_r + d_i + beta_xi*rt0_i + log_alpha_i,   d_i = sigma_r*z_r_i
//   kappa_i = (1+delta) + gamma_in*d_i + beta_k0*rt0_i + beta_k1*rt1_i + zeta_i
//   theta_i(t) = xi_i + kappa_i*log((t-s)/a0) + item terms
//   lwl_log_rt[n]    ~ N(tau_s + rt0_i + (psi_s + rt1_i)*lz_n, sigma_lwl)
//   log_input_obs[v] ~ N(mu_r_s[study_i] + d_i, sigma_meas[instr_v])  # estimated
//
// Residuals (log_alpha, zeta) are INDEPENDENT of input and of RT, so
// Var(xi) = sigma_r^2 (input) + beta_xi^2 sigma_rt0^2 (processing)
//           + sigma_alpha^2 (residual) is a clean orthogonal partition.
//
// Ladder via prior SDs (tiny -> pinned at 0):
//   D'0 = {gamma_in on, betas off};  D'1 +beta_xi;  D'2 +beta_k0;  D'3 +beta_k1.

functions {
  // Parallelized CDI log-likelihood (reduce_sum). Per-obs work is intentionally
  // minimal -- the heavy per-admin (admin_base) and per-item (item_offset) terms
  // are precomputed in transformed parameters -- so this scales near-linearly
  // with threads_per_chain.
  real partial_sum_lpmf(array[] int y_slice, int start, int end,
                        array[] int aa, array[] int jj,
                        vector admin_base, vector item_offset, vector lambda) {
    int ns = end - start + 1;
    vector[ns] eta;
    for (i in 1:ns) {
      int o = start + i - 1;
      eta[i] = lambda[jj[o]] * (admin_base[aa[o]] + item_offset[jj[o]]);
    }
    return bernoulli_logit_lpmf(y_slice | eta);
  }
}

data {
  // ---- CDI / vocab side ---- //
  int<lower=1> N;
  int<lower=1> grainsize;                       // reduce_sum chunk hint (TBB tunes)
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

  real mu_r;                                    // xi input-intercept anchor (data)
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

  // RT priors carry informative MEANS (frank_etal_2026, centered at a0=21mo)
  real mu_rt_prior_mean;                        // pop log-RT level @ a0  (~6.84)
  real<lower=0> mu_rt_prior_sd;
  real psi_prior_mean;                          // pop RT age-slope        (~-0.35)
  real<lower=0> mu_rtslope_prior_sd;
  real sigma_rt0_prior_mean;                    // between-child RT level SD (~0.143)
  real<lower=0> sigma_rt0_prior_sd;
  real sigma_rt1_prior_mean;                    // between-child RT slope SD (~0.26)
  real<lower=0> sigma_rt1_prior_sd;
  real sigma_lwl_prior_mean;                    // RT residual            (~0.24)
  real<lower=0> sigma_lwl_prior_sd;

  // ---- Observed-input channel: MEASUREMENT MODEL on RAW per-recording input ---- //
  int<lower=0> V_obs;                           // per-recording input observations
  array[V_obs] int<lower=1, upper=I> rec_to_child;
  vector[V_obs] log_input_obs;                  // RAW log input rate (NOT standardized)
  int<lower=1> n_instr;                         // measurement instruments (head-cam, LENA)
  array[V_obs] int<lower=1, upper=n_instr> instr_of_rec;
  real sigma_r_prior_mean;                      // input-rate literature anchor (~0.44)
  real<lower=0> sigma_r_prior_sd;
  vector[S] mu_r_s_prior_mean;                  // per-study observed-input level anchor
  real<lower=0> mu_r_s_prior_sd;
  real<lower=0> sigma_meas_prior_sd;            // per-instrument input-noise prior scale

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

  // Input measurement model: per-child position z_r (d_i = sigma_r*z_r_within),
  // between-child input SD sigma_r (ESTIMATED), per-study level, per-instr noise
  vector[I] z_r;
  real<lower=0> sigma_r;
  vector[S] mu_r_s;
  vector<lower=0>[n_instr] sigma_meas;

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

  // Input deviation, centered WITHIN study: between-study/instrument level lives
  // in mu_r_s (input likelihood), so d_i carries only within-study variation.
  vector[I] z_r_c;
  {
    vector[S] ssum = rep_vector(0, S);
    vector[S] scnt = rep_vector(0, S);
    for (i in 1:I) { ssum[study_of_child[i]] += z_r[i]; scnt[study_of_child[i]] += 1; }
    for (i in 1:I) z_r_c[i] = z_r[i] - ssum[study_of_child[i]] / scnt[study_of_child[i]];
  }
  vector[I] log_r_dev = sigma_r * z_r_c;        // d_i

  // Ability intercept and slope per child
  vector[I] xi    = mu_r + log_r_dev + beta_xi * rt0 + log_alpha;
  vector[I] kappa = 1 + delta + gamma_in * log_r_dev
                    + beta_k0 * rt0 + beta_k1 * rt1 + zeta;

  vector[J] delta_j;
  for (j in 1:J) delta_j[j] = mu_c[cc[j]] + tau_c[cc[j]] * delta_j_raw[j];
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);

  // Precompute per-admin and per-item terms so the N-obs CDI likelihood does
  // minimal per-obs work: eta[n] = lambda[j] * (admin_base[a] + item_offset[j]).
  // This collapses the autodiff graph from O(N) heavy nodes to O(A+J), and is
  // what reduce_sum then parallelizes. (A=admins, J=items, both << N.)
  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + kappa[ch] * log(fmax(admin_age[a] - s, 0.01) / a0);
  }
  vector[J] item_offset;
  for (j in 1:J) item_offset[j] = beta_c[cc[j]] * log_p[j] - delta_j[j];
}

model {
  // Latent priors
  to_vector(z_rt) ~ std_normal();
  sigma_rt[1] ~ normal(sigma_rt0_prior_mean, sigma_rt0_prior_sd);  // <lower=0> => truncated
  sigma_rt[2] ~ normal(sigma_rt1_prior_mean, sigma_rt1_prior_sd);
  L_rt ~ lkj_corr_cholesky(2);

  log_alpha_raw ~ std_normal();
  zeta_raw      ~ std_normal();
  sigma_alpha   ~ normal(0, sigma_alpha_prior_sd);
  sigma_zeta    ~ normal(0, sigma_zeta_prior_sd);

  // Input measurement-model priors
  z_r        ~ std_normal();
  sigma_r    ~ normal(sigma_r_prior_mean, sigma_r_prior_sd);   // <lower=0> => truncated
  mu_r_s     ~ normal(mu_r_s_prior_mean, mu_r_s_prior_sd);
  sigma_meas ~ normal(0, sigma_meas_prior_sd);

  // RT priors (frank_etal_2026 informative means)
  tau_s       ~ normal(mu_rt_prior_mean, mu_rt_prior_sd);
  psi_s       ~ normal(psi_prior_mean, mu_rtslope_prior_sd);

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
  sigma_lwl ~ normal(sigma_lwl_prior_mean, sigma_lwl_prior_sd);  // <lower=0> => truncated

  // Regression coefficients (tiny prior sd -> pinned at 0)
  gamma_in ~ normal(0, gamma_in_prior_sd);
  beta_xi  ~ normal(0, beta_xi_prior_sd);
  beta_k0  ~ normal(0, beta_k0_prior_sd);
  beta_k1  ~ normal(0, beta_k1_prior_sd);

  // ---- CDI likelihood (parallelized across cores via reduce_sum) ---- //
  target += reduce_sum(partial_sum_lpmf, y, grainsize,
                       aa, jj, admin_base, item_offset, lambda);

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

  // ---- Observed-input measurement likelihood (RAW per-recording log input) ---- //
  if (V_obs > 0) {
    vector[V_obs] in_mean;
    vector[V_obs] in_sd;
    for (v in 1:V_obs) {
      int ch = rec_to_child[v];
      in_mean[v] = mu_r_s[study_of_child[ch]] + log_r_dev[ch];  // study level + d_i
      in_sd[v]   = sigma_meas[instr_of_rec[v]];                 // per-instrument noise
    }
    log_input_obs ~ normal(in_mean, in_sd);
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

  // (Per-obs log_lik for LOO removed -- LOO is skipped for these big fits, and
  //  the N-vector cost it added to every saved draw was pure waste.)
}
