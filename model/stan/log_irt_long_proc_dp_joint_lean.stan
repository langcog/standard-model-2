// Longitudinal log-linear IRT + LWL processing REGRESSION ladder -- "io-proc LEAN".
//
// Lean re-spec of log_irt_long_proc_dp_joint_mm.stan. Two simplifications:
//   (1) RT measurement model is LEVEL-ONLY: each child has a single RT-level deviation
//       rt0_i (a detrended residual), with per-study intercepts tau_s and ONE global
//       age slope psi. This matches the glmer detrend log(rt) ~ log(age) + (1|dataset),
//       done jointly (so rt0's uncertainty propagates). Dropped: per-child RT slope
//       rt1_i (1-3 sessions/kid can't identify it) and per-study psi_s.
//   (2) RAW coefficients (gamma_in, beta_xi, beta_k0), NOT divided by an estimated SD.
//       Common-scale priors are set in the driver as prior SD = tau / sigma_ref (a FIXED
//       reference scale), giving symmetric per-SD-effect priors across channels WITHOUT
//       the divide-by-estimated-SD funnel that broke the per-SD reparam (the _mm variant).
//
//   xi_i    = mu_r + d_i + beta_xi*rt0_i + log_alpha_i,    d_i = sigma_r*z_r_c_i
//   kappa_i = (1+delta) + gamma_in*d_i + beta_k0*rt0_i + zeta_i
//   theta_i(t) = xi_i + kappa_i*log((t-s)/a0) + item terms
//   lwl_log_rt[n]    ~ N(tau_s + rt0_i + psi*lz_n, sigma_lwl)        # single global psi
//   log_input_obs[v] ~ N(mu_r_s[study_i] + d_i, sigma_meas[instr_v])
//
// Ladder via prior SDs (tiny -> pinned at 0):
//   D'0 = {gamma_in on, betas off};  D'1 +beta_xi;  D'2 +beta_k0.  (no D'3 -- no rt1)

functions {
  // Parallelized CDI log-likelihood (reduce_sum); heavy per-admin/per-item terms are
  // precomputed in transformed parameters so this scales near-linearly with threads.
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
  int<lower=1> grainsize;
  int<lower=1> A;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> C;
  int<lower=1> S;

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

  // ---- LWL processing side (LEVEL-ONLY: rt0 + global psi) ---- //
  int<lower=0> N_lwl;
  array[N_lwl] int<lower=1, upper=I> lwl_to_child;
  vector[N_lwl] lwl_log_age;                     // log(lwl_age / a0)
  vector[N_lwl] lwl_log_rt;                       // observed log(rt_ms)

  real mu_rt_prior_mean;                          // pop log-RT level @ a0  (~6.84)
  real<lower=0> mu_rt_prior_sd;
  real psi_prior_mean;                            // global RT age-slope     (~-0.35)
  real<lower=0> mu_rtslope_prior_sd;
  real sigma_rt0_prior_mean;                      // between-child RT level SD (~0.143)
  real<lower=0> sigma_rt0_prior_sd;
  real sigma_lwl_prior_mean;                      // RT residual            (~0.24)
  real<lower=0> sigma_lwl_prior_sd;

  // ---- Observed-input channel: MEASUREMENT MODEL on RAW per-recording input ---- //
  int<lower=0> V_obs;
  array[V_obs] int<lower=1, upper=I> rec_to_child;
  vector[V_obs] log_input_obs;
  int<lower=1> n_instr;
  array[V_obs] int<lower=1, upper=n_instr> instr_of_rec;
  real sigma_r_prior_mean;
  real<lower=0> sigma_r_prior_sd;
  vector[S] mu_r_s_prior_mean;
  real<lower=0> mu_r_s_prior_sd;
  real<lower=0> sigma_meas_prior_sd;

  // ---- Ladder regression coefficients (RAW; driver sets prior SD = tau/sigma_ref
  //      for a common per-SD scale across channels). tiny sd -> pinned at 0. ---- //
  real<lower=0> gamma_in_prior_sd;               // D'0 input -> kappa
  real<lower=0> beta_xi_prior_sd;                // D'1 rt0   -> xi
  real<lower=0> beta_k0_prior_sd;                // D'2 rt0   -> kappa
}

parameters {
  // Per-child RT LEVEL latent (univariate; per-study tau_s + global psi carry the rest)
  vector[I] z_rt0;
  real<lower=0> sigma_rt0;

  // Per-child residuals (independent of input & RT)
  vector[I] log_alpha_raw;
  vector[I] zeta_raw;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_zeta;

  // Input measurement model
  vector[I] z_r;
  real<lower=0> sigma_r;
  vector[S] mu_r_s;
  vector<lower=0>[n_instr] sigma_meas;

  // Per-study RT intercept + ONE global age slope
  vector[S] tau_s;
  real psi;

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

  // Regression coefficients (ladder) -- RAW (multiply un-standardized predictors)
  real gamma_in;   // input (d_i) -> kappa
  real beta_xi;    // rt0          -> xi
  real beta_k0;    // rt0          -> kappa
}

transformed parameters {
  // Per-child RT level deviation (sum-to-zero; per-study level lives in tau_s)
  vector[I] rt0 = sigma_rt0 * z_rt0;
  rt0 = rt0 - mean(rt0);

  vector[I] log_alpha = sigma_alpha * log_alpha_raw;
  log_alpha = log_alpha - mean(log_alpha);
  vector[I] zeta = sigma_zeta * zeta_raw;
  zeta = zeta - mean(zeta);

  // Input deviation, centered WITHIN study
  vector[I] z_r_c;
  {
    vector[S] ssum = rep_vector(0, S);
    vector[S] scnt = rep_vector(0, S);
    for (i in 1:I) { ssum[study_of_child[i]] += z_r[i]; scnt[study_of_child[i]] += 1; }
    for (i in 1:I) z_r_c[i] = z_r[i] - ssum[study_of_child[i]] / scnt[study_of_child[i]];
  }
  vector[I] log_r_dev = sigma_r * z_r_c;          // d_i

  // Ability intercept and slope per child (RAW coefficients)
  vector[I] xi    = mu_r + log_r_dev + beta_xi * rt0 + log_alpha;
  vector[I] kappa = 1 + delta + gamma_in * log_r_dev + beta_k0 * rt0 + zeta;

  vector[J] delta_j;
  for (j in 1:J) delta_j[j] = mu_c[cc[j]] + tau_c[cc[j]] * delta_j_raw[j];
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);

  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + kappa[ch] * log(fmax(admin_age[a] - s, 0.01) / a0);
  }
  vector[J] item_offset;
  for (j in 1:J) item_offset[j] = beta_c[cc[j]] * log_p[j] - delta_j[j];
}

model {
  // RT level latent + priors
  z_rt0     ~ std_normal();
  sigma_rt0 ~ normal(sigma_rt0_prior_mean, sigma_rt0_prior_sd);  // truncated <lower=0>

  log_alpha_raw ~ std_normal();
  zeta_raw      ~ std_normal();
  sigma_alpha   ~ normal(0, sigma_alpha_prior_sd);
  sigma_zeta    ~ normal(0, sigma_zeta_prior_sd);

  // Input measurement-model priors
  z_r        ~ std_normal();
  sigma_r    ~ normal(sigma_r_prior_mean, sigma_r_prior_sd);
  mu_r_s     ~ normal(mu_r_s_prior_mean, mu_r_s_prior_sd);
  sigma_meas ~ normal(0, sigma_meas_prior_sd);

  // RT priors (frank_etal_2026 informative means)
  tau_s ~ normal(mu_rt_prior_mean, mu_rt_prior_sd);
  psi   ~ normal(psi_prior_mean, mu_rtslope_prior_sd);

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
  sigma_lwl ~ normal(sigma_lwl_prior_mean, sigma_lwl_prior_sd);

  // Raw regression coefficient priors (driver sets SDs = tau/sigma_ref; tiny -> pinned)
  gamma_in ~ normal(0, gamma_in_prior_sd);
  beta_xi  ~ normal(0, beta_xi_prior_sd);
  beta_k0  ~ normal(0, beta_k0_prior_sd);

  // ---- CDI likelihood (reduce_sum) ---- //
  target += reduce_sum(partial_sum_lpmf, y, grainsize,
                       aa, jj, admin_base, item_offset, lambda);

  // ---- LWL measurement likelihood (level-only: rt0 + global psi) ---- //
  if (N_lwl > 0) {
    vector[N_lwl] lwl_mean;
    for (n in 1:N_lwl) {
      int ch = lwl_to_child[n];
      int st = study_of_child[ch];
      lwl_mean[n] = tau_s[st] + rt0[ch] + psi * lwl_log_age[n];
    }
    lwl_log_rt ~ normal(lwl_mean, sigma_lwl);
  }

  // ---- Observed-input measurement likelihood ---- //
  if (V_obs > 0) {
    vector[V_obs] in_mean;
    vector[V_obs] in_sd;
    for (v in 1:V_obs) {
      int ch = rec_to_child[v];
      in_mean[v] = mu_r_s[study_of_child[ch]] + log_r_dev[ch];
      in_sd[v]   = sigma_meas[instr_of_rec[v]];
    }
    log_input_obs ~ normal(in_mean, in_sd);
  }
}

generated quantities {
  // Intercept-side variance partition (raw coefs x predictor SDs)
  real var_input_xi = square(sigma_r);
  real var_proc_xi  = square(beta_xi) * square(sigma_rt0);
  real var_resid_xi = square(sigma_alpha);
  real sigma_xi = sqrt(var_input_xi + var_proc_xi + var_resid_xi);
  real share_input_xi = var_input_xi / square(sigma_xi);
  real share_proc_xi  = var_proc_xi  / square(sigma_xi);
  real share_resid_xi = var_resid_xi / square(sigma_xi);

  // Slope-side variance partition
  real var_input_k = square(gamma_in) * square(sigma_r);
  real var_proc_k  = square(beta_k0) * square(sigma_rt0);
  real var_resid_k = square(sigma_zeta);
  real sigma_kappa = sqrt(var_input_k + var_proc_k + var_resid_k);

  // Per-SD effects (interpretable, common scale) for reporting / glmer comparison
  real eff_input_k = gamma_in * sigma_r;    // input -> kappa per 1 SD input (~glmer a_in)
  real eff_proc_xi = beta_xi  * sigma_rt0;  // rt0   -> xi    per 1 SD RT level
  real eff_proc_k  = beta_k0  * sigma_rt0;  // rt0   -> kappa per 1 SD RT level
}
