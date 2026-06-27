// Longitudinal log-linear IRT + LWL processing ladder -- the "io-proc LEAN" model.
//
// Measurement core is a RASCH accumulator (common item discrimination). Three latent child
// traits, each with its own measurement model:
//   vocabulary (CDI):  y ~ Bernoulli_logit( xi_i + log_H + kappa_i*log(age/a0) - delta_j )
//   processing (LWL):  log_rt ~ N( tau_s[study] + rt0_i + psi*log(age/a0), sigma_lwl )   [level-only]
//   input (recordings): log_input ~ N( mu_r_s[study] + d_i, sigma_meas[instr] )
// Structural model linking the traits (RAW coefficients; driver sets per-SD prior SDs):
//   xi_i    = mu_r + d_i + beta_xi*rt0_i + log_alpha_i
//   kappa_i = (1+delta) + gamma_in*d_i + beta_k0*rt0_i + zeta_i
//
// Phase-0 verified (50x50 subsample, 21/21 scalars match the prior _mm lean): item
// discrimination lambda removed (was pinned ~1 -> exact Rasch), developmental start time s
// removed (was pinned ~0; the (s,delta) ridge is gone), word-frequency offset log_p removed
// (item_offset = -delta_j). delta (acceleration) prior widened by the driver. Lexical-class
// difficulty hierarchy (mu_c/tau_c, C from bundle) retained; C=1 collapse verified-equivalent
// but not adopted.

functions {
  real partial_sum_lpmf(array[] int y_slice, int start, int end,
                        array[] int aa, array[] int jj,
                        vector admin_base, vector item_offset) {
    int ns = end - start + 1;
    vector[ns] eta;
    for (i in 1:ns) {
      int o = start + i - 1;
      eta[i] = admin_base[aa[o]] + item_offset[jj[o]];   // Rasch: no lambda
    }
    return bernoulli_logit_lpmf(y_slice | eta);
  }
}

data {
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

  real log_H;
  real<lower=0> a0;

  real mu_r;
  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> sigma_alpha_prior_sd;

  // ---- LWL processing side (LEVEL-ONLY: rt0 + global psi) ---- //
  int<lower=0> N_lwl;
  array[N_lwl] int<lower=1, upper=I> lwl_to_child;
  vector[N_lwl] lwl_log_age;
  vector[N_lwl] lwl_log_rt;

  real mu_rt_prior_mean;
  real<lower=0> mu_rt_prior_sd;
  real psi_prior_mean;
  real<lower=0> mu_rtslope_prior_sd;
  real sigma_rt0_prior_mean;
  real<lower=0> sigma_rt0_prior_sd;
  real sigma_lwl_prior_mean;
  real<lower=0> sigma_lwl_prior_sd;

  // ---- Observed-input measurement model ---- //
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

  // ---- Ladder regression coefficient prior SDs (driver sets tau/sigma_ref) ---- //
  real<lower=0> gamma_in_prior_sd;
  real<lower=0> beta_xi_prior_sd;
  real<lower=0> beta_k0_prior_sd;
}

parameters {
  vector[I] z_rt0;
  real<lower=0> sigma_rt0;

  vector[I] log_alpha_raw;
  vector[I] zeta_raw;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_zeta;

  vector[I] z_r;
  real<lower=0> sigma_r;
  vector[S] mu_r_s;
  vector<lower=0>[n_instr] sigma_meas;

  vector[S] tau_s;
  real psi;

  vector[J] delta_j_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real delta;
  real<lower=0> sigma_lwl;

  real gamma_in;
  real beta_xi;
  real beta_k0;
}

transformed parameters {
  vector[I] rt0 = sigma_rt0 * z_rt0;
  rt0 = rt0 - mean(rt0);

  vector[I] log_alpha = sigma_alpha * log_alpha_raw;
  log_alpha = log_alpha - mean(log_alpha);
  vector[I] zeta = sigma_zeta * zeta_raw;
  zeta = zeta - mean(zeta);

  vector[I] z_r_c;
  {
    vector[S] ssum = rep_vector(0, S);
    vector[S] scnt = rep_vector(0, S);
    for (i in 1:I) { ssum[study_of_child[i]] += z_r[i]; scnt[study_of_child[i]] += 1; }
    for (i in 1:I) z_r_c[i] = z_r[i] - ssum[study_of_child[i]] / scnt[study_of_child[i]];
  }
  vector[I] log_r_dev = sigma_r * z_r_c;

  vector[I] xi    = mu_r + log_r_dev + beta_xi * rt0 + log_alpha;
  vector[I] kappa = 1 + delta + gamma_in * log_r_dev + beta_k0 * rt0 + zeta;

  vector[J] delta_j;
  for (j in 1:J) delta_j[j] = mu_c[cc[j]] + tau_c[cc[j]] * delta_j_raw[j];

  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + kappa[ch] * log(fmax(admin_age[a], 0.01) / a0);
  }
  vector[J] item_offset;
  for (j in 1:J) item_offset[j] = -delta_j[j];          // no frequency term
}

model {
  z_rt0     ~ std_normal();
  sigma_rt0 ~ normal(sigma_rt0_prior_mean, sigma_rt0_prior_sd);

  log_alpha_raw ~ std_normal();
  zeta_raw      ~ std_normal();
  sigma_alpha   ~ normal(0, sigma_alpha_prior_sd);
  sigma_zeta    ~ normal(0, sigma_zeta_prior_sd);

  z_r        ~ std_normal();
  sigma_r    ~ normal(sigma_r_prior_mean, sigma_r_prior_sd);
  mu_r_s     ~ normal(mu_r_s_prior_mean, mu_r_s_prior_sd);
  sigma_meas ~ normal(0, sigma_meas_prior_sd);

  tau_s ~ normal(mu_rt_prior_mean, mu_rt_prior_sd);
  psi   ~ normal(psi_prior_mean, mu_rtslope_prior_sd);

  delta_j_raw ~ std_normal();
  mu_c  ~ normal(mu_mu_c, sigma_mu_c);
  tau_c ~ normal(0, 1);

  delta     ~ normal(delta_prior_mean, delta_prior_sd);
  sigma_lwl ~ normal(sigma_lwl_prior_mean, sigma_lwl_prior_sd);

  gamma_in ~ normal(0, gamma_in_prior_sd);
  beta_xi  ~ normal(0, beta_xi_prior_sd);
  beta_k0  ~ normal(0, beta_k0_prior_sd);

  target += reduce_sum(partial_sum_lpmf, y, grainsize,
                       aa, jj, admin_base, item_offset);

  if (N_lwl > 0) {
    vector[N_lwl] lwl_mean;
    for (n in 1:N_lwl) {
      int ch = lwl_to_child[n];
      int st = study_of_child[ch];
      lwl_mean[n] = tau_s[st] + rt0[ch] + psi * lwl_log_age[n];
    }
    lwl_log_rt ~ normal(lwl_mean, sigma_lwl);
  }

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
  real var_input_xi = square(sigma_r);
  real var_proc_xi  = square(beta_xi) * square(sigma_rt0);
  real var_resid_xi = square(sigma_alpha);
  real sigma_xi = sqrt(var_input_xi + var_proc_xi + var_resid_xi);
  real share_input_xi = var_input_xi / square(sigma_xi);
  real share_proc_xi  = var_proc_xi  / square(sigma_xi);
  real share_resid_xi = var_resid_xi / square(sigma_xi);

  real var_input_k = square(gamma_in) * square(sigma_r);
  real var_proc_k  = square(beta_k0) * square(sigma_rt0);
  real var_resid_k = square(sigma_zeta);
  real sigma_kappa = sqrt(var_input_k + var_proc_k + var_resid_k);

  real eff_input_k = gamma_in * sigma_r;
  real eff_proc_xi = beta_xi  * sigma_rt0;
  real eff_proc_k  = beta_k0  * sigma_rt0;
}
