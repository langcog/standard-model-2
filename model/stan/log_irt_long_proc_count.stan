// PROC-ONLY factorization arm: the bi-lean model (log_irt_long_proc_bilingual.stan)
// with the OBSERVED-INPUT channel stripped -- the mirror image of
// log_irt_long_io_count.stan (which strips processing). Together the three form the
// step-wise factorization: enct (io+proc) / io_count (io only) / proc_count (proc only),
// all fit to the SAME bundle so the samples are identical across arms.
//
// Removed relative to bi-lean: the input measurement model (V_obs recordings,
// mu_r_s, sigma_meas, z_r/log_r_dev, sigma_r) and the input->kappa coefficient
// gamma_in. Kept identical: the Rasch vocabulary core, lexical-class difficulty
// hierarchy, LWL processing measurement (tau_s, rt0, psi), the sumscore count
// branch, and all centering conventions. xi keeps the fixed mu_r intercept (a data
// constant) so theta stays on the same scale as the other two arms; level
// identification runs through mu_c exactly as there.
//
//   vocabulary (CDI):  y ~ Bernoulli_logit( xi_i + log_H + kappa_i*log(age/a0) - delta_j )
//   processing (LWL):  log_rt ~ N( tau_s[study] + rt0_i + psi*log(age/a0), sigma_lwl )
//   xi_i    = mu_r + beta_xi*rt0_i + log_alpha_i
//   kappa_i = (1+delta) + beta_k0*rt0_i + zeta_i

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

  // ---- Ladder regression coefficient prior SDs (driver sets per-SD scale) ---- //
  real<lower=0> beta_xi_prior_sd;
  real<lower=0> beta_k0_prior_sd;

  // ---- Sumscore count likelihood (Binomial over a form's items) ---- //
  int<lower=0> n_sum;
  array[n_sum] int<lower=1, upper=I> sum_child;
  vector[n_sum] sum_log_age;
  array[n_sum] int<lower=0> sum_k;
  array[n_sum] int<lower=1> sum_form;
  int<lower=0> n_forms;
  int<lower=0> n_form_items;
  array[n_form_items] int<lower=1, upper=J> form_item;
  array[n_forms] int<lower=1> form_start;
  array[n_forms] int<lower=0> form_len;
}

parameters {
  vector[I] z_rt0;
  real<lower=0> sigma_rt0;

  vector[I] log_alpha_raw;
  vector[I] zeta_raw;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_zeta;

  vector[S] tau_s;
  real psi;

  vector[J] delta_j_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real delta;
  real<lower=0> sigma_lwl;

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

  vector[I] xi    = mu_r + beta_xi * rt0 + log_alpha;
  vector[I] kappa = 1 + delta + beta_k0 * rt0 + zeta;

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

  tau_s ~ normal(mu_rt_prior_mean, mu_rt_prior_sd);
  psi   ~ normal(psi_prior_mean, mu_rtslope_prior_sd);

  delta_j_raw ~ std_normal();
  mu_c  ~ normal(mu_mu_c, sigma_mu_c);
  tau_c ~ normal(0, 1);

  delta     ~ normal(delta_prior_mean, delta_prior_sd);
  sigma_lwl ~ normal(sigma_lwl_prior_mean, sigma_lwl_prior_sd);

  beta_xi ~ normal(0, beta_xi_prior_sd);
  beta_k0 ~ normal(0, beta_k0_prior_sd);

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

  // ---- Sumscore count likelihood: k ~ Binomial(n, mean produced rate over the form's items) ---- //
  if (n_sum > 0) {
    for (a in 1:n_sum) {
      int ch = sum_child[a];
      int nf = form_len[sum_form[a]];
      real th = xi[ch] + log_H + kappa[ch] * sum_log_age[a];
      real mu_k = 0;
      for (m in 1:nf) {
        int j = form_item[form_start[sum_form[a]] + m - 1];
        mu_k += inv_logit(th - delta_j[j]);
      }
      sum_k[a] ~ binomial(nf, mu_k / nf);
    }
  }
}

generated quantities {
  real var_proc_xi  = square(beta_xi) * square(sigma_rt0);
  real var_resid_xi = square(sigma_alpha);
  real sigma_xi = sqrt(var_proc_xi + var_resid_xi);
  real share_proc_xi  = var_proc_xi  / square(sigma_xi);
  real share_resid_xi = var_resid_xi / square(sigma_xi);

  real var_proc_k  = square(beta_k0) * square(sigma_rt0);
  real var_resid_k = square(sigma_zeta);
  real sigma_kappa = sqrt(var_proc_k + var_resid_k);

  real eff_proc_xi = beta_xi * sigma_rt0;
  real eff_proc_k  = beta_k0 * sigma_rt0;
}
