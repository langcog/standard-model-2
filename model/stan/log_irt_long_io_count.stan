// Longitudinal log-linear IRT -- the INPUT-ONLY ("io") lean model + sumscore count branch.
// = the bi-lean model with the PROCESSING (LWL) measurement model removed: no rt0, no beta_xi,
// no beta_k0, no LWL channel. This is the paper's STEP 1 (input + vocabulary); processing is
// layered in as a second measurement model in step 2 (the io-proc model). Keeps:
//   vocabulary (CDI):  y ~ Bernoulli_logit( xi_i + log_H + kappa_i*log(age/a0) - delta_j )
//   input (recordings): log_input ~ N( mu_r_s[study] + d_i, sigma_meas[instr] )
//   sumscore (count):  k ~ Binomial( n_form, mean_j inv_logit(theta_i - delta_j) )
// Structural model (input channel only):
//   xi_i    = mu_r + d_i + log_alpha_i              (input -> efficiency via the coeff-1 identity)
//   kappa_i = (1+delta) + gamma_in*d_i + zeta_i     (input -> acceleration)

functions {
  real partial_sum_lpmf(array[] int y_slice, int start, int end,
                        array[] int aa, array[] int jj,
                        vector admin_base, vector item_offset) {
    int ns = end - start + 1;
    vector[ns] eta;
    for (i in 1:ns) {
      int o = start + i - 1;
      eta[i] = admin_base[aa[o]] + item_offset[jj[o]];
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

  real<lower=0> gamma_in_prior_sd;

  // ---- Sumscore count likelihood ---- //
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
  vector[I] log_alpha_raw;
  vector[I] zeta_raw;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_zeta;

  vector[I] z_r;
  real<lower=0> sigma_r;
  vector[S] mu_r_s;
  vector<lower=0>[n_instr] sigma_meas;

  vector[J] delta_j_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real delta;
  real gamma_in;
}

transformed parameters {
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

  vector[I] xi    = mu_r + log_r_dev + log_alpha;
  vector[I] kappa = 1 + delta + gamma_in * log_r_dev + zeta;

  vector[J] delta_j;
  for (j in 1:J) delta_j[j] = mu_c[cc[j]] + tau_c[cc[j]] * delta_j_raw[j];

  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + kappa[ch] * log(fmax(admin_age[a], 0.01) / a0);
  }
  vector[J] item_offset;
  for (j in 1:J) item_offset[j] = -delta_j[j];
}

model {
  log_alpha_raw ~ std_normal();
  zeta_raw      ~ std_normal();
  sigma_alpha   ~ normal(0, sigma_alpha_prior_sd);
  sigma_zeta    ~ normal(0, sigma_zeta_prior_sd);

  z_r        ~ std_normal();
  sigma_r    ~ normal(sigma_r_prior_mean, sigma_r_prior_sd);
  mu_r_s     ~ normal(mu_r_s_prior_mean, mu_r_s_prior_sd);
  sigma_meas ~ normal(0, sigma_meas_prior_sd);

  delta_j_raw ~ std_normal();
  mu_c  ~ normal(mu_mu_c, sigma_mu_c);
  tau_c ~ normal(0, 1);

  delta    ~ normal(delta_prior_mean, delta_prior_sd);
  gamma_in ~ normal(0, gamma_in_prior_sd);

  target += reduce_sum(partial_sum_lpmf, y, grainsize, aa, jj, admin_base, item_offset);

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
  real var_input_xi = square(sigma_r);
  real var_resid_xi = square(sigma_alpha);
  real sigma_xi = sqrt(var_input_xi + var_resid_xi);
  real share_input_xi = var_input_xi / square(sigma_xi);
  real share_resid_xi = var_resid_xi / square(sigma_xi);

  real var_input_k = square(gamma_in) * square(sigma_r);
  real var_resid_k = square(sigma_zeta);
  real sigma_kappa = sqrt(var_input_k + var_resid_k);
  real eff_input_k = gamma_in * sigma_r;
}
