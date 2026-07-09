// M3-LINEAR -- the M3 headline model with age entering LINEARLY instead of as
// log(age/a0). The Bayesian analog of the glmer ladder's D_lin (vs D_log): random
// per-child intercept + slope, correlated, differing from m3_full.stan ONLY in the
// age transform. Used for the log-vs-linear age model-form comparison (LOO).
//
// Note: there is NO "slope = 1" pure-accumulator null on the linear scale (that
// meaning is specific to log-age), so the population slope is a free `beta`
// (= delta here), not `1 + delta`. sigma_b is the SD of the per-child LINEAR slope
// (per month), on a different scale than the log model's sigma_b -- fine, the LOO
// comparison is scale-free.
//
//   (a_i, b_i) ~ MVN(0, Sigma)             correlated intercept / slope
//   xi_i    = mu_xi + a_i                  efficiency (intercept at age a0)
//   slope_i = delta + b_i                  per-child linear slope on standardized age
//   y ~ Bernoulli_logit( xi_i + log_H + slope_i * (age - a0)/sd_age - delta_j )

functions {
  real partial_sum_lpmf(array[] int y_slice, int start, int end,
                        array[] int aa, array[] int jj,
                        vector admin_base, vector item_offset) {
    int ns = end - start + 1;
    vector[ns] eta;
    for (i in 1:ns) { int o = start + i - 1; eta[i] = admin_base[aa[o]] + item_offset[jj[o]]; }
    return bernoulli_logit_lpmf(y_slice | eta);
  }
}

data {
  int<lower=1> N;
  int<lower=1> grainsize;
  int<lower=1> A;
  int<lower=1> I;
  int<lower=1> J;
  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[A] int<lower=1, upper=I> admin_to_child;
  array[N] int<lower=0, upper=1> y;
  vector[A] admin_age;
  real log_H;
  real<lower=0> a0;
  real mu_xi_prior_mean;
  real<lower=0> mu_xi_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_a_prior_sd;
  real<lower=0> sigma_b_prior_sd;
  real<lower=0> tau_delta_prior_sd;
}

transformed data {
  // Standardize the linear age predictor: (age - a0)/sd. Raw (age - a0) spans
  // ~+-10-18 months, forcing the slope + sigma_b to a tiny scale (~0.2) that
  // funnels the non-centered parameterization (severe non-convergence, worst on
  // wide-age Norwegian). Scaling puts slope/sigma_b at ~unit scale. The likelihood
  // (hence LOO/elpd) is INVARIANT to this linear rescaling -- only the geometry
  // improves. beta_pop is then the population slope per SD-of-age.
  real age_scale = sd(admin_age - a0);
  vector[A] age_z = (admin_age - a0) / age_scale;
}

parameters {
  matrix[2, I] z_child;                        // non-centered (intercept, slope)
  vector<lower=0>[2] sigma_child;              // (sigma_a, sigma_b)
  cholesky_factor_corr[2] L_child;
  real mu_xi;
  real delta;
  vector[J] delta_j_raw;
  real<lower=0> tau_delta;
}

transformed parameters {
  matrix[I, 2] child_eff = (diag_pre_multiply(sigma_child, L_child) * z_child)';
  vector[I] a_i = child_eff[, 1] - mean(child_eff[, 1]);   // mu_xi carries intercept mean
  vector[I] b_i = child_eff[, 2] - mean(child_eff[, 2]);   // delta carries slope mean
  vector[I] xi    = mu_xi + a_i;
  vector[I] slope = delta + b_i;                           // linear slope (no accumulator +1)

  vector[J] delta_j = tau_delta * delta_j_raw;
  delta_j = delta_j - mean(delta_j);

  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + slope[ch] * age_z[a];
  }
  vector[J] item_offset = -delta_j;
}

model {
  to_vector(z_child) ~ std_normal();
  sigma_child[1] ~ normal(0, sigma_a_prior_sd);
  sigma_child[2] ~ normal(0, sigma_b_prior_sd);
  L_child ~ lkj_corr_cholesky(2);
  mu_xi       ~ normal(mu_xi_prior_mean, mu_xi_prior_sd);
  delta       ~ normal(delta_prior_mean, delta_prior_sd);
  delta_j_raw ~ std_normal();
  tau_delta   ~ normal(0, tau_delta_prior_sd);
  target += reduce_sum(partial_sum_lpmf, y, grainsize, aa, jj, admin_base, item_offset);
}

generated quantities {
  real sigma_a  = sigma_child[1];
  real sigma_b  = sigma_child[2];
  real rho_ab   = multiply_lower_tri_self_transpose(L_child)[2, 1];
  real beta_pop = delta;                                   // population linear slope
}
