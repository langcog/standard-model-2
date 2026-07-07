// M0 -- pure accumulator baseline.  kappa = 1, NO between-child variation.
// One global learning curve: every child shares the same trajectory; only item
// difficulties vary. This is the null the paper exists to falsify -- and the
// model a language model instantiates (fixed exponent, no acceleration, no
// between-learner variation).
//
//   y ~ Bernoulli_logit( mu_xi + log_H + log(age/a0) - delta_j )

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
  int<lower=1> J;
  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[N] int<lower=0, upper=1> y;
  vector[A] admin_age;
  real log_H;
  real<lower=0> a0;
  real mu_xi_prior_mean;
  real<lower=0> mu_xi_prior_sd;
  real<lower=0> tau_delta_prior_sd;
}

parameters {
  real mu_xi;
  vector[J] delta_j_raw;
  real<lower=0> tau_delta;
}

transformed parameters {
  vector[J] delta_j = tau_delta * delta_j_raw;
  delta_j = delta_j - mean(delta_j);
  vector[A] admin_base;
  for (a in 1:A) admin_base[a] = mu_xi + log_H + log(fmax(admin_age[a], 0.01) / a0);
  vector[J] item_offset = -delta_j;
}

model {
  mu_xi       ~ normal(mu_xi_prior_mean, mu_xi_prior_sd);
  delta_j_raw ~ std_normal();
  tau_delta   ~ normal(0, tau_delta_prior_sd);
  target += reduce_sum(partial_sum_lpmf, y, grainsize, aa, jj, admin_base, item_offset);
}
// LOO reconstructed in R from admin_base + item_offset (kept small; see 01_fit.R)
