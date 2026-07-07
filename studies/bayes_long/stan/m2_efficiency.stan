// M2 -- accelerating accumulator + between-child EFFICIENCY variation.
// Adds a per-child intercept (xi_i); acceleration is still shared (kappa = 1 + delta).
// The M1 -> M2 step adds the (obvious) fact that children differ in overall level.
//
//   xi_i = mu_xi + a_i,    a_i ~ N(0, sigma_a)
//   y ~ Bernoulli_logit( xi_i + log_H + (1 + delta) * log(age/a0) - delta_j )

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
  real<lower=0> tau_delta_prior_sd;
}

parameters {
  real mu_xi;
  real delta;
  vector[I] a_raw;
  real<lower=0> sigma_a;
  vector[J] delta_j_raw;
  real<lower=0> tau_delta;
}

transformed parameters {
  vector[I] a_i = sigma_a * a_raw;
  a_i = a_i - mean(a_i);                       // mu_xi carries the intercept mean
  vector[I] xi = mu_xi + a_i;

  vector[J] delta_j = tau_delta * delta_j_raw;
  delta_j = delta_j - mean(delta_j);

  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + (1 + delta) * log(fmax(admin_age[a], 0.01) / a0);
  }
  vector[J] item_offset = -delta_j;
}

model {
  mu_xi       ~ normal(mu_xi_prior_mean, mu_xi_prior_sd);
  delta       ~ normal(delta_prior_mean, delta_prior_sd);
  a_raw       ~ std_normal();
  sigma_a     ~ normal(0, sigma_a_prior_sd);
  delta_j_raw ~ std_normal();
  tau_delta   ~ normal(0, tau_delta_prior_sd);
  target += reduce_sum(partial_sum_lpmf, y, grainsize, aa, jj, admin_base, item_offset);
}

generated quantities {
  real kappa_pop = 1 + delta;
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = bernoulli_logit_lpmf(y[n] | admin_base[aa[n]] + item_offset[jj[n]]);
}
