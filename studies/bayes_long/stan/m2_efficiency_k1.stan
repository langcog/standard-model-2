// M2k1 ("M2_0") -- per-child EFFICIENCY, acceleration pinned off (kappa = 1).
//
// This is M2 with delta deleted, i.e. M0's pure-accumulator exponent combined with M2's
// per-child intercept:
//
//   xi_i = mu_xi + a_i,    a_i ~ N(0, sigma_a)
//   y ~ Bernoulli_logit( xi_i + log_H + 1 * log(age/a0) - delta_j )
//
// WHY THIS RUNG EXISTS. The forward-CV campaign so far compares M2 vs M3, which asks
// whether per-child acceleration VARIANCE forecasts -- not the paper's headline claim.
// The headline (kappa >> 1) sits on the M0 -> M1 step, but neither M0 nor M1 has a
// per-child intercept, so both must predict every child's held-out administration from
// population values alone. That comparison is dominated by between-child level variance,
// and the kappa it rewards is fitted largely to CROSS-SECTIONAL age structure -- exactly
// the confound forward CV exists to escape.
//
// Pairing this model against M2 fixes that. Both carry xi_i, so each child is anchored to
// their own level using their training administrations, and the only thing left for the
// acceleration exponent to explain is WITHIN-child growth shape. It also needs no
// per-child slope, so none of M3's identification problems apply: kappa is a population
// parameter and does not require deep per-child data to estimate.
//
// The data block is deliberately identical to m2_efficiency.stan, including the now-unused
// delta_prior_* entries, so the same bundles work as a drop-in with no harness changes.

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
  real delta_prior_mean;          // unused: kappa is pinned to 1
  real<lower=0> delta_prior_sd;   // unused: kappa is pinned to 1
  real<lower=0> sigma_a_prior_sd;
  real<lower=0> tau_delta_prior_sd;
}

parameters {
  real mu_xi;
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
    admin_base[a] = xi[ch] + log_H + log(fmax(admin_age[a], 0.01) / a0);
  }
  vector[J] item_offset = -delta_j;
}

model {
  mu_xi       ~ normal(mu_xi_prior_mean, mu_xi_prior_sd);
  a_raw       ~ std_normal();
  sigma_a     ~ normal(0, sigma_a_prior_sd);
  delta_j_raw ~ std_normal();
  tau_delta   ~ normal(0, tau_delta_prior_sd);
  target += reduce_sum(partial_sum_lpmf, y, grainsize, aa, jj, admin_base, item_offset);
}

generated quantities {
  // Constant by construction. Emitted anyway so the scoring harness's non-M3 branch --
  // which reads kappa_pop and broadcasts it across children -- works unchanged.
  real kappa_pop = 1;
}
