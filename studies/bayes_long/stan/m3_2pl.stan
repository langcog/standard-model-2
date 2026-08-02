// M3-2PL -- M3 (accelerating accumulator + per-child efficiency AND acceleration)
// with a per-item DISCRIMINATION lambda_j, i.e. a two-parameter logistic item model
// instead of the Rasch/1PL used in the main text.
//
//   (a_i, b_i) ~ MVN(0, Sigma)
//   xi_i    = mu_xi + a_i
//   kappa_i = 1 + delta + b_i
//   theta_i(t) = xi_i + log_H + kappa_i * log(t/a0)
//   y ~ Bernoulli_logit( lambda_j * ( theta_i(t) - delta_j ) )
//
// WHY. Kachergis et al. (2022, JSLHR) found the 2PL preferred over the 1PL for CDI
// data by AIC and BIC in three of four datasets, so equal discrimination is a real
// assumption to check. It is not merely cosmetic here: lambda_j MULTIPLIES theta, so it
// rescales the whole ability axis. If late-acquired words carry systematically higher
// discrimination, a 1PL has to explain the steepening vocabulary curve with a LARGER
// kappa -- trajectory misfit absorbed into the parameter of interest, the same failure
// mode that inflated sigma_b in earlier diagnostics. So this asks directly: does the
// acceleration estimate survive relaxing equal discrimination?
//
// IDENTIFICATION (the crux). The likelihood is invariant under
//     lambda -> c*lambda,   (theta - delta) -> (theta - delta)/c,
// which means kappa -> kappa/c. Note log_H being fixed data does NOT break this: the
// free intercept mu_xi simply absorbs it. So kappa is identified only relative to the
// scale of lambda, and WITHOUT a normalization a 2PL kappa is not comparable to the 1PL
// kappa at all -- and "kappa = 1 is pure accumulation" loses its meaning.
//
// We therefore pin the GEOMETRIC MEAN of lambda to exactly 1 by construction:
//     log_lambda = sigma_lambda * (log_lambda_raw - mean(log_lambda_raw))
// so mean(log lambda) = 0 identically. The 1PL is then the sigma_lambda -> 0 limit and
// kappa stays on the same scale, which is what makes the comparison interpretable.
// Location is pinned as in M3, by centring delta_j at 0 with mu_xi carrying the level.

functions {
  // eta = lambda_j * (admin_base_a - delta_j)
  //     = lambda_j * admin_base_a + item_neg_j,   item_neg_j = -lambda_j * delta_j
  // Keeping the per-item scale and offset precomputed preserves the O(A + J) then O(N)
  // structure of the 1PL models: one multiply and two lookups per observation.
  real partial_sum_lpmf(array[] int y_slice, int start, int end,
                        array[] int aa, array[] int jj,
                        vector admin_base, vector lambda, vector item_neg) {
    int ns = end - start + 1;
    vector[ns] eta;
    for (i in 1:ns) { int o = start + i - 1;
      eta[i] = lambda[jj[o]] * admin_base[aa[o]] + item_neg[jj[o]]; }
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
  real<lower=0> sigma_lambda_prior_sd;   // half-normal sd for the SD of log discrimination
}

parameters {
  matrix[2, I] z_child;
  vector<lower=0>[2] sigma_child;
  cholesky_factor_corr[2] L_child;
  real mu_xi;
  real delta;
  vector[J] delta_j_raw;
  real<lower=0> tau_delta;
  vector[J] log_lambda_raw;              // non-centred log discrimination
  real<lower=0> sigma_lambda;            // spread of log discrimination (0 => 1PL)
}

transformed parameters {
  matrix[I, 2] child_eff = (diag_pre_multiply(sigma_child, L_child) * z_child)';
  vector[I] a_i = child_eff[, 1] - mean(child_eff[, 1]);
  vector[I] b_i = child_eff[, 2] - mean(child_eff[, 2]);
  vector[I] xi    = mu_xi + a_i;
  vector[I] kappa = 1 + delta + b_i;

  vector[J] delta_j = tau_delta * delta_j_raw;
  delta_j = delta_j - mean(delta_j);

  // geometric mean of lambda fixed to 1 -- see the identification note above
  vector[J] log_lambda = sigma_lambda * (log_lambda_raw - mean(log_lambda_raw));
  vector[J] lambda = exp(log_lambda);

  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + kappa[ch] * log(fmax(admin_age[a], 0.01) / a0);
  }
  vector[J] item_neg = -lambda .* delta_j;
}

model {
  to_vector(z_child) ~ std_normal();
  sigma_child[1] ~ normal(0, sigma_a_prior_sd);
  sigma_child[2] ~ normal(0, sigma_b_prior_sd);
  L_child ~ lkj_corr_cholesky(2);
  mu_xi          ~ normal(mu_xi_prior_mean, mu_xi_prior_sd);
  delta          ~ normal(delta_prior_mean, delta_prior_sd);
  delta_j_raw    ~ std_normal();
  tau_delta      ~ normal(0, tau_delta_prior_sd);
  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);
  target += reduce_sum(partial_sum_lpmf, y, grainsize, aa, jj, admin_base, lambda, item_neg);
}

generated quantities {
  real sigma_a   = sigma_child[1];
  real sigma_b   = sigma_child[2];
  real rho_ab    = multiply_lower_tri_self_transpose(L_child)[2, 1];
  real kappa_pop = 1 + delta;
  // spread of discrimination on the natural scale; ~1 for all items means the 1PL is fine
  real lambda_sd = sd(lambda);
  real lambda_p10 = quantile(lambda, 0.10);
  real lambda_p90 = quantile(lambda, 0.90);
  // does discrimination covary with difficulty? If late/hard words discriminate more,
  // that is exactly the pattern a 1PL would have to absorb into kappa.
  real cor_lambda_delta = (mean(log_lambda .* delta_j) - mean(log_lambda) * mean(delta_j))
                          / (sd(log_lambda) * sd(delta_j));
}
