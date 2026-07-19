// m_pool -- MEGA-MODEL: one accelerating-accumulator fit across ALL datasets, with
// PARTIAL POOLING of the child-level variance so sparse (two-wave) datasets borrow
// strength from the deep ones. This is the "best unified estimate" model; the
// per-dataset M0-M3 LOO ladder is kept separately as independent replication.
//
// Structure (DRAFT for MCF review -- see the design notes / open questions below):
//   * Dataset-specific MEANS (free): efficiency mu_xi[d], acceleration 1 + delta[d].
//     Datasets differ in language/instrument, so the mean levels legitimately differ.
//   * Per-child deviations (a_i, b_i) ~ MVN(0, Sigma[d]) with a DATASET-SPECIFIC
//     covariance -- but the SCALES are PARTIALLY POOLED across datasets on the log
//     scale:  log sigma_a[d] ~ N(m_a, s_a),  log sigma_b[d] ~ N(m_b, s_b).
//     This is the whole point: a sparse dataset's sigma_b[d] shrinks toward exp(m_b)
//     (set mostly by the deep datasets), instead of the likelihood over-fitting it.
//     If the sigma_b[d] come out similar, that similarity is itself a finding.
//   * Intercept-slope correlation rho[d] per dataset (LKJ) -- it genuinely varied
//     (Norwegian ~ -0.6 vs others ~ +0.1-0.3), so not pooled here.
//   * Items nested within dataset: delta_j ~ N(0, tau_delta[d(j)]). Instruments differ.
//   * Shared interpretive offsets log_H, a0.
//
//   xi_i    = mu_xi[d] + a_i
//   kappa_i = 1 + delta[d] + b_i
//   y ~ Bernoulli_logit( xi_i + log_H + kappa_i * log(age/a0) - delta_j )
//
// reduce_sum within-chain threading over observations (same as m3_full); admin_base
// and item_offset are exposed so LOO can be reconstructed in R (per-obs log_lik would
// be enormous at ~9M obs). Non-centered throughout.

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
  int<lower=1> N;                          // observations (child x age x item)
  int<lower=1> grainsize;
  int<lower=1> A;                          // administrations (child x age)
  int<lower=1> I;                          // children (pooled across datasets)
  int<lower=1> J;                          // items (pooled; each belongs to one dataset)
  int<lower=1> D;                          // datasets (5)
  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[N] int<lower=0, upper=1> y;
  array[A] int<lower=1, upper=I> admin_to_child;
  vector[A] admin_age;
  array[I] int<lower=1, upper=D> child_ds; // dataset of each child
  array[J] int<lower=1, upper=D> item_ds;  // dataset of each item
  real log_H;
  real<lower=0> a0;
  // priors (dataset means)
  real mu_xi_prior_mean;   real<lower=0> mu_xi_prior_sd;
  real delta_prior_mean;   real<lower=0> delta_prior_sd;
  // hyperpriors on the pooled log-scales (mean level + how much datasets vary)
  real log_sa_prior_mean;  real<lower=0> log_sa_prior_sd;   // m_a, prior on it
  real log_sb_prior_mean;  real<lower=0> log_sb_prior_sd;   // m_b, prior on it
  real<lower=0> s_scale_prior_sd;                            // half-N sd for s_a, s_b
  real<lower=0> tau_delta_prior_sd;
}

parameters {
  // dataset-level means
  vector[D] mu_xi;
  vector[D] delta;
  // partial-pooled child scales (log space, non-centered)
  real m_a; real<lower=0> s_a; vector[D] z_log_sa;
  real m_b; real<lower=0> s_b; vector[D] z_log_sb;
  // per-dataset intercept-slope correlation
  array[D] cholesky_factor_corr[2] L_child;
  // per-child standardized deviations (non-centered)
  matrix[2, I] z_child;
  // item difficulties (non-centered), per-dataset spread
  vector[J] delta_j_raw;
  vector<lower=0>[D] tau_delta;
}

transformed parameters {
  vector<lower=0>[D] sigma_a = exp(m_a + s_a * z_log_sa);   // pooled scales
  vector<lower=0>[D] sigma_b = exp(m_b + s_b * z_log_sb);

  vector[I] xi;
  vector[I] kappa;
  {
    for (i in 1:I) {
      int d = child_ds[i];
      // (a_i, b_i) = diag(sigma_a[d], sigma_b[d]) * L_child[d] * z_child[,i]
      vector[2] ab = diag_pre_multiply([sigma_a[d], sigma_b[d]]', L_child[d]) * z_child[, i];
      xi[i]    = mu_xi[d] + ab[1];
      kappa[i] = 1 + delta[d] + ab[2];
    }
  }

  vector[J] delta_j;
  for (j in 1:J) delta_j[j] = tau_delta[item_ds[j]] * delta_j_raw[j];

  vector[A] admin_base;
  for (a in 1:A) {
    int ch = admin_to_child[a];
    admin_base[a] = xi[ch] + log_H + kappa[ch] * log(fmax(admin_age[a], 0.01) / a0);
  }
  vector[J] item_offset = -delta_j;
}

model {
  mu_xi ~ normal(mu_xi_prior_mean, mu_xi_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);
  m_a ~ normal(log_sa_prior_mean, log_sa_prior_sd);
  m_b ~ normal(log_sb_prior_mean, log_sb_prior_sd);
  s_a ~ normal(0, s_scale_prior_sd);
  s_b ~ normal(0, s_scale_prior_sd);
  z_log_sa ~ std_normal();
  z_log_sb ~ std_normal();
  for (d in 1:D) L_child[d] ~ lkj_corr_cholesky(2);
  to_vector(z_child) ~ std_normal();
  delta_j_raw ~ std_normal();
  tau_delta ~ normal(0, tau_delta_prior_sd);

  target += reduce_sum(partial_sum_lpmf, y, grainsize, aa, jj, admin_base, item_offset);
}

generated quantities {
  // per-dataset acceleration mean + scales, and the pooled (population) hyper-means
  vector[D] kappa_pop;
  vector[D] rho_ab;
  for (d in 1:D) {
    kappa_pop[d] = 1 + delta[d];
    rho_ab[d] = multiply_lower_tri_self_transpose(L_child[d])[2, 1];
  }
  real sigma_b_pop = exp(m_b);   // pooled (typical) acceleration SD
  real sigma_a_pop = exp(m_a);
}
