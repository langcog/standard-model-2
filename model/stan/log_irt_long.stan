// Longitudinal log-linear IRT accumulator model.
//
// Extends log_irt.stan to properly handle longitudinal data:
// - Observations are indexed by (admin, item); each admin has its own age.
// - Children (I) are indexed by admin (A), via admin_to_child.
// - Per-child latents xi_i (intercept) and zeta_i (slope deviation) are
//   modeled with a BIVARIATE prior (LKJ on correlation), since
//   longitudinal LMM shows corr(xi, zeta) ~ 0.7.
// - 2PL discrimination lambda_j per word, hierarchical.
// - Global start time s, age rate-change exponent delta.
//
// Linear predictor:
//   eta = lambda_j * [xi_i + log p_j + log H
//                     + (1 + delta + zeta_i) * log((age_a - s)/a0) - psi_j]

data {
  int<lower=1> N;                       // observations
  int<lower=1> A;                       // admins
  int<lower=1> I;                       // children
  int<lower=1> J;                       // words
  int<lower=1> C;                       // lexical classes

  array[N] int<lower=1, upper=A> aa;    // admin index per obs
  array[N] int<lower=1, upper=J> jj;    // word index per obs
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;
  array[N] int<lower=0, upper=1> y;

  vector[A] admin_age;                  // age per admin (months)
  vector[J] log_p;

  real log_H;
  real<lower=0> a0;

  real mu_r;
  real<lower=0> sigma_r;

  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  // Toggles (near-zero SD disables feature)
  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
}

parameters {
  // Child-level bivariate (xi, zeta) using non-centered Cholesky
  matrix[2, I] z_child;                 // standardized
  vector<lower=0>[2] sigma_child;       // (sigma_xi_effective, sigma_zeta)
  cholesky_factor_corr[2] L_child;      // Cholesky of correlation

  real<lower=0> sigma_alpha;

  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real<lower=0, upper=15> s;
  real delta;

  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;
}

transformed parameters {
  // Enforce sigma_child[1]^2 = sigma_r^2 + sigma_alpha^2 so mu_r and sigma_r
  // pinning remain meaningful. We do this via a deterministic constraint:
  // don't declare sigma_xi as a parameter; compute it.
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  // Build child-level effects with the required marginal SDs.
  // Note: sigma_child[1] is NOT used as a free parameter for the xi SD
  // (we replace it with sigma_xi below), but we keep sigma_child[2] as
  // sigma_zeta freely estimated.
  matrix[I, 2] child_effs;
  {
    matrix[2, 2] L_scaled;
    L_scaled[1, 1] = sigma_xi * L_child[1, 1];
    L_scaled[1, 2] = 0;
    L_scaled[2, 1] = sigma_child[2] * L_child[2, 1];
    L_scaled[2, 2] = sigma_child[2] * L_child[2, 2];
    child_effs = (L_scaled * z_child)';
  }
  // Sum-to-zero centering on each random-effect column. Without this
  // the (delta, mean(zeta)) split is partially unidentified: the
  // likelihood only sees (1 + delta + zeta_i) so the random-effect
  // mean can absorb part of the population slope. Centering forces
  // delta to carry the full population slope and zeta_i to be
  // deviations centered at 0. Same for xi: its mean should equal mu_r.
  vector[I] xi   = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta = child_effs[, 2] - mean(child_effs[, 2]);
  real<lower=0> sigma_zeta = sigma_child[2];

  vector[J] psi;
  for (j in 1:J) {
    psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];
  }
  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);
}

model {
  // Priors
  to_vector(z_child) ~ std_normal();
  // sigma_child[1] is replaced by sigma_xi (deterministic; see
  // transformed parameters), so its raw N(0,1) prior is harmless.
  // sigma_child[2] = sigma_zeta is the actual slopes-toggle param;
  // route the variant-grammar prior here so `baseline` (tight prior
  // ~0.001) actually pins zeta off, and `slopes` (~1) frees it.
  sigma_child[1] ~ normal(0, 1);
  sigma_child[2] ~ normal(0, sigma_zeta_prior_sd);
  L_child       ~ lkj_corr_cholesky(2);  // mild prior toward 0 corr
  sigma_alpha   ~ normal(0, 1);

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);

  // Likelihood
  vector[N] eta;
  {
    vector[N] ae;
    for (n in 1:N) ae[n] = fmax(admin_age[aa[n]] - s, 0.01);
    vector[N] log_age = log(ae / a0);
    // map admin -> child for xi and zeta
    vector[N] xi_per_obs;
    vector[N] zeta_per_obs;
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      xi_per_obs[n]   = xi[ch];
      zeta_per_obs[n] = zeta[ch];
    }
    vector[N] slope_per_obs = 1 + delta + zeta_per_obs;
    vector[N] base = xi_per_obs + log_p[jj] + log_H
                   + slope_per_obs .* log_age - psi[jj];
    eta = lambda[jj] .* base;
  }
  y ~ bernoulli_logit(eta);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  real rho_xi_zeta = L_child[2, 1];  // since L_child is 2x2 Cholesky
  // posterior-expected log_alpha given xi
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
