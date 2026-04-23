// Alternative parameterization of log_irt.stan.
//
// The likelihood depends on xi_i = log r_i + log alpha_i only through
// its sum. The original parameterization samples both latents per child
// (2I parameters where only I are identified), which can confuse Stan's
// sampler. This version collapses to a single xi_i per child with
//   xi_i ~ N(mu_r, sigma_r^2 + sigma_alpha^2)
// so sigma_alpha is identified purely from Var(xi_i) after subtracting
// the externally-pinned sigma_r^2.

data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> C;

  array[N] int<lower=1, upper=I> ii;
  array[N] int<lower=1, upper=J> jj;
  array[J] int<lower=1, upper=C> cc;
  array[N] int<lower=0, upper=1> y;

  vector[I] age;
  vector[J] log_p;

  real log_H;
  real<lower=0> a0;

  real mu_r;
  real<lower=0> sigma_r;

  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  // Diagnostic hyperpriors: allow fixing delta or s via near-zero SD.
  // Default (baseline): s_prior=(4.5, 2), delta_prior=(0, 0.5).
  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
}

parameters {
  vector[I] xi_raw;                 // standardized child-level latent
  real<lower=0> sigma_alpha;

  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  real<lower=0, upper=15> s;
  real delta;
}

transformed parameters {
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));
  vector[I] xi = mu_r + sigma_xi * xi_raw;
  vector[J] psi;
  for (j in 1:J) {
    psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];
  }
}

model {
  xi_raw      ~ std_normal();
  sigma_alpha ~ normal(0, 1);

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);

  s     ~ normal(s_prior_mean, s_prior_sd);
  delta ~ normal(delta_prior_mean, delta_prior_sd);

  vector[N] eta;
  {
    vector[N] ae;
    for (n in 1:N) ae[n] = fmax(age[ii[n]] - s, 0.01);
    eta = xi[ii] + log_p[jj] + log_H
        + (1 + delta) * log(ae / a0) - psi[jj];
  }
  y ~ bernoulli_logit(eta);
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  // Posterior-expected log_alpha given xi (closed form).
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
  // Per-observation log-likelihood for LOO-CV.
  vector[N] log_lik;
  {
    vector[N] ae;
    for (n in 1:N) ae[n] = fmax(age[ii[n]] - s, 0.01);
    for (n in 1:N) {
      real eta_n = xi[ii[n]] + log_p[jj[n]] + log_H
                 + (1 + delta) * log(ae[n] / a0) - psi[jj[n]];
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
    }
  }
}
