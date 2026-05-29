// Pooled hierarchical IO model + MULTIPLICATIVE input-on-slope term (gamma).
//
// Identical to log_irt_io_pooled.stan EXCEPT the per-child slope is SCALED
// by the input deviation: a child's input deviation modulates the whole
// accumulation rate proportionally (a "bootstrapping" form). This is the
// multiplicative companion to log_irt_io_pooled_gamma_add.stan.
//
//   slope_i = (1 + delta + study_delta + zeta_i) * (1 + gamma * log_r_dev_i)
//
// gamma is a PROPORTION (fractional slope change per unit log_r_dev), so its
// natural scale is ~ gamma_add / kappa (kappa ~ 10). The runner sets a
// correspondingly smaller gamma_prior_sd. zeta_i is the multiplicative-base
// slope deviation; the input modulation acts on top of it.
//
// Identification: as in the additive form, the intercept (level) vs slope
// (age-interaction) signatures separate gamma from the intercept's input
// term. Watch the gamma x sigma_r posterior correlation.

functions {
  // Partial log-likelihood over a slice [start:end] of CDI observations.
  real cdi_partial_lpmf(array[] int y_slice, int start, int end,
                        array[] int aa, array[] int jj,
                        array[] int admin_to_child, vector admin_age,
                        vector xi, vector slope_child, vector item_offset,
                        vector lambda, real s, real a0, real log_H) {
    int M = end - start + 1;
    vector[M] eta;
    for (i in 1:M) {
      int n  = start + i - 1;
      int a  = aa[n];
      int ch = admin_to_child[a];
      int j  = jj[n];
      real ae = fmax(admin_age[a] - s, 0.01);
      real base = xi[ch] + log_H + slope_child[ch] * log(ae / a0)
                  + item_offset[j];
      eta[i] = lambda[j] * base;
    }
    return bernoulli_logit_lpmf(y_slice | eta);
  }
}

data {
  int<lower=1> N;            // CDI observations
  int<lower=1> A;            // admins
  int<lower=1> I;            // children
  int<lower=1> J;            // items
  int<lower=1> C;            // lexical classes
  int<lower=1> S;            // studies
  int<lower=1> n_meas;       // measurement types (head-cam, LENA)

  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[A] int<lower=1, upper=I> admin_to_child;
  array[I] int<lower=1, upper=S> study_of_child;
  array[S] int<lower=1, upper=n_meas> meas_of_study;
  array[J] int<lower=1, upper=C> cc;
  array[N] int<lower=0, upper=1> y;

  vector[A] admin_age;
  vector[J] log_p;
  real log_H;
  real<lower=0> a0;

  // input channel
  int<lower=1> V;
  array[V] int<lower=1, upper=I> recording_to_child;
  array[V] int<lower=1, upper=S> study_of_recording;
  vector[V] log_r_obs;

  // priors
  real mu_r_prior_mean;
  real<lower=0> mu_r_prior_sd;
  real<lower=0> sigma_r_prior_sd;
  real<lower=0> sigma_within_prior_sd;
  real<lower=0> sigma_study_delta_prior_sd;
  real s_prior_mean;
  real<lower=0> s_prior_sd;
  real delta_prior_mean;
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> beta_c_prior_sd;
  real beta_c_prior_mean;

  // delta_j anchor
  vector[J] delta_j_prior_mean;
  vector<lower=0>[J] delta_j_prior_sd;

  // reduce_sum partition size for the CDI likelihood
  int<lower=1> grainsize;

  // input-on-slope prior (gamma is a fractional/proportional modulation)
  real<lower=0> gamma_prior_sd;
}

parameters {
  // study-level intercepts (free; identified by data + anchored delta_j)
  vector[S] study_input_mean;     // per-study input level
  vector[S] study_ability_mean;   // per-study vocab/ability level

  // per-child latents (non-centered)
  vector[I] log_r_dev_raw;        // input deviation, SD sigma_r
  vector[I] log_alpha_raw;        // efficiency, SD sigma_alpha
  vector[I] zeta_raw;             // slope deviation, SD sigma_zeta
  real<lower=0> sigma_r;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_zeta;

  // acceleration
  real delta;                     // shared population acceleration
  vector[S] study_delta_raw;      // per-study deviation (sum-to-zero)
  real<lower=0> sigma_study_delta;
  real gamma;                     // input-on-slope modulation (shared)

  // word level
  vector[J] delta_j;              // anchored via prior
  vector[J] log_lambda_raw;
  real<lower=0> sigma_lambda;
  vector[C] beta_c;

  // global onset
  real<lower=0, upper=15> s;

  // measurement noise by type
  vector<lower=0>[n_meas] sigma_within;
}

transformed parameters {
  vector[I] log_r_dev = sigma_r * log_r_dev_raw;
  vector[I] log_alpha = sigma_alpha * log_alpha_raw;
  vector[I] zeta      = sigma_zeta * zeta_raw;
  // per-study delta deviation, sum-to-zero so `delta` is the mean
  vector[S] study_delta = study_delta_raw - mean(study_delta_raw);

  vector[J] log_lambda = sigma_lambda * log_lambda_raw;
  vector[J] lambda = exp(log_lambda);

  // ability intercept: study level + shared input deviation + efficiency
  vector[I] xi;
  for (i in 1:I)
    xi[i] = study_ability_mean[study_of_child[i]] + log_r_dev[i] + log_alpha[i];
}

model {
  // priors
  study_input_mean   ~ normal(mu_r_prior_mean, mu_r_prior_sd);
  study_ability_mean ~ normal(0, 5);
  log_r_dev_raw ~ std_normal();
  log_alpha_raw ~ std_normal();
  zeta_raw      ~ std_normal();
  sigma_r     ~ normal(0, sigma_r_prior_sd);
  sigma_alpha ~ normal(0, 1);
  sigma_zeta  ~ normal(0, sigma_zeta_prior_sd);

  delta            ~ normal(delta_prior_mean, delta_prior_sd);
  study_delta_raw  ~ normal(0, sigma_study_delta);
  sigma_study_delta ~ normal(0, sigma_study_delta_prior_sd);
  gamma            ~ normal(0, gamma_prior_sd);

  delta_j        ~ normal(delta_j_prior_mean, delta_j_prior_sd);
  log_lambda_raw ~ std_normal();
  sigma_lambda   ~ normal(0, sigma_lambda_prior_sd);
  beta_c         ~ normal(beta_c_prior_mean, beta_c_prior_sd);
  s              ~ normal(s_prior_mean, s_prior_sd);
  sigma_within   ~ normal(0, sigma_within_prior_sd);

  // input channel: observed log rate = study mean + kid deviation + noise
  {
    vector[V] mu_obs;
    vector[V] sd_obs;
    for (v in 1:V) {
      mu_obs[v] = study_input_mean[study_of_recording[v]]
                  + log_r_dev[recording_to_child[v]];
      sd_obs[v] = sigma_within[meas_of_study[study_of_recording[v]]];
    }
    log_r_obs ~ normal(mu_obs, sd_obs);
  }

  // CDI likelihood (threaded over observations via reduce_sum)
  {
    // per-child slope (MULTIPLICATIVE input term) and per-item offset
    vector[I] slope_child;
    for (i in 1:I)
      slope_child[i] = (1 + delta + study_delta[study_of_child[i]] + zeta[i])
                       * (1 + gamma * log_r_dev[i]);
    vector[J] item_offset;
    for (j in 1:J)
      item_offset[j] = beta_c[cc[j]] * log_p[j] - delta_j[j];

    target += reduce_sum(cdi_partial_lpmf, y, grainsize,
                         aa, jj, admin_to_child, admin_age,
                         xi, slope_child, item_offset, lambda,
                         s, a0, log_H);
  }
}

generated quantities {
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
  // per-study acceleration at MEAN input (log_r_dev = 0): 1 + delta + study_delta
  vector[S] kappa_study;
  for (st in 1:S) kappa_study[st] = 1 + delta + study_delta[st];
  real kappa_pop = 1 + delta;
}
