// Baseline log-linear IRT accumulator model
// See notes/model_explainer.pdf for full specification
//
// Linear predictor:
//   eta_ij = log r_i + log alpha_i + log p_j + log H + (1+delta)*log((a_i - s)/a0) - psi_j
// Likelihood:
//   y_ij ~ Bernoulli(inv_logit(eta_ij))   [logit link; see note in recovery_test.R]
//
// External pinning of r_i is the identifiability lever for RQ4.

data {
  int<lower=1> N;                       // total observations
  int<lower=1> I;                       // children
  int<lower=1> J;                       // words
  int<lower=1> C;                       // lexical classes

  array[N] int<lower=1, upper=I> ii;    // child index per obs
  array[N] int<lower=1, upper=J> jj;    // word index per obs
  array[J] int<lower=1, upper=C> cc;    // class index per word
  array[N] int<lower=0, upper=1> y;     // CDI response

  vector[I] age;                        // age in months
  vector[J] log_p;                      // log word probability (CHILDES)

  // fixed constants
  real log_H;                           // log waking hrs / month (~log 365)
  real<lower=0> a0;                     // reference age (months)

  // externally pinned prior on input rate
  real mu_r;
  real<lower=0> sigma_r;

  // weak hyperprior for class means
  real mu_mu_c;
  real<lower=0> sigma_mu_c;
}

parameters {
  // child-level (non-centered)
  vector[I] log_r_raw;
  vector[I] log_alpha_raw;
  real<lower=0> sigma_alpha;

  // word-level (non-centered, class-hierarchical)
  vector[J] psi_raw;
  vector[C] mu_c;
  vector<lower=0>[C] tau_c;

  // global
  real<lower=0, upper=15> s;            // start time (months)
  real delta;                            // age rate-change exponent
}

transformed parameters {
  vector[I] log_r     = mu_r + sigma_r * log_r_raw;
  vector[I] log_alpha = sigma_alpha * log_alpha_raw;
  vector[J] psi;
  for (j in 1:J) {
    psi[j] = mu_c[cc[j]] + tau_c[cc[j]] * psi_raw[j];
  }
}

model {
  // priors
  log_r_raw     ~ std_normal();
  log_alpha_raw ~ std_normal();
  sigma_alpha   ~ normal(0, 1);           // half-normal by positivity

  psi_raw ~ std_normal();
  mu_c    ~ normal(mu_mu_c, sigma_mu_c);
  tau_c   ~ normal(0, 1);                 // half-normal

  s     ~ normal(4.5, 2);                 // truncated [0,15] by decl
  delta ~ normal(0, 0.5);

  // linear predictor
  vector[N] eta;
  {
    vector[N] ae;
    for (n in 1:N) {
      ae[n] = fmax(age[ii[n]] - s, 0.01); // floor to keep log finite
    }
    eta = log_r[ii] + log_alpha[ii]
        + log_p[jj] + log_H
        + (1 + delta) * log(ae / a0)
        - psi[jj];
  }

  y ~ bernoulli_logit(eta);
}

generated quantities {
  // Variance-decomposition (RQ4): proportion of child-level variance from efficiency
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));
}
