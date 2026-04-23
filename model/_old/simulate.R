## Simulator for the baseline log-linear IRT accumulator model.
## Returns both observed data and ground-truth parameters for recovery testing.

simulate_data <- function(
  I = 200,                      # children
  J = 100,                      # words
  C = 3,                        # lexical classes
  mu_r = log(1198),             # mean log tokens/hour (Sperry/HR/WF)
  sigma_r = 0.4,                # sd of log tokens/hour (to-be-pinned)
  sigma_alpha_true = 0.5,       # sd of log learning efficiency
  mu_c_true = c(7, 8, 9),       # class-level log-threshold means
  tau_c_true = c(0.6, 0.8, 0.8),# class-level log-threshold sds
  s_true = 4.5,                 # start time (months)
  delta_true = 0,               # age rate-change exponent
  log_H = log(365),             # log waking hrs/month (12 hr/day)
  a0 = 20,                      # reference age
  age_range = c(12, 30),        # uniform age in months
  log_p_range = c(log(1e-5), log(1e-3)),  # spread of word log-probs
  seed = 42
) {
  set.seed(seed)
  stopifnot(length(mu_c_true) == C, length(tau_c_true) == C)

  # Word classes: stratified so each class has at least 5 words
  cc <- sort(rep(seq_len(C), length.out = J))
  cc <- sample(cc)

  # Word log-probabilities: uniform in log space (covers 2 orders of magnitude)
  log_p <- runif(J, log_p_range[1], log_p_range[2])

  # Word thresholds, class-hierarchical
  psi <- rnorm(J, mu_c_true[cc], tau_c_true[cc])

  # Child-level params
  log_r     <- rnorm(I, mu_r, sigma_r)
  log_alpha <- rnorm(I, 0, sigma_alpha_true)
  ages      <- runif(I, age_range[1], age_range[2])

  # Build all (child, word) observations
  idx <- expand.grid(ii = seq_len(I), jj = seq_len(J),
                     KEEP.OUT.ATTRS = FALSE)
  ae <- pmax(ages[idx$ii] - s_true, 0.01)
  eta <- log_r[idx$ii] + log_alpha[idx$ii] +
         log_p[idx$jj] + log_H +
         (1 + delta_true) * log(ae / a0) -
         psi[idx$jj]
  p_y <- plogis(eta)
  y   <- rbinom(length(p_y), 1, p_y)

  list(
    obs = data.frame(ii = idx$ii, jj = idx$jj,
                     age = ages[idx$ii],
                     log_p = log_p[idx$jj],
                     eta = eta, p_y = p_y, y = y),
    true = list(
      log_r = log_r, log_alpha = log_alpha, ages = ages,
      psi = psi, cc = cc, log_p = log_p,
      sigma_alpha = sigma_alpha_true,
      mu_c = mu_c_true, tau_c = tau_c_true,
      s = s_true, delta = delta_true,
      mu_r = mu_r, sigma_r = sigma_r
    ),
    constants = list(log_H = log_H, a0 = a0, I = I, J = J, C = C)
  )
}

## Package up into the list Stan expects.
build_stan_data <- function(sim,
                            pin_sigma_r = TRUE,
                            mu_mu_c = 8, sigma_mu_c = 3) {
  obs <- sim$obs
  list(
    N = nrow(obs),
    I = sim$constants$I,
    J = sim$constants$J,
    C = sim$constants$C,
    ii = obs$ii,
    jj = obs$jj,
    cc = sim$true$cc,
    y  = obs$y,
    age   = sim$true$ages,
    log_p = sim$true$log_p,
    log_H = sim$constants$log_H,
    a0    = sim$constants$a0,
    mu_r     = sim$true$mu_r,
    sigma_r  = if (pin_sigma_r) sim$true$sigma_r else 1.0,
    mu_mu_c  = mu_mu_c,
    sigma_mu_c = sigma_mu_c
  )
}
