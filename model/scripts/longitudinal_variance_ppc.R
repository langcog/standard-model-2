## Proper MARGINAL posterior-predictive check for within-age variance.
##
## For each posterior draw, we sample NEW hypothetical children from the
## fitted population distribution -- MVN((mu_r, 0), (sigma_xi, sigma_zeta),
## correlation rho) -- and simulate their vocabulary at each observed
## admin age. Then compute SD of logit(vocab/J) within age bins.
##
## This is distinct from a conditional PPC, which would re-simulate the
## same observed children using their posterior-inferred abilities and
## therefore closely track the observed SDs by construction.

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(MASS)
})

fit    <- readRDS(file.path(PATHS$fits_dir, "long_2pl_slopes.rds"))
bundle <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))
draws  <- as_draws_df(fit)
sd_    <- bundle$stan_data
J <- sd_$J

# Thin to 100 draws
set.seed(42)
keep <- round(seq(1, nrow(draws), length.out = 100))
d100 <- draws[keep, ]

psi_cols <- grep("^psi\\[",    names(d100), value = TRUE)
lam_cols <- grep("^lambda\\[", names(d100), value = TRUE)

psi_mat <- as.matrix(d100[, psi_cols])
lam_mat <- as.matrix(d100[, lam_cols])
s_d      <- d100$s
del_d    <- d100$delta
sigma_xi <- d100$sigma_xi
sigma_zeta <- d100$sigma_zeta
rho      <- d100$rho_xi_zeta
mu_r     <- sd_$mu_r

log_p_obs <- sd_$log_p[sd_$jj]
jj_vec    <- sd_$jj
aa_vec    <- sd_$aa
age_per_obs <- sd_$admin_age[aa_vec]
A <- sd_$A

# For each draw d:
#   For each admin a, sample (xi_new, zeta_new) ~ MVN(...) under this draw's
#   population parameters.
#   Compute eta and simulate y_ij for that admin's items.
# Aggregate to admin-level vocab counts.

message("Marginal PPC: simulating hypothetical children at each admin age...")

admin_vocab_sims <- matrix(NA_real_,
                           nrow = nrow(d100), ncol = A)

for (d in seq_len(nrow(d100))) {
  Sigma <- matrix(c(
    sigma_xi[d]^2,                       rho[d] * sigma_xi[d] * sigma_zeta[d],
    rho[d] * sigma_xi[d] * sigma_zeta[d], sigma_zeta[d]^2
  ), 2, 2)
  new_child <- mvrnorm(A, mu = c(mu_r, 0), Sigma = Sigma)   # A x 2

  # Per-admin new xi and zeta
  xi_new_per_admin   <- new_child[, 1]
  zeta_new_per_admin <- new_child[, 2]

  # Map admin -> new child's xi/zeta
  xi_o   <- xi_new_per_admin[aa_vec]
  zeta_o <- zeta_new_per_admin[aa_vec]
  psi_o  <- psi_mat[d, jj_vec]
  lam_o  <- lam_mat[d, jj_vec]

  ae <- pmax(age_per_obs - s_d[d], 0.01)
  base <- xi_o + log_p_obs + MODEL_CONSTANTS$log_H +
          (1 + del_d[d] + zeta_o) * log(ae / MODEL_CONSTANTS$a0) - psi_o
  eta <- lam_o * base
  y_sim <- rbinom(length(eta), 1, plogis(eta))
  admin_vocab_sims[d, ] <- tapply(y_sim, aa_vec, sum)
}

logit_p <- function(v) qlogis(pmin(pmax((v + 0.5) / (J + 1), 1e-4), 1 - 1e-4))

# Observed
obs_vocab <- tapply(sd_$y, aa_vec, sum)
obs_logit <- logit_p(obs_vocab)

# Marginal sim: every draw x admin
sim_logit <- apply(admin_vocab_sims, 1, logit_p)  # A x n_draws

# SD per age bin
age_bin <- round(sd_$admin_age)
age_grp <- split(seq_along(age_bin), age_bin)
n_per_bin <- lengths(age_grp)

observed_sd <- sapply(age_grp, function(idx) {
  if (length(idx) < 5) NA else sd(obs_logit[idx])
})

sim_sd_per_draw <- sapply(age_grp, function(idx) {
  if (length(idx) < 5) rep(NA, nrow(d100))
  else apply(sim_logit[idx, , drop = FALSE], 2, sd)
})
sim_sd_med <- apply(sim_sd_per_draw, 2, median, na.rm = TRUE)
sim_sd_lo  <- apply(sim_sd_per_draw, 2, quantile, 0.025, na.rm = TRUE)
sim_sd_hi  <- apply(sim_sd_per_draw, 2, quantile, 0.975, na.rm = TRUE)

# Also: smooth theoretical curve by MC over xi/zeta population
ages_smooth <- seq(15, 31, by = 0.25)
smooth_df <- purrr::map_dfr(ages_smooth, function(a) {
  # For each draw, sample many new kids at this age, compute predicted mean
  # vocab count (expected over items), then SD of logit across kids.
  sds <- sapply(seq_len(nrow(d100)), function(d) {
    Sigma <- matrix(c(
      sigma_xi[d]^2,                       rho[d] * sigma_xi[d] * sigma_zeta[d],
      rho[d] * sigma_xi[d] * sigma_zeta[d], sigma_zeta[d]^2
    ), 2, 2)
    n_new <- 300
    xz <- mvrnorm(n_new, mu = c(mu_r, 0), Sigma = Sigma)
    xi_new <- xz[, 1]; zeta_new <- xz[, 2]
    ae <- max(a - s_d[d], 0.01)
    logL <- log(ae / MODEL_CONSTANTS$a0)
    # Expected vocab per child by summing P(y=1) over items j=1..J
    vocab_new <- sapply(seq_len(n_new), function(k) {
      eta_j <- lam_mat[d, ] * (xi_new[k] + sd_$log_p + MODEL_CONSTANTS$log_H +
                               (1 + del_d[d] + zeta_new[k]) * logL - psi_mat[d, ])
      sum(rbinom(J, 1, plogis(eta_j)))
    })
    sd(logit_p(vocab_new))
  })
  tibble::tibble(age = a,
                 sd_med = median(sds, na.rm = TRUE),
                 sd_lo = quantile(sds, 0.025, na.rm = TRUE),
                 sd_hi = quantile(sds, 0.975, na.rm = TRUE))
})

bin_df <- tibble::tibble(
  age = as.integer(names(age_grp)),
  n   = n_per_bin,
  obs_sd = observed_sd
) |> dplyr::filter(n >= 5, !is.na(obs_sd))

cat("\nObserved vs marginal PPC summary:\n")
print(bin_df, digits = 3)

p <- ggplot(smooth_df, aes(age)) +
  geom_ribbon(aes(ymin = sd_lo, ymax = sd_hi), alpha = 0.25,
              fill = "steelblue") +
  geom_line(aes(y = sd_med), colour = "steelblue", linewidth = 0.9) +
  geom_point(data = bin_df, aes(age, obs_sd, size = n),
             alpha = 0.7, show.legend = FALSE) +
  scale_size_area(max_size = 4) +
  labs(title = "Marginal PPC: within-age SD of logit(vocab/J)",
       subtitle = paste0("Blue ribbon: smooth prediction from fitted population ",
                         "distribution (new hypothetical children at each age).\n",
                         "Dots: observed SD per 1-mo bin — expected to scatter ",
                         "around the ribbon if the model is right."),
       x = "Age (months)", y = "SD of logit(vocab/J)") +
  theme_minimal(base_size = 11)

ggsave(file.path(PATHS$figs_dir, "long_2pl_slopes_ppc_variance_marginal.png"),
       p, width = 7.5, height = 5, dpi = 150)

cat("\nSaved outputs/figs/long_2pl_slopes_ppc_variance_marginal.png\n")
