## Analyze a longitudinal fit.
##
## Usage:
##   Rscript model/scripts/analyze_longitudinal.R <variant> [dataset]
## Examples:
##   Rscript model/scripts/analyze_longitudinal.R long_2pl_slopes
##   Rscript model/scripts/analyze_longitudinal.R long_2pl_slopes norwegian
##
## Writes figures under model/figs/<variant>[_<dataset>]_*.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(MASS)
})

args    <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args) >= 1) args[1] else "long_2pl_slopes"
dataset <- if (length(args) >= 2) args[2] else "english"

tag     <- if (dataset == "english") { variant } else { sprintf("%s_%s", variant, dataset) }
fit_path <- file.path(PATHS$fits_dir, sprintf("%s.rds", tag))
if (!file.exists(fit_path))
  stop(sprintf("Fit not found at %s. Run the fit first.", fit_path))

fit    <- readRDS(fit_path)
bundle <- load_dataset_bundle(dataset)
draws  <- as_draws_df(fit)
sd_    <- bundle$stan_data
I <- sd_$I; A <- sd_$A; J <- sd_$J

fig_prefix <- function(n) file.path(PATHS$figs_dir, sprintf("%s_%s", tag, n))

## ---- 1. Headline scalar posteriors ----
pars <- c("sigma_alpha", "sigma_zeta", "rho_xi_zeta", "pi_alpha",
          "s", "delta", "sigma_lambda")
pars <- intersect(pars, names(draws))
scalars <- tibble::tibble(
  param  = pars,
  median = sapply(pars, function(p) median(draws[[p]])),
  lo     = sapply(pars, function(p) quantile(draws[[p]], .025)),
  hi     = sapply(pars, function(p) quantile(draws[[p]], .975))
)
cat(sprintf("\n--- Scalars (%s / %s) ---\n", variant, dataset))
print(scalars, digits = 3)

p1 <- draws[, pars] |>
  tidyr::pivot_longer(everything(), names_to = "param", values_to = "val") |>
  ggplot(aes(val)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~param, scales = "free", ncol = 3) +
  labs(title = sprintf("%s on %s: posteriors", variant, dataset),
       x = NULL, y = "density") +
  theme_minimal(base_size = 11)
ggsave(fig_prefix("1_scalars.png"), p1, width = 10, height = 5, dpi = 150)

## ---- 2. (xi, zeta) scatter ----
if ("rho_xi_zeta" %in% names(draws) && "zeta[1]" %in% names(draws)) {
  xi_cols   <- grep("^xi\\[",   names(draws), value = TRUE)
  zeta_cols <- grep("^zeta\\[", names(draws), value = TRUE)
  xz <- tibble::tibble(xi = colMeans(draws[, xi_cols]),
                       zeta = colMeans(draws[, zeta_cols]))
  r_obs <- cor(xz$xi, xz$zeta)
  p2 <- ggplot(xz, aes(xi, zeta)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, colour = "firebrick",
                linewidth = 0.6) +
    labs(title = sprintf("(xi, zeta) per child — %s / %s", variant, dataset),
         subtitle = sprintf("r(xi, zeta) = %.2f  |  posterior rho = %.2f",
                            r_obs, median(draws$rho_xi_zeta)),
         x = expression(xi[i]),
         y = expression(zeta[i])) +
    theme_minimal(base_size = 11)
  ggsave(fig_prefix("2_xi_zeta_scatter.png"), p2,
         width = 5.5, height = 5, dpi = 150)
}

## ---- 3. Observed trajectories (>=3 admins) ----
admins_per_child <- tapply(seq_len(sd_$A), sd_$admin_to_child, length)
kids_3plus <- which(admins_per_child >= 3)
obs_admin_vocab <- tapply(sd_$y, sd_$aa, sum)
obs_df <- tibble::tibble(aa = seq_len(sd_$A),
                         ii = sd_$admin_to_child,
                         age = sd_$admin_age,
                         vocab = obs_admin_vocab) |>
  dplyr::filter(ii %in% kids_3plus)
logit_p <- function(v) qlogis(pmin(pmax((v + 0.5) / (J + 1), 1e-4), 1 - 1e-4))
obs_df$logit_v <- logit_p(obs_df$vocab)

p3 <- ggplot(obs_df, aes(age, logit_v, group = ii)) +
  geom_line(alpha = 0.35, colour = "black", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.8) +
  labs(title = sprintf("Observed trajectories — %s", dataset),
       subtitle = sprintf("%d children with ≥3 admins (%d total admins)",
                          length(kids_3plus), nrow(obs_df)),
       x = "Age (months)", y = "logit(vocab / J)") +
  theme_minimal(base_size = 11)
ggsave(fig_prefix("3_observed_trajectories.png"),
       p3, width = 6, height = 5, dpi = 150)

## ---- 4. Marginal PPC variance-by-age ----
thin_keep <- round(seq(1, nrow(draws), length.out = 100))
d100 <- draws[thin_keep, ]

psi_cols <- grep("^psi\\[",    names(d100), value = TRUE)
lam_cols <- grep("^lambda\\[", names(d100), value = TRUE)

psi_mat  <- as.matrix(d100[, psi_cols])
lam_mat  <- if (length(lam_cols) > 0) { as.matrix(d100[, lam_cols]) } else { matrix(1, nrow = nrow(d100), ncol = J) }
s_d      <- d100$s
del_d    <- d100$delta
sig_xi_d <- d100$sigma_xi
sig_z_d  <- if ("sigma_zeta" %in% names(d100)) { d100$sigma_zeta } else { rep(0, nrow(d100)) }
rho_d    <- if ("rho_xi_zeta" %in% names(d100)) { d100$rho_xi_zeta } else { rep(0, nrow(d100)) }
mu_r     <- sd_$mu_r

age_per_obs <- sd_$admin_age[sd_$aa]
log_p_obs   <- sd_$log_p[sd_$jj]
jj_vec      <- sd_$jj
aa_vec      <- sd_$aa

admin_vocab_sims <- matrix(NA_real_, nrow = nrow(d100), ncol = A)
set.seed(42)
for (d in seq_len(nrow(d100))) {
  Sigma <- matrix(c(sig_xi_d[d]^2,
                    rho_d[d] * sig_xi_d[d] * sig_z_d[d],
                    rho_d[d] * sig_xi_d[d] * sig_z_d[d],
                    sig_z_d[d]^2), 2, 2)
  new_child <- mvrnorm(A, mu = c(mu_r, 0), Sigma = Sigma)
  xi_new_pa   <- new_child[, 1]
  zeta_new_pa <- new_child[, 2]
  xi_o   <- xi_new_pa[aa_vec]
  zeta_o <- zeta_new_pa[aa_vec]
  psi_o  <- psi_mat[d, jj_vec]
  lam_o  <- lam_mat[d, jj_vec]
  ae <- pmax(age_per_obs - s_d[d], 0.01)
  base <- xi_o + log_p_obs + MODEL_CONSTANTS$log_H +
          (1 + del_d[d] + zeta_o) * log(ae / MODEL_CONSTANTS$a0) - psi_o
  eta <- lam_o * base
  y_sim <- rbinom(length(eta), 1, plogis(eta))
  admin_vocab_sims[d, ] <- tapply(y_sim, aa_vec, sum)
}

obs_vocab <- tapply(sd_$y, aa_vec, sum)
obs_logit <- logit_p(obs_vocab)
sim_logit <- apply(admin_vocab_sims, 1, logit_p)

age_bin <- round(sd_$admin_age)
age_grp <- split(seq_along(age_bin), age_bin)
n_per_bin <- lengths(age_grp)
observed_sd <- sapply(age_grp, function(idx) {
  if (length(idx) < 5) NA else sd(obs_logit[idx])
})

ages_smooth <- seq(min(sd_$admin_age), max(sd_$admin_age), by = 0.25)
smooth_df <- purrr::map_dfr(ages_smooth, function(a) {
  sds <- sapply(seq_len(nrow(d100)), function(d) {
    Sigma <- matrix(c(sig_xi_d[d]^2,
                      rho_d[d] * sig_xi_d[d] * sig_z_d[d],
                      rho_d[d] * sig_xi_d[d] * sig_z_d[d],
                      sig_z_d[d]^2), 2, 2)
    n_new <- 300
    xz <- mvrnorm(n_new, mu = c(mu_r, 0), Sigma = Sigma)
    logL <- log(max(a - s_d[d], 0.01) / MODEL_CONSTANTS$a0)
    vocab_new <- sapply(seq_len(n_new), function(k) {
      eta_j <- lam_mat[d, ] * (xz[k, 1] + sd_$log_p + MODEL_CONSTANTS$log_H +
                               (1 + del_d[d] + xz[k, 2]) * logL - psi_mat[d, ])
      sum(rbinom(J, 1, plogis(eta_j)))
    })
    sd(logit_p(vocab_new))
  })
  tibble::tibble(age = a,
                 sd_med = median(sds, na.rm = TRUE),
                 sd_lo = quantile(sds, 0.025, na.rm = TRUE),
                 sd_hi = quantile(sds, 0.975, na.rm = TRUE))
})
bin_df <- tibble::tibble(age = as.integer(names(age_grp)),
                         n = n_per_bin, obs_sd = observed_sd) |>
  dplyr::filter(n >= 5, !is.na(obs_sd))

p4 <- ggplot(smooth_df, aes(age)) +
  geom_ribbon(aes(ymin = sd_lo, ymax = sd_hi), alpha = 0.25,
              fill = "steelblue") +
  geom_line(aes(y = sd_med), colour = "steelblue", linewidth = 0.9) +
  geom_point(data = bin_df, aes(age, obs_sd, size = n),
             alpha = 0.7, show.legend = FALSE) +
  scale_size_area(max_size = 4) +
  labs(title = sprintf("Marginal PPC: SD of logit(vocab/J) — %s / %s",
                       variant, dataset),
       subtitle = "Blue = model prediction for new children; dots = observed per 1-mo bin",
       x = "Age (months)", y = "SD of logit(vocab/J)") +
  theme_minimal(base_size = 11)
ggsave(fig_prefix("4_ppc_variance_marginal.png"),
       p4, width = 7.5, height = 5, dpi = 150)

## ---- 5. Marginal spaghetti PPC ----
target <- c(sigma_xi = median(draws$sigma_xi),
            sigma_zeta = if ("sigma_zeta" %in% names(draws)) { median(draws$sigma_zeta) } else { 0 },
            rho = if ("rho_xi_zeta" %in% names(draws)) { median(draws$rho_xi_zeta) } else { 0 },
            s = median(draws$s),
            delta = median(draws$delta))
dist2 <- with(draws,
              (sigma_xi - target["sigma_xi"])^2 +
              ((if ("sigma_zeta" %in% names(draws)) sigma_zeta else 0) -
                  target["sigma_zeta"])^2 +
              ((if ("rho_xi_zeta" %in% names(draws)) rho_xi_zeta else 0) -
                  target["rho"])^2 +
              ((s - target["s"]) / 2)^2 +
              ((delta - target["delta"]) / 2)^2)
d_idx <- which.min(dist2)

psi_v <- as.numeric(draws[d_idx, psi_cols])
lam_v <- if (length(lam_cols) > 0) { as.numeric(draws[d_idx, lam_cols]) } else { rep(1, J) }
s_v   <- draws$s[d_idx]
del_v <- draws$delta[d_idx]
sx    <- draws$sigma_xi[d_idx]
sz    <- if ("sigma_zeta" %in% names(draws)) { draws$sigma_zeta[d_idx] } else { 0 }
rh    <- if ("rho_xi_zeta" %in% names(draws)) { draws$rho_xi_zeta[d_idx] } else { 0 }

schedules <- obs_df |>
  dplyr::arrange(ii, age) |>
  dplyr::group_by(ii) |>
  dplyr::summarise(ages = list(age), .groups = "drop")

Sigma <- matrix(c(sx^2, rh * sx * sz, rh * sx * sz, sz^2), 2, 2)
if (sz == 0) Sigma <- matrix(c(sx^2, 0, 0, 1e-8), 2, 2)
new_kids <- mvrnorm(nrow(schedules), mu = c(mu_r, 0), Sigma = Sigma)

sim_rows <- list()
for (k in seq_len(nrow(schedules))) {
  for (a in schedules$ages[[k]]) {
    ae <- max(a - s_v, 0.01)
    eta_j <- lam_v * (new_kids[k, 1] + sd_$log_p + MODEL_CONSTANTS$log_H +
                      (1 + del_v + new_kids[k, 2]) * log(ae / MODEL_CONSTANTS$a0) -
                      psi_v)
    sim_rows[[length(sim_rows) + 1]] <- data.frame(
      ii = schedules$ii[k], age = a,
      vocab = sum(rbinom(J, 1, plogis(eta_j))))
  }
}
sim_df <- do.call(rbind, sim_rows)
sim_df$logit_v <- logit_p(sim_df$vocab)

shared_y <- ggplot2::ylim(min(c(obs_df$logit_v, sim_df$logit_v)) - 0.2,
                          max(c(obs_df$logit_v, sim_df$logit_v)) + 0.2)
shared_x <- ggplot2::xlim(min(obs_df$age) - 0.5, max(obs_df$age) + 0.5)

p_obs <- ggplot(obs_df, aes(age, logit_v, group = ii)) +
  geom_line(alpha = 0.35, colour = "black", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.8) +
  shared_y + shared_x +
  labs(title = "Observed", x = "Age (mo)", y = "logit(vocab / J)") +
  theme_minimal(base_size = 11)
p_sim <- ggplot(sim_df, aes(age, logit_v, group = ii)) +
  geom_line(alpha = 0.35, colour = "steelblue", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.8, colour = "steelblue") +
  shared_y + shared_x +
  labs(title = "Marginal PPC (new children)",
       x = "Age (mo)", y = "logit(vocab / J)") +
  theme_minimal(base_size = 11)

combined <- p_obs | p_sim
ggsave(fig_prefix("5_ppc_spaghetti_marginal.png"),
       combined, width = 12, height = 5.5, dpi = 150)

cat(sprintf("\nFigures saved under model/figs/%s_*.png\n", tag))
