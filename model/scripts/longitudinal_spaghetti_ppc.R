## Spaghetti PPC for the longitudinal fit — MARGINAL / population-level.
##
## Left panel:  observed child trajectories (kids with >=3 admins).
## Right panel: for each of those kids' AGE SCHEDULES, we sample a
##              NEW hypothetical child from the fitted population
##              distribution (MVN over (xi, zeta)) and simulate their
##              vocabulary at those same ages. So the right panel shows
##              what "kids we might have seen" would look like under
##              the fitted model, matched to the observed age structure
##              but NOT using any real child's inferred ability.

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

## ---- 1. Observed trajectories for kids with >=3 admins
admins_per_child <- tapply(seq_len(sd_$A), sd_$admin_to_child, length)
kids_3plus <- which(admins_per_child >= 3)

obs_admin_vocab <- tapply(sd_$y, sd_$aa, sum)
obs_df <- tibble::tibble(
  aa = seq_len(sd_$A),
  ii = sd_$admin_to_child,
  age = sd_$admin_age,
  vocab = obs_admin_vocab
) |> dplyr::filter(ii %in% kids_3plus)

logit_p <- function(v) qlogis(pmin(pmax((v + 0.5) / (J + 1), 1e-4), 1 - 1e-4))
obs_df$logit_v <- logit_p(obs_df$vocab)

## ---- 2. Marginal PPC trajectories
## For each kid's age schedule, sample a new (xi, zeta) from the
## fitted population distribution under a single posterior draw, then
## simulate y_ij at each admin age. Use several posterior draws to show
## both population and posterior uncertainty. Keep the number of draws
## modest so the figure stays readable.

psi_cols <- grep("^psi\\[",    names(draws), value = TRUE)
lam_cols <- grep("^lambda\\[", names(draws), value = TRUE)

set.seed(42)
n_sim_kids_per_draw <- length(kids_3plus)   # match observed count
# Pick a single representative draw near the posterior center for the
# main display; 1 draw gives the cleanest spaghetti. (Using many draws
# piles every simulation on top of each other and obscures the picture.)
# Choose the draw closest to posterior medians on the structural params.
target <- c(sigma_xi = median(draws$sigma_xi),
            sigma_zeta = median(draws$sigma_zeta),
            rho = median(draws$rho_xi_zeta),
            s = median(draws$s),
            delta = median(draws$delta))
dist2 <- with(draws,
              (sigma_xi - target["sigma_xi"])^2 +
              (sigma_zeta - target["sigma_zeta"])^2 +
              (rho_xi_zeta - target["rho"])^2 +
              ((s - target["s"]) / 2)^2 +
              ((delta - target["delta"]) / 2)^2)
d_idx <- which.min(dist2)

psi_v  <- as.numeric(draws[d_idx, psi_cols])
lam_v  <- as.numeric(draws[d_idx, lam_cols])
s_v    <- draws$s[d_idx]
del_v  <- draws$delta[d_idx]
sx     <- draws$sigma_xi[d_idx]
sz     <- draws$sigma_zeta[d_idx]
rho    <- draws$rho_xi_zeta[d_idx]
mu_r   <- sd_$mu_r

# Precompute each real kid's age schedule and map to item indices
schedules <- obs_df |>
  dplyr::arrange(ii, age) |>
  dplyr::group_by(ii) |>
  dplyr::summarise(ages = list(age), .groups = "drop")

# Pre-build item info once
log_p_j <- sd_$log_p
# For each admin age, simulate vocab for a new hypothetical kid
Sigma <- matrix(c(sx^2, rho * sx * sz,
                  rho * sx * sz, sz^2), 2, 2)
new_kids <- mvrnorm(nrow(schedules), mu = c(mu_r, 0), Sigma = Sigma)

sim_rows <- list()
for (k in seq_len(nrow(schedules))) {
  ages_k <- schedules$ages[[k]]
  xi_k   <- new_kids[k, 1]
  zeta_k <- new_kids[k, 2]
  for (a in ages_k) {
    ae <- max(a - s_v, 0.01)
    eta_j <- lam_v * (xi_k + log_p_j + MODEL_CONSTANTS$log_H +
                      (1 + del_v + zeta_k) * log(ae / MODEL_CONSTANTS$a0) -
                      psi_v)
    vocab_k <- sum(rbinom(J, 1, plogis(eta_j)))
    sim_rows[[length(sim_rows) + 1]] <- data.frame(
      ii = schedules$ii[k], age = a,
      vocab = vocab_k, logit_v = logit_p(vocab_k))
  }
}
sim_df <- do.call(rbind, sim_rows)

## ---- 3. Plot
shared_y <- ggplot2::ylim(
  min(c(obs_df$logit_v, sim_df$logit_v)) - 0.2,
  max(c(obs_df$logit_v, sim_df$logit_v)) + 0.2
)
shared_x <- ggplot2::xlim(
  min(obs_df$age) - 0.5,
  max(obs_df$age) + 0.5
)

p_obs <- ggplot(obs_df, aes(age, logit_v, group = ii)) +
  geom_line(alpha = 0.35, colour = "black", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.8) +
  shared_y + shared_x +
  labs(title = "Observed",
       subtitle = sprintf("%d children with ≥3 admins (%d total admins)",
                          length(kids_3plus), nrow(obs_df)),
       x = "Age (months)", y = "logit(vocab / J)") +
  theme_minimal(base_size = 11)

p_sim <- ggplot(sim_df, aes(age, logit_v, group = ii)) +
  geom_line(alpha = 0.35, colour = "steelblue", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.8, colour = "steelblue") +
  shared_y + shared_x +
  labs(title = "Marginal PPC (new hypothetical children)",
       subtitle = paste0("Same # kids and same age schedules; ",
                         "new (xi, zeta) sampled from fitted MVN"),
       x = "Age (months)", y = "logit(vocab / J)") +
  theme_minimal(base_size = 11)

combined <- p_obs | p_sim
ggsave(file.path(PATHS$figs_dir,
                 "long_2pl_slopes_ppc_spaghetti_marginal.png"),
       combined, width = 12, height = 5.5, dpi = 150)

cat("\nSaved outputs/figs/long_2pl_slopes_ppc_spaghetti_marginal.png\n")
