## Posterior-predictive + theta-space comparison across the three pooled
## IO model variants (baseline, gamma-additive, gamma-multiplicative).
##
## Each panel: per-kid model trajectory (spaghetti, posterior MEDIAN),
## colored by study. PPC panel adds observed admin proportions as dots.
## Theta panel shows the implied latent ability trajectory.
##
## Output: figs/io/pooled_ppc_theta.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

cat("Loading bundle + fits ...\n")
b <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))
sd <- b$stan_data
fits <- list(
  baseline = readRDS(file.path(PATHS$fits_dir, "io_pooled_widedelta.rds")),
  add      = readRDS(file.path(PATHS$fits_dir, "io_pooled_gamma_widedelta_add.rds"))
  # mult excluded: wide-delta multiplicative fit had identifiability
  # collapse (Rhat up to 1.13, ess down to 20-70 for many params); the
  # tight delta prior was implicitly regularizing the gamma*A ridge.
  # See experiments.md §30.
)

# ---- extract median posterior params -----------------------------------
median_params <- function(fit, model_type) {
  d <- fit$draws(format = "df")
  pull <- function(pat) {
    cols <- grep(pat, colnames(d), value = TRUE)
    apply(d[, cols, drop = FALSE], 2, median)
  }
  p <- list(
    xi          = pull("^xi\\["),
    log_r_dev   = pull("^log_r_dev\\["),
    zeta        = pull("^zeta\\["),
    study_delta = pull("^study_delta\\["),
    lambda      = pull("^lambda\\["),
    delta_j     = pull("^delta_j\\["),
    beta_c      = pull("^beta_c\\["),
    delta       = median(d$delta),
    s           = median(d$s)
  )
  if (model_type != "baseline") p$gamma <- median(d$gamma)
  p
}

# ---- model-specific slope per kid --------------------------------------
slope_per_kid <- function(p, sd, model_type) {
  s_of_c <- sd$study_of_child
  if (model_type == "baseline") {
    1 + p$delta + p$study_delta[s_of_c] + p$zeta
  } else if (model_type == "add") {
    1 + p$delta + p$study_delta[s_of_c] + p$gamma * p$log_r_dev + p$zeta
  } else {                                          # mult
    (1 + p$delta + p$study_delta[s_of_c] + p$zeta) * (1 + p$gamma * p$log_r_dev)
  }
}

# ---- predictions over age grid -----------------------------------------
pred_for_model <- function(fit, model_type, sd, age_grid) {
  cat(sprintf("  computing %s ...\n", model_type))
  p <- median_params(fit, model_type)
  slope_c  <- slope_per_kid(p, sd, model_type)
  item_off <- p$beta_c[sd$cc] * sd$log_p - p$delta_j        # length J
  out <- expand.grid(ii = seq_len(sd$I), age = age_grid)
  out$pred  <- NA_real_
  out$theta <- NA_real_
  for (r in seq_len(nrow(out))) {
    i  <- out$ii[r]; a <- out$age[r]
    ae <- max(a - p$s, 0.01); la <- log(ae / sd$a0)
    theta_i <- p$xi[i] + slope_c[i] * la
    out$theta[r] <- theta_i
    eta <- p$lambda * (p$xi[i] + sd$log_H + slope_c[i] * la + item_off)
    out$pred[r]  <- mean(plogis(eta))
  }
  out$model <- model_type
  out
}

age_grid <- seq(8, 30, by = 0.5)
preds <- bind_rows(lapply(names(fits),
                          function(m) pred_for_model(fits[[m]], m, sd, age_grid)))

# attach study labels (1=BabyView, 2=SEEDLingS, 3=AM2018, 4=FMW2013)
study_map <- setNames(b$studies$study, b$studies$s)
preds$study <- factor(study_map[as.character(b$child_info$s[preds$ii])],
                      levels = b$studies$study)

# observed admin-level proportions
obs <- b$df |>
  group_by(study, ii, aa, age) |>
  summarise(prop = mean(produces), .groups = "drop") |>
  mutate(study = factor(study, levels = b$studies$study))

# ---- plot --------------------------------------------------------------
model_labels <- c(baseline = "baseline (no γ)",
                  add      = "γ additive")
preds$model <- factor(preds$model, levels = names(model_labels), labels = model_labels)
study_pal <- c(BabyView = "#E69F00", SEEDLingS = "#56B4E9",
               AM2018   = "#009E73", FMW2013   = "#D55E00")

ppc <- ggplot() +
  geom_point(data = obs, aes(age, prop, color = study), alpha = 0.45, size = 0.9) +
  geom_line(data = preds, aes(age, pred, color = study, group = interaction(model, ii)),
            alpha = 0.22, linewidth = 0.35) +
  facet_grid(. ~ model) +
  scale_color_manual(values = study_pal) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "age (months)", y = "vocabulary proportion produced",
       color = "study", title = "Posterior-predictive curves vs. observed admins") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

theta <- ggplot(preds, aes(age, theta, color = study, group = interaction(model, ii))) +
  geom_line(alpha = 0.22, linewidth = 0.35) +
  facet_grid(. ~ model) +
  scale_color_manual(values = study_pal) +
  labs(x = "age (months)", y = expression(theta == xi + kappa[i]*log(age/a[0])),
       color = "study", title = "Implied θ trajectories (latent ability)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

out_path <- file.path("figs", "io", "pooled_ppc_theta.png")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(out_path, ppc / theta, width = 13, height = 8.5, dpi = 150)
cat(sprintf("Wrote %s\n", out_path))
