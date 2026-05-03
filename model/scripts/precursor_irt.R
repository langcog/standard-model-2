## Precursor analysis: naive longitudinal IRT, no accumulator dynamics.
##
## Fits the simplest possible 1-PL Rasch model with per-admin theta_at
## (NOT per-child xi_i):
##
##     y_aj ~ Bernoulli(sigmoid(theta_a - beta_j))
##
## via glmer with random intercepts on admin and item. Extracts per-admin
## theta posterior estimates, joins with admin_info for (child, age),
## then quantifies between-child variation in two interpretable units:
##
##  1. logit-scale standard deviation of theta_a at each age bin
##  2. *months-of-development equivalent*: sigma_theta divided by the
##     population age slope of theta. Tells us "how many months ahead
##     or behind do +/-1 SD kids look".
##
## Outputs:
##   outputs/figs/longitudinal/precursor_trajectories.png  spaghetti plot of
##                                                       per-child theta_a
##                                                       trajectories.
##   outputs/figs/longitudinal/precursor_variation.png     sigma_theta vs.
##                                                       age + months-of-dev
##                                                       interpretation.

source("model/R/config.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(lme4)
})

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

bundle <- load_dataset_bundle("english")
df <- bundle$df %>% as_tibble()
admin_info <- bundle$admin_info %>% as_tibble()

cat(sprintf("Bundle: %d obs, %d admins, %d kids, %d items\n",
            nrow(df), nrow(admin_info), max(admin_info$ii),
            bundle$stan_data$J))

# ---- Fit naive longitudinal IRT ---- #
# y ~ produces, indexed by (admin aa, item jj). Per-admin random
# intercept = theta_a (ability of that admin). Per-item random
# intercept = -beta_j (negated, since glmer uses '+').
# No covariates; we want to see raw variation.
cat("\nFitting glmer y ~ 1 + (1|aa) + (1|jj) ...\n")
t0 <- Sys.time()
fit <- glmer(produces ~ 1 + (1 | aa) + (1 | jj),
             data = df, family = binomial,
             control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5)))
cat(sprintf("Fit time: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
print(summary(fit), correlation = FALSE)

# ---- Extract per-admin theta ---- #
re <- ranef(fit)$aa
theta_a <- tibble(aa = as.integer(rownames(re)),
                  theta = re[, "(Intercept)"] + fixef(fit)["(Intercept)"]) %>%
  left_join(admin_info, by = "aa")
cat(sprintf("\nExtracted theta for %d admins.\n", nrow(theta_a)))

# ---- Per-age summary ---- #
age_summary <- theta_a %>%
  mutate(age_int = round(age)) %>%
  group_by(age_int) %>%
  summarise(n = n(),
            mean_theta = mean(theta),
            sd_theta = sd(theta),
            .groups = "drop") %>%
  filter(n >= 5)

cat("\nPer-age summary (n >= 5 only):\n")
print(as.data.frame(age_summary), digits = 3, row.names = FALSE)

# ---- Population age slope of theta ---- #
slope_fit <- lm(theta ~ age, data = theta_a)
pop_slope <- coef(slope_fit)["age"]
cat(sprintf("\nPopulation age slope: %.3f logits / mo\n", pop_slope))

# Pooled sigma_theta (across all admins, ignoring age)
sigma_theta_pooled <- sd(theta_a$theta)
# Within-age sigma_theta (residual after regressing out age)
sigma_theta_within <- sd(residuals(slope_fit))

cat(sprintf("Pooled SD(theta_a):       %.3f logits\n", sigma_theta_pooled))
cat(sprintf("Within-age SD(theta_a):   %.3f logits\n", sigma_theta_within))
cat(sprintf("Months-of-dev equiv (within-age): %.1f mo\n",
            sigma_theta_within / pop_slope))
cat(sprintf("Months-of-dev equiv (pooled):     %.1f mo\n",
            sigma_theta_pooled / pop_slope))

# ---- Plot 1: spaghetti of per-child theta trajectories ---- #
# Connect admins from the same child. Highlight a few example kids.
n_kids <- length(unique(theta_a$ii))
set.seed(20260502)
example_ii <- sample(unique(theta_a$ii), 6)

theta_plot <- theta_a %>%
  mutate(highlight = ifelse(ii %in% example_ii, as.character(ii), "other"))

p_traj <- ggplot(theta_a, aes(x = age, y = theta, group = ii)) +
  geom_line(color = "gray70", alpha = 0.35, linewidth = 0.3) +
  geom_line(data = filter(theta_plot, highlight != "other"),
            aes(color = factor(highlight)),
            linewidth = 1.0, alpha = 0.9) +
  geom_point(data = filter(theta_plot, highlight != "other"),
             aes(color = factor(highlight)),
             size = 1.5, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm",
              color = "black", linewidth = 1, fill = "gray50",
              alpha = 0.18, se = TRUE) +
  scale_color_brewer(palette = "Set2", guide = "none") +
  labs(x = "age (months)",
       y = expression(theta[a] ~ "  (per-admin ability, logit)"),
       title = "(1) Per-child trajectories of theta_a (Wordbank longitudinal)",
       subtitle = sprintf("Each gray line = one of %d kids (median %d admins).  Thick colored = 6 example kids.  Black = pop. age slope.",
                          n_kids,
                          as.integer(median(table(theta_a$ii))))) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_FIGS, "precursor_trajectories.png"),
       p_traj, width = 10, height = 6, dpi = 200)
cat(sprintf("\nWrote %s\n",
            file.path(OUT_FIGS, "precursor_trajectories.png")))

# ---- Plot 2: sigma_theta vs age + months interpretation ---- #
# Build a panel showing sigma at each age, plus a translation chart.
p_sigma <- ggplot(age_summary, aes(x = age_int, y = sd_theta)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 0.7) +
  geom_hline(yintercept = sigma_theta_pooled,
             linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  annotate("text", x = max(age_summary$age_int),
           y = sigma_theta_pooled, hjust = 1, vjust = -0.3,
           label = sprintf("pooled SD = %.2f logits", sigma_theta_pooled),
           color = "firebrick", size = 3) +
  labs(x = "age (months, integer-binned)",
       y = expression("SD(" ~ theta[a] ~ ") within age bin (logits)"),
       title = "(2a) Within-age SD of theta_a vs. age",
       subtitle = sprintf("Spread of ability across kids at each age.  Constant or rising = real individual differences, not noise.")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Translation: at each age, what does +/-1 SD mean in months?
# Take the within-age SD at that bin, divide by the pop slope to get
# months-of-development equivalent for kids 1 SD above/below the mean.
age_summary_months <- age_summary %>%
  mutate(months_equiv = sd_theta / pop_slope)

p_months <- ggplot(age_summary_months, aes(x = age_int, y = months_equiv)) +
  geom_col(fill = "darkorange", alpha = 0.75, width = 0.7) +
  geom_hline(yintercept = sigma_theta_within / pop_slope,
             linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  annotate("text", x = max(age_summary_months$age_int),
           y = sigma_theta_within / pop_slope, hjust = 1, vjust = -0.3,
           label = sprintf("avg = %.1f mo / SD",
                            sigma_theta_within / pop_slope),
           color = "firebrick", size = 3) +
  labs(x = "age (months)",
       y = "months-of-development equivalent of 1 SD",
       title = "(2b) Translation: 1 SD across kids = how many months of dev?",
       subtitle = sprintf("SD(theta) / pop_slope.  Pop. slope estimated as %.3f logits/mo.",
                          pop_slope)) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_combo <- p_sigma / p_months
ggsave(file.path(OUT_FIGS, "precursor_variation.png"),
       p_combo, width = 10, height = 8, dpi = 200)
cat(sprintf("Wrote %s\n",
            file.path(OUT_FIGS, "precursor_variation.png")))

# ---- Save the per-admin theta for downstream ---- #
saveRDS(theta_a, file.path(PATHS$fits_dir, "precursor_theta_admin.rds"))
cat("\nSaved per-admin theta to fits/precursor_theta_admin.rds\n")
