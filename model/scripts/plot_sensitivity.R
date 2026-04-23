## Plot sigma_r sensitivity results: pi_alpha vs sigma_r with CrI.
## Usage: Rscript model/scripts/plot_sensitivity.R [variant]

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args) >= 1) args[1] else "2pl"

tbl <- readRDS(file.path(PATHS$fits_dir,
                         sprintf("sensitivity_sigma_r_%s.rds", variant)))

pi_tbl <- tbl %>% filter(param == "pi_alpha")
alp_tbl <- tbl %>% filter(param == "sigma_alpha")

p1 <- ggplot(pi_tbl, aes(sigma_r, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.25,
              fill = "firebrick") +
  geom_line(colour = "firebrick", linewidth = 1) +
  geom_point(colour = "firebrick", size = 2.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = sprintf("Sensitivity of pi_alpha to sigma_r  (%s)", variant),
       subtitle = "Proportion of child-level variance from learning efficiency",
       x = expression(sigma[r]~"(pinned from external data)"),
       y = expression(pi[alpha])) +
  theme_minimal(base_size = 12)
ggsave(file.path(PATHS$figs_dir,
                 sprintf("sensitivity_pi_alpha_%s.png", variant)),
       p1, width = 6, height = 4, dpi = 150)

p2 <- ggplot(alp_tbl, aes(sigma_r, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.25,
              fill = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_point(colour = "steelblue", size = 2.5) +
  labs(title = sprintf("sigma_alpha vs sigma_r  (%s)", variant),
       subtitle = "Note sigma_xi² = sigma_r² + sigma_alpha² stays constant",
       x = expression(sigma[r]), y = expression(sigma[alpha])) +
  theme_minimal(base_size = 12)
ggsave(file.path(PATHS$figs_dir,
                 sprintf("sensitivity_sigma_alpha_%s.png", variant)),
       p2, width = 6, height = 4, dpi = 150)

cat(sprintf("Figures saved to %s/sensitivity_*_%s.png\n",
            PATHS$figs_dir, variant))
