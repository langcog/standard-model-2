## Model-free variability figure for slide 3b (English Wordbank WS).
##
## Within-age SD of theta_admin grows with age — the variability
## signature. Same visual style as the cross-language overlay so the
## two slides feel parallel.
##
## Output: outputs/figs/longitudinal/precursor_variability.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(lme4); library(wordbankr)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
CACHE   <- file.path(PATHS$fits_dir, "precursor", "English_American_.rds")

# Use the existing cached English fit if present; else fit fresh.
if (file.exists(CACHE)) {
  cat("Using cached English (American) precursor fit\n")
  r <- readRDS(CACHE)
  age_summary <- r$age_summary
} else {
  cat("Pulling English (American) WS data and fitting glmer ...\n")
  raw <- get_instrument_data(language = "English (American)", form = "WS",
                             administration_info = TRUE, item_info = TRUE)
  df <- raw |>
    filter(item_kind == "word", age >= 12, age <= 30) |>
    mutate(produces = as.integer(value == "produces")) |>
    filter(!is.na(produces))
  admin_idx <- df |> distinct(child_id, age) |> mutate(aa = row_number())
  item_idx <- df |> distinct(item_definition) |> mutate(jj = row_number())
  df <- df |> left_join(admin_idx, by = c("child_id","age")) |>
    left_join(item_idx, by = "item_definition")
  fit <- glmer(produces ~ 1 + (1|aa) + (1|jj), data = df, family = binomial,
               control = glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 2e5)))
  re <- ranef(fit)$aa
  theta_a <- tibble::tibble(aa = as.integer(rownames(re)),
                            theta = re[,"(Intercept)"] + fixef(fit)[1]) |>
    inner_join(admin_idx, by = "aa")
  age_summary <- theta_a |> mutate(age_int = round(age)) |>
    group_by(age_int) |>
    summarise(n = n(), sd_theta = sd(theta), .groups = "drop") |>
    filter(n >= 5)
}

# Restrict 16-30 to match cross-language slide
age_summary <- age_summary |> filter(age_int >= 16, age_int <= 30)

p <- ggplot(age_summary, aes(age_int, sd_theta)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.9,
              colour = "black", linewidth = 1.6) +
  geom_point(size = 2.5, colour = "grey25") +
  scale_x_continuous(breaks = c(16, 18, 20, 22, 24, 26, 28, 30)) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Age (months)",
       y = expression("Within-age SD of " * theta[admin] * " (logits)"),
       title = "Variability: within-age spread grows with age",
       subtitle = "Wordbank English (American) WS form. Per-admin theta from 1-PL Rasch IRT.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey25"))

ggsave(file.path(OUT_DIR, "precursor_variability.png"), p,
       width = 8, height = 5.5, dpi = 200)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_variability.png")))
