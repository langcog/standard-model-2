## Data-variance pilot: per-word Chang & Bergen sigmoid slopes for two DISJOINT
## 10M-word CHILDES samples (chunk A vs chunk B), trained with identical seed,
## tokenizer, eval set, and epoch count. Data identity is the only variable.
##
## Slope = 0.434 / ParamScale (logit per natural-log experience), the same
## statistic compared to children's per-child kappa_i elsewhere in this repo.
##
## Produces:
##   outputs/figs/longitudinal/marlowe_data_variance_pilot.png
##   outputs/feng_eval/pilot_slope_summary.csv
## and prints the paired statistics quoted in outputs/marlowe_pilot_results.md.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(patchwork)
})

DATA_DIR <- "data/feng_2026"
OUT_FIG  <- "outputs/figs/longitudinal/marlowe_data_variance_pilot.png"
OUT_CSV  <- "outputs/feng_eval/pilot_slope_summary.csv"

## C&B filters: drop degenerate fits with no meaningful transition.
load_slopes <- function(path, col) {
  read_tsv(path, show_col_types = FALSE) |>
    mutate(slope = 0.434 / ParamScale,
           surp_range = ParamUpper - ParamLower) |>
    filter(is.finite(slope), ParamScale > 0.01, ParamScale < 10, surp_range > 1.0) |>
    select(Token, !!col := slope)
}

A <- load_slopes(file.path(DATA_DIR, "gpt2_childes_chunkA_seed42_sigmoids.txt"), "slope_A")
B <- load_slopes(file.path(DATA_DIR, "gpt2_childes_chunkB_seed42_sigmoids.txt"), "slope_B")
paired <- inner_join(A, B, by = "Token")
n <- nrow(paired)

pear  <- cor(paired$slope_A, paired$slope_B, method = "pearson")
spear <- cor(paired$slope_A, paired$slope_B, method = "spearman")
d     <- paired$slope_A - paired$slope_B
medA  <- median(paired$slope_A); medB <- median(paired$slope_B)

cat(sprintf("n shared words            = %d\n", n))
cat(sprintf("Pearson r(A,B)            = %.3f\n", pear))
cat(sprintf("Spearman rho(A,B)         = %.3f\n", spear))
cat(sprintf("median A / median B       = %.3f / %.3f  (marginal gap %.3f)\n",
            medA, medB, abs(medA - medB)))
cat(sprintf("paired diff (A-B)         = median %+.3f, mean %+.3f, SD %.3f\n",
            median(d), mean(d), sd(d)))
cat(sprintf("seed-only floor (24.5M)   ~ 0.01   |  kids sigma_kappa ~ 3.5\n"))

## ---- summary CSV ----
summ <- tibble(
  metric = c("n_shared", "pearson_r", "spearman_rho",
             "median_A", "median_B", "marginal_median_gap",
             "paired_diff_median", "paired_diff_mean", "paired_diff_sd"),
  value  = c(n, pear, spear, medA, medB, abs(medA - medB),
             median(d), mean(d), sd(d))
)
write_csv(summ, OUT_CSV)

## ---- Panel 1: per-word reproducibility scatter ----
lab <- sprintf("Pearson r = %.2f\nSpearman ρ = %.2f\nn = %d", pear, spear, n)
p1 <- ggplot(paired, aes(slope_A, slope_B)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_point(alpha = 0.30, size = 0.9, color = "#2c7fb8") +
  annotate("text", x = 0.05, y = 2.9, hjust = 0, vjust = 1, label = lab, size = 3.1) +
  coord_equal(xlim = c(0, 3), ylim = c(0, 3)) +
  labs(x = "slope (chunk A)", y = "slope (chunk B)") +
  theme_minimal(base_size = 11)

## ---- Panel 2: slope distributions ----
dens <- bind_rows(
  transmute(A, slope = slope_A, chunk = "chunk A"),
  transmute(B, slope = slope_B, chunk = "chunk B")
)
seed_med <- mean(c(0.72, 0.74, 0.74))  # 24.5M three-seed medians (Sherlock)
p2 <- ggplot(dens, aes(slope, fill = chunk, color = chunk)) +
  geom_density(alpha = 0.25) +
  geom_vline(xintercept = 1, color = "grey50") +
  annotate("text", x = 1.03, y = 0.02, label = "κ = 1", hjust = 0, size = 2.8) +
  geom_vline(xintercept = seed_med, linetype = "dotted") +
  annotate("text", x = seed_med - 0.03, y = 1.2, angle = 90, hjust = 1,
           label = "24.5M seed median ~0.73", size = 2.6) +
  scale_fill_manual(values = c("chunk A" = "#2c7fb8", "chunk B" = "#d95f0e")) +
  scale_color_manual(values = c("chunk A" = "#2c7fb8", "chunk B" = "#d95f0e")) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(x = "per-word slope (logit per natural-log experience)", y = "density") +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))

g <- (p1 | p2) + plot_annotation(
  title = "Data-identity variance is negligible across two disjoint 10M-word CHILDES samples",
  subtitle = paste0("Left: per-word slopes reproduce (Pearson 0.76, Spearman 0.88).  ",
                    "Right: slope distributions near-identical; both far below children (median ~10, σ_κ ≈ 3.5)."),
  theme = theme(plot.title = element_text(face = "bold", size = 12),
                plot.subtitle = element_text(size = 9, color = "grey30")))

ggsave(OUT_FIG, g, width = 10, height = 4.4, dpi = 150)
cat(sprintf("\nwrote %s\nwrote %s\n", OUT_FIG, OUT_CSV))
