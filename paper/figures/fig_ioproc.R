## io-proc figure (single column, single panel): the direct input + processing
## study's coefficients, laid out like the demographics meta panel so the two
## read the same way: rows = predictor (Input / Processing, strips tinted with
## CHANNEL_PAL), columns = trait (Efficiency / Acceleration, strips tinted
## with TRAIT_PAL, points in the trait hue), and within each cell the main
## analysis above its check -- JOINT model (filled circle) over the channel
## fit ALONE (io-only / proc-only arm, open diamond), exactly as
## cross-sectional / longitudinal in the demographics figure. All three arms
## fit the same english_count bundle, so the comparison is exact. Per-SD
## effects in the model's own units for both traits (logits for efficiency,
## kappa units for acceleration -- the same units as the demographics figure
## and the text; multiply acceleration by log(30/21) = 0.36 for the
## level-equivalent effect at 30 months). 90% CIs. No panel letter.
## The +/-2 SD trajectory demo that used to be panel A is panel B of the
## graphical-model figure now (fig_model_trajectories.R); the imputed-
## population sweep lives in fig_imputed_meta.R.
##
## Reads: paper/cache/fig_io_imputed_proc.rds (panelB_arms, built by
##        paper/build_fig_io_cache.R).
## Run:   Rscript paper/figures/fig_ioproc.R
source(here::here("paper", "figures", "_common.R"))

io_fig <- readRDS(file.path(CACHE, "fig_io_imputed_proc.rds"))

SHP <- c("Joint model" = 16, "Separate model" = 23)     # main = filled circle, check = open diamond
pb <- io_fig$panelB_arms |>
  mutate(component = factor(channel, levels = names(TRAIT_PAL)),
         predictor = factor(factor, levels = names(CHANNEL_PAL)),
         method    = factor(model, levels = names(SHP)),
         row       = factor(method, levels = rev(names(SHP))))   # joint on top
p <- ggplot(pb, aes(med, row, colour = component, shape = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.45) +
  geom_point(fill = "white", size = 2, stroke = 0.7) +
  scale_colour_manual(values = TRAIT_PAL, guide = "none") +
  scale_shape_manual(values = SHP, guide = "none") +
  ggh4x::facet_grid2(predictor ~ component, strip = trait_strips(row_pal = CHANNEL_PAL)) +
  labs(x = "effect of +1 SD (logits), 90% CI", y = NULL) +
  theme_fig()

save_fig(p, "ioproc", width = PNAS_1COL, height = PNAS_1COL * 0.78)
