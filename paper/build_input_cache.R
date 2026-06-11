## build_input_cache.R — data for the input-experiments triptych (Fig 3).
## Produces paper/cache/fig3_input.rds with three panels, all on one
## y-axis ("input share of between-child variance, 1 - pi"):
##   (A) input share vs AGE (continuous; intercept channel)
##   (B) ANALYTIC sensitivity of the input share to the sigma_r pin, with
##       refit anchors overlaid (see Supplemental Methods derivation)
##   (C) input share per model (EN D, NO D, IO wide-delta) vs the
##       meta-analytic band
##
## INTERCEPT CHANNEL ONLY for now: the slope/gamma channel (pi_zeta) needs
## the D' fits. EN/NO read the post-dedup D summaries; NO is currently the
## OLD pre-dedup summary as a placeholder until the salvage lands. Rebuild
## when new NO D + the D' fits arrive.

suppressPackageStartupMessages({ library(dplyr); library(here); library(tibble); library(posterior) })
CACHE <- here("paper", "cache"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

## ---- model scalars from the slim summaries ----
scal <- function(tag) {
  s <- as.data.frame(readRDS(here("fits", "summaries", paste0(tag, ".summary.rds"))))
  g <- function(v, col = "median") {
    if (!col %in% names(s)) return(NA_real_)
    i <- which(s$variable == v); if (length(i)) s[[col]][i] else NA_real_
  }
  sx <- g("sigma_xi"); sa <- g("sigma_alpha")
  list(sigma_xi = sx, sigma_alpha = sa, sigma_zeta = g("sigma_zeta"),
       rho = g("rho_xi_zeta"), sigma_r = sqrt(max(sx^2 - sa^2, 1e-6)),
       pi_alpha = g("pi_alpha"), pi_lo = g("pi_alpha", "q025"),
       pi_hi = g("pi_alpha", "q975"))
}
en <- scal("long_no_freq_slopes")             # post-dedup EN D
no <- scal("long_no_freq_slopes_norwegian")   # OLD pre-dedup placeholder

## ---- IO wide-delta (from the fit draws) ----
iofit <- readRDS(here("fits", "io_pooled_widedelta.rds"))
iod <- iofit$draws(format = "df")
io_pi <- iod$sigma_alpha^2 / (iod$sigma_alpha^2 + iod$sigma_r^2)
io <- list(sigma_xi = sqrt(median(iod$sigma_alpha)^2 + median(iod$sigma_r)^2),
           sigma_alpha = median(iod$sigma_alpha), sigma_zeta = median(iod$sigma_zeta),
           rho = 0, sigma_r = median(iod$sigma_r),
           pi_alpha = median(io_pi),
           pi_lo = unname(quantile(io_pi, 0.025)),
           pi_hi = unname(quantile(io_pi, 0.975)))

MODELS <- list(`English (D)` = en, `Norwegian (D)` = no,
               `Input-observed` = io)
MODEL_LEVELS <- names(MODELS)

a0 <- 18   # anchor age (kept for the saveRDS payload / fan helpers)

## ---- Panel A (headline): factor x channel variance partition ----
## The share of between-child variance in EFFICIENCY (xi) and ACCELERATION
## (kappa) attributable to INPUT and to PROCESSING, with everything from
## observed-input fits (no imputation):
##   xi_i    = mu_r + sigma_r z_r + beta_xi rt0 + log_alpha
##   kappa_i = (1+delta) + gamma sigma_r z_r + beta_k0 rt0 + beta_k1 rt1 + zeta
## INPUT channel  -> io_pooled_gamma_widedelta_add (4 datasets, gamma on slope).
## PROC channel   -> proc_dp3 (all betas free, so the acceleration share is
##                   ESTIMATED, not assumed 0; the selected rung D'1 agrees on xi).
## sigma_r differs between the two fits (io 0.37, proc-pinned 0.53) -- each
## factor is taken from the model best suited to it; noted in the caption.
qfac3 <- function(x) c(med = median(x), lo = quantile(x, .05, names = FALSE),
                       hi = quantile(x, .95, names = FALSE))
io_g <- as_draws_df(readRDS(here("fits", "io_pooled_gamma_widedelta_add.rds"))$draws(
          c("sigma_r", "sigma_alpha", "sigma_zeta", "gamma")))
in_eff <- io_g$sigma_r^2 / (io_g$sigma_alpha^2 + io_g$sigma_r^2)
in_acc <- (io_g$gamma^2 * io_g$sigma_r^2) / (io_g$gamma^2 * io_g$sigma_r^2 + io_g$sigma_zeta^2)
pd  <- as.data.frame(readRDS(here("fits", "summaries", "proc_dp3_all.draws.rds")))
srP <- readRDS(here("fits", "proc_dp_all_subset_data.rds"))$stan_data$sigma_r
rho <- if ("rho_rt" %in% names(pd)) pd$rho_rt else 0
Vxi  <- srP^2 + pd$beta_xi^2 * pd$sigma_rt0^2 + pd$sigma_alpha^2
Vkap <- pd$gamma_in^2 * srP^2 +
        (pd$beta_k0^2 * pd$sigma_rt0^2 + pd$beta_k1^2 * pd$sigma_rt1^2 +
         2 * pd$beta_k0 * pd$beta_k1 * rho * pd$sigma_rt0 * pd$sigma_rt1) + pd$sigma_zeta^2
pr_eff <- pd$beta_xi^2 * pd$sigma_rt0^2 / Vxi
pr_acc <- (pd$beta_k0^2 * pd$sigma_rt0^2 + pd$beta_k1^2 * pd$sigma_rt1^2 +
           2 * pd$beta_k0 * pd$beta_k1 * rho * pd$sigma_rt0 * pd$sigma_rt1) / Vkap
panel_partition <- bind_rows(
  tibble(factor = "Input",      channel = "Efficiency",   !!!qfac3(in_eff)),
  tibble(factor = "Input",      channel = "Acceleration", !!!qfac3(in_acc)),
  tibble(factor = "Processing", channel = "Efficiency",   !!!qfac3(pr_eff)),
  tibble(factor = "Processing", channel = "Acceleration", !!!qfac3(pr_acc))
) |> mutate(factor  = factor(factor,  levels = c("Input", "Processing")),
            channel = factor(channel, levels = c("Efficiency", "Acceleration")))

## ---- Panel B: analytic sensitivity (EN, NO) + refit anchors ----
## Holding sigma_xi^2 at its fitted (data-determined) value, the input
## share is the algebraic split  sigma_r^2 / sigma_xi^2.  Refit anchors
## show the small second-order drift (sigma_xi^2 shrinks slightly as the
## pin grows) -- see Supplemental Methods.
sr_grid <- seq(0.2, 1.0, by = 0.01)
panelB_curve <- bind_rows(lapply(c("English (D)", "Norwegian (D)"), function(nm) {
  p <- MODELS[[nm]]
  tibble(model = nm, sigma_r = sr_grid, input_share = sr_grid^2 / p$sigma_xi^2)
})) |> mutate(model = factor(model, levels = MODEL_LEVELS))
anch <- function(tag, sr, model, note) {
  fp <- here("fits", "summaries", paste0(tag, ".summary.rds"))
  if (!file.exists(fp)) return(NULL)
  s <- as.data.frame(readRDS(fp)); pa <- s$median[s$variable == "pi_alpha"]
  tibble(model = model, sigma_r = sr, input_share = 1 - pa, note = note)
}
panelB_anchors <- bind_rows(
  anch("long_no_freq_slopes", en$sigma_r, "English (D)", "post-dedup (on pin)"),
  anch("long_no_freq_slopes_sigmaR_0p80", 0.80, "English (D)", "pre-dedup refit"))

## Meta-analytic input-quantity share of vocabulary variance, Coffey 2026:
## 4-7% (reference band for the INPUT shares in the partition panel).
meta <- tibble(source = "Coffey 2026 (meta)", lo = 0.04, hi = 0.07, mid = 0.055)

## ---- Panels D/E: input & processing fans (model curves + per-child spaghetti) ----
## D = io_pooled input fan; E = proc_dp RT fan. Children are split into quartiles
## on the model's PREDICTOR (io log_r_dev / proc rt0 = centred log-RT) and labelled
## by percentile. Curves = model-predicted vocab at each quartile's median
## predictor; spaghetti = per-child empirical trajectories (drawn grey in the chunk).
PCT_LEVS <- c("0–25%", "25–50%", "50–75%", "75–100%")
qfac <- function(x) factor(PCT_LEVS[cut(x, quantile(x, 0:4/4), include.lowest = TRUE,
                                        labels = FALSE)], levels = PCT_LEVS)

## (E) proc_dp processing fan (D'1, the selected rung)
pb   <- readRDS(here("fits", "proc_dp_all_subset_data.rds")); psd <- pb$stan_data
pdr  <- readRDS(here("fits", "summaries", "proc_dp1_all.draws.rds"))
ppsi <- read.csv(here("fits", "summaries", "proc_dp1_all_psi.csv"))
pwi  <- pb$word_info[order(pb$word_info$jj), ]
plp  <- log(pwi$prob); pdj <- ppsi$delta_j[order(ppsi$jj)]
pbx  <- median(pdr$beta_xi); pdel <- median(pdr$delta)
prt  <- pb$lwl |> group_by(ii, dataset_name) |>
  summarise(m = mean(lwl_log_rt), .groups = "drop") |>
  group_by(dataset_name) |> mutate(rt0 = m - mean(m)) |> ungroup()   # centred log-RT = predictor
prt$q <- qfac(prt$rt0)
pages <- seq(13, 30, by = 0.5)
panelE_curves <- bind_rows(lapply(PCT_LEVS, function(L) {
  rt0q <- median(prt$rt0[prt$q == L])
  tibble(q = factor(L, levels = PCT_LEVS), age = pages,
         vocab = vapply(pages, function(t)
           mean(plogis((psd$mu_r + pbx * rt0q) + (1 + pdel) * log(t / psd$a0) + plp + psd$log_H - pdj)),
           numeric(1))) }))
panelE_spag <- pb$df |> group_by(ii, age) |> summarise(vocab = mean(produces), .groups = "drop") |>
  inner_join(prt |> select(ii, q), by = "ii")

## (D) io_pooled input fan (reuse iofit loaded above)
iob   <- readRDS(here("fits", "io_pooled_subset_data.rds")); isd <- iob$stan_data
iosum <- function(v) iofit$summary(v, ~quantile(.x, .5, names = FALSE))[[2]]
ilrd  <- iosum("log_r_dev"); idj <- iosum("delta_j")
iref  <- mean(iosum("study_ability_mean")); idel <- iosum("delta")[1]; ibc <- iosum("beta_c")[1]
iwi   <- iob$word_info[order(iob$word_info$jj), ]; ilp <- log(iwi$prob)
ich   <- tibble(ii = seq_along(ilrd), lrd = ilrd, q = qfac(ilrd))
iages <- seq(floor(min(isd$admin_age)), ceiling(max(isd$admin_age)), by = 0.5)
panelD_curves <- bind_rows(lapply(PCT_LEVS, function(L) {
  lrdq <- median(ich$lrd[ich$q == L])
  tibble(q = factor(L, levels = PCT_LEVS), age = iages,
         vocab = vapply(iages, function(t)
           mean(plogis((iref + lrdq) + (1 + idel) * log(t / isd$a0) + ibc * ilp + isd$log_H - idj)),
           numeric(1))) }))
panelD_spag <- iob$df |> group_by(ii, age) |> summarise(vocab = mean(produces), .groups = "drop") |>
  inner_join(ich |> select(ii, q), by = "ii")
cat(sprintf("Fans: proc beta_xi=%.2f; io beta_c=%.2f. proc RT quartiles n=%s\n",
            pbx, ibc, paste(as.integer(table(prt$q)), collapse = "/")))

## ---- (top row, "io-imputed") EN/NO model-IMPLIED input fans ----
## Input is imputed for EN/NO (no LENA), so we do NOT split kids by it (that
## redraws the efficiency gradient). Instead: quartile curves at
## xi = mu_r + sigma_r * z_q (normal quartile midpoints) over grey per-child
## trajectories -- the model-implied input slice. Narrative: io-imputed -> io -> proc.
zq <- qnorm(c(.125, .375, .625, .875))
mi_fan <- function(bundle, summ, psi_csv, draws, n_spag = 150, seed = 1) {
  b  <- readRDS(here("fits", bundle)); sd <- b$stan_data
  s  <- as.data.frame(readRDS(here("fits", "summaries", summ)))
  gv <- function(v) s$median[s$variable == v]
  sigma_r <- sqrt(max(gv("sigma_xi")^2 - gv("sigma_alpha")^2, 1e-6))
  kappa   <- 1 + median(as.data.frame(readRDS(here("fits", "summaries", draws)))$delta)
  psi <- read.csv(here("fits", "summaries", psi_csv))
  wi  <- b$word_info[order(b$word_info$jj), ]
  dj  <- psi$delta_j_median[match(wi$jj, psi$jj)]; dj <- dj[!is.na(dj)]
  ages <- seq(8, 30, by = 0.5)
  curves <- bind_rows(lapply(seq_along(zq), function(i) {
    xi <- sd$mu_r + sigma_r * zq[i]
    tibble(q = factor(PCT_LEVS[i], levels = PCT_LEVS), age = ages,
           vocab = vapply(ages, function(t) mean(plogis(xi + sd$log_H + kappa * log(t / sd$a0) - dj)),
                          numeric(1))) }))
  set.seed(seed)
  kids <- sample(unique(b$df$child_id), min(n_spag, n_distinct(b$df$child_id)))
  spag <- b$df |> filter(child_id %in% kids) |> group_by(child_id, age) |>
    summarise(vocab = mean(produces), .groups = "drop")
  list(curves = curves, spag = spag, sigma_r = sigma_r)
}
en_mi <- mi_fan("long_subset_data.rds", "long_no_freq_slopes.summary.rds",
                "long_no_freq_slopes_psi.csv", "long_no_freq_slopes.draws.rds")
no_mi <- mi_fan("long_subset_data_nor.rds", "long_no_freq_slopes_norwegian.summary.rds",
                "long_no_freq_slopes_norwegian_psi.csv", "long_no_freq_slopes_norwegian.draws.rds")
panelEN_curves <- en_mi$curves; panelEN_spag <- en_mi$spag
panelNO_curves <- no_mi$curves; panelNO_spag <- no_mi$spag
cat(sprintf("EN/NO model-implied input fans: sigma_r EN=%.2f NO=%.2f\n", en_mi$sigma_r, no_mi$sigma_r))

## ---- Supplement: empirical proc splits (median RT / median input), per dataset ----
## Cached here (rather than read live in the supplement) so the supplement chunk
## needs no fit bundle. Median split on the predictor; spaghetti coloured by dataset.
DL <- c(adams_marchman_2018 = "AM2018", fernald_marchman_2012 = "FM2012", fmw_2013 = "FMW2013")
svoc <- pb$df |> group_by(ii, dataset_name, age) |> summarise(vocab = mean(produces), .groups = "drop") |>
  mutate(ds = factor(DL[dataset_name], levels = DL))
srt_grp  <- prt |> transmute(ii, grp = factor(ifelse(m < median(m), "Slower processors", "Faster processors"),
                                              levels = c("Slower processors", "Faster processors")))
sinp_grp <- pb$lena |> group_by(ii) |> summarise(z = mean(z_lena), .groups = "drop") |>
  mutate(grp = factor(ifelse(z < median(z), "Lower input", "Higher input"),
                      levels = c("Lower input", "Higher input")))
panelSupp_rt    <- svoc |> inner_join(srt_grp, by = "ii")
panelSupp_input <- svoc |> inner_join(sinp_grp |> select(ii, grp), by = "ii")

saveRDS(list(panel_partition = panel_partition, meta = meta,
             panelB_curve = panelB_curve, panelB_anchors = panelB_anchors,
             panelD_curves = panelD_curves, panelD_spag = panelD_spag,
             panelE_curves = panelE_curves, panelE_spag = panelE_spag,
             panelEN_curves = panelEN_curves, panelEN_spag = panelEN_spag,
             panelNO_curves = panelNO_curves, panelNO_spag = panelNO_spag,
             panelSupp_rt = panelSupp_rt, panelSupp_input = panelSupp_input,
             sr_main = 0.53, sr_range = c(0.40, 0.70),  # plausible band: MCF to set
             a0 = a0, intercept_only = TRUE,
             no_is_placeholder = TRUE),
        file.path(CACHE, "fig3_input.rds"))
cat("Wrote paper/cache/fig3_input.rds\n\n")
cat("Headline partition (share of between-child variance):\n")
print(as.data.frame(panel_partition |> mutate(across(c(med, lo, hi), ~ sprintf("%.1f%%", 100*.x)))))
