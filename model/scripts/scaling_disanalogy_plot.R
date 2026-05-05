## Scaling disanalogy figure: kids (SM2) vs LLM scaling laws.
##
## Two-panel composition:
##   D1_scaling_disanalogy.png
##     Left panel  — log-log "scaling diagram." Per-child trajectories
##                   (slope = 1 + delta + zeta_i, drawn from posterior of
##                   the named fit) over their D range, alongside Chinchilla
##                   data-bound and compute-optimal lines over the LLM D
##                   range. Slopes annotated.
##     Right panel — distribution of per-instance scaling exponents:
##                   posterior over (1 + delta + zeta_i) marginal across
##                   draws and kids, vs. LLM exponents as point lines.
##
## Pulls global posterior summaries from
##   fits/summaries/<fit_name>.draws.rds
## and the externally pinned input-rate constants from
##   fits/<stan_data_name>.rds$stan_data
##
## Defaults: long_slopes (English n=200). Pass a different fit_name to
## re-render on production fits without changing the script.
##
## Run:
##   Rscript model/scripts/scaling_disanalogy_plot.R                # default
##   Rscript model/scripts/scaling_disanalogy_plot.R long_2pl_slopes
##
## Output: outputs/figs/schematic/D1_scaling_disanalogy.png
##         outputs/figs/schematic/D1_scaling_disanalogy_data.rds
##         (the simulated trajectories + LLM reference data, for re-styling
##          without re-simulating)

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PATHS$figs_dir, "schematic")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- LLM scaling-law reference values -------------------------------
## Chinchilla (Hoffmann et al. 2022) data-bound exponent: ~0.28
## Chinchilla compute-optimal exponent on aggregate loss: ~0.155
## Kaplan et al. (2020) data-bound exponent: ~0.095
## (These are population values; within-architecture dispersion is ~0.)
LLM_REF <- list(
  chinchilla_data    = list(beta = 0.28,  label = "Chinchilla (data-bound)",
                            colour = "#1f78b4"),
  chinchilla_compute = list(beta = 0.155, label = "Chinchilla (compute-optimal)",
                            colour = "#a6cee3"),
  kaplan_data        = list(beta = 0.095, label = "Kaplan (data-bound)",
                            colour = "#999999")
)

# Typical LLM scaling-law D range (in log10): ~10^9 to ~10^12
LLM_D_RANGE_LOG10 <- c(9, 12.3)

# ---- Pull fit parameters --------------------------------------------
get_kid_scaling_params <- function(
  fit_name = "long_slopes",
  summaries_dir = file.path(PATHS$fits_dir, "summaries"),
  stan_data_paths = list(
    long_slopes      = file.path(PATHS$fits_dir, "long_subset_data.rds"),
    long_2pl_slopes  = file.path(PATHS$fits_dir, "long_subset_data.rds"),
    long_proc_slopes = file.path(PATHS$fits_dir, "stanford_linked_subset_data.rds"),
    io_slopes        = file.path(PATHS$fits_dir, "babyview_subset_data.rds"),
    io_slopes_seedlings = file.path(PATHS$fits_dir, "seedlings_subset_data.rds")
  )
) {
  draws_path <- file.path(summaries_dir, paste0(fit_name, ".draws.rds"))
  if (!file.exists(draws_path)) {
    stop("Posterior draws not found: ", draws_path,
         "\nAvailable fits: ",
         paste(sub("\\.draws\\.rds$", "",
                   list.files(summaries_dir, pattern = "\\.draws\\.rds$")),
               collapse = ", "))
  }
  draws <- readRDS(draws_path)

  # stan_data: try a fit-name-specific path, fall back to long_subset_data
  stan_path <- stan_data_paths[[fit_name]]
  if (is.null(stan_path) || !file.exists(stan_path)) {
    stan_path <- file.path(PATHS$fits_dir, "long_subset_data.rds")
  }
  sd <- readRDS(stan_path)$stan_data

  # Some fits don't include rho_xi_zeta or sigma_zeta (e.g., baseline w/o
  # slopes). Default those to 0 with a warning.
  has_rho   <- "rho_xi_zeta" %in% names(draws)
  has_zeta  <- "sigma_zeta"  %in% names(draws)
  if (!has_zeta) {
    warning(fit_name, ": no sigma_zeta in posterior; treating as 0 ",
            "(no per-child slope variance).")
  }

  # Coerce posterior columns to plain numeric vectors — draws_df columns
  # otherwise carry metadata that breaks downstream tibble construction
  # (length-recycling surprises).
  as_num <- function(x) as.numeric(unlist(x))

  delta_v       <- as_num(draws$delta)
  sigma_alpha_v <- as_num(draws$sigma_alpha)
  sigma_xi_v    <- if ("sigma_xi" %in% names(draws)) as_num(draws$sigma_xi)
                   else sqrt(sigma_alpha_v^2 + sd$sigma_r^2)
  sigma_zeta_v  <- if (has_zeta) as_num(draws$sigma_zeta) else rep(0, length(delta_v))
  rho_v         <- if (has_rho)  as_num(draws$rho_xi_zeta) else rep(0, length(delta_v))
  s_v           <- if ("s" %in% names(draws)) as_num(draws$s) else rep(0.5, length(delta_v))

  list(
    fit_name        = fit_name,
    n_draws         = length(delta_v),
    delta           = delta_v,
    sigma_xi        = sigma_xi_v,
    sigma_zeta      = sigma_zeta_v,
    sigma_alpha     = sigma_alpha_v,
    rho_xi_zeta     = rho_v,
    s               = s_v,
    mu_r            = sd$mu_r,
    sigma_r         = sd$sigma_r,
    a0              = sd$a0,
    log_H           = sd$log_H,
    age_obs_range   = range(sd$admin_age)
  )
}

# ---- Simulate per-child trajectories from posterior ----------------
##
## For each posterior draw, sample n_kids_per_draw kids' (xi, zeta) from
## the population MVN and trace a trajectory across age. Returns long-form
## tibble with one row per (draw, kid, age).
##
## "Sample-from-posterior-of-pop" rather than per-child posteriors keeps
## the script general — it doesn't need to load the full per-child draws.
simulate_kid_trajectories <- function(
  P,
  n_kids_per_draw = 30,
  n_draws_use     = 60,
  age_grid        = seq(12, 30, by = 0.5),
  seed            = 42
) {
  set.seed(seed)
  draw_idx <- sample.int(P$n_draws, size = min(n_draws_use, P$n_draws))

  rows <- vector("list", length(draw_idx))
  for (k in seq_along(draw_idx)) {
    d <- draw_idx[k]
    sx <- P$sigma_xi[d];  sz <- P$sigma_zeta[d];  rho <- P$rho_xi_zeta[d]
    delta_d <- P$delta[d]; s_d <- P$s[d]

    # MVN draws of (xi, zeta) per kid
    Sigma <- matrix(c(sx^2,            rho * sx * sz,
                      rho * sx * sz,   sz^2 + 1e-12), 2, 2)
    L <- chol(Sigma)
    z <- matrix(rnorm(2 * n_kids_per_draw), 2, n_kids_per_draw)
    effs <- t(L) %*% z
    xi   <- P$mu_r + effs[1, ]
    zeta <- effs[2, ]

    # For each kid, compute log-D(t) and theta(t) on the age grid.
    # log r_i is part of xi; we use it to put per-kid trajectories at the
    # right point on the x-axis. (xi is log r_i + log alpha_i; the log r_i
    # piece lives in xi - log alpha_i, and we approximate with the marginal
    # mean of log r_i conditioned on xi: same conditional mean as in the
    # explainer, used here only for x-axis placement, not interpretation.)
    sigma_alpha_d <- P$sigma_alpha[d]
    cond_log_r <- P$mu_r + (P$sigma_r^2 / sx^2) * (xi - P$mu_r)

    age_eff <- pmax(age_grid - s_d, 0.01)
    log_age_rel <- log(age_eff / P$a0)

    # log10(D_i(t)) = (log r_i + log H + log(t - s)) / log(10)
    log_D <- outer(rep(1, length(age_grid)),
                   cond_log_r + P$log_H) +
             outer(log(age_eff), rep(1, n_kids_per_draw))
    log10_D <- log_D / log(10)

    theta <- outer(rep(1, length(age_grid)), xi) +
             outer(log_age_rel, 1 + delta_d + zeta)
    # Anchor population intercept at 0: subtract mu_r so y-axis is
    # "ability beyond population log-rate baseline."
    theta_centered <- theta - P$mu_r

    # Compute all per-row vectors *outside* the tibble call — tibble lets
    # later columns reference earlier ones, which would shadow the local
    # `zeta` (length 30) with the just-assigned `zeta` column (length 1110).
    n_age <- length(age_grid)
    slope_vec <- 1 + delta_d + zeta
    rows[[k]] <- tibble::tibble(
      draw    = d,
      kid     = rep(seq_len(n_kids_per_draw), each = n_age),
      age     = rep(age_grid, n_kids_per_draw),
      xi_kid  = rep(xi,        each = n_age),
      zeta_kid = rep(zeta,     each = n_age),
      slope   = rep(slope_vec, each = n_age),
      log10_D = as.vector(log10_D),
      theta   = as.vector(theta_centered)
    )
  }

  bind_rows(rows)
}

# ---- LLM reference lines --------------------------------------------
##
## Chinchilla et al. lines as straight segments on log-log axes.
## We pin the y-axis offset for visual placement only — the y-axis is
## "log-likelihood units (offset)" and the structurally meaningful quantity
## is the SLOPE.
build_llm_lines <- function(
  llm_d_range_log10 = LLM_D_RANGE_LOG10,
  y_anchor_at_logD  = 11,            # where lines pass through y = 0
  refs              = LLM_REF
) {
  lines <- lapply(names(refs), function(nm) {
    r <- refs[[nm]]
    logD <- seq(llm_d_range_log10[1], llm_d_range_log10[2], length.out = 50)
    # y = beta * (logD - anchor)  →  passes through (anchor, 0)
    y    <- r$beta * (logD - y_anchor_at_logD) * log(10)  # convert log10 → ln
    tibble::tibble(name = nm, label = r$label,
                   colour = r$colour, beta = r$beta,
                   log10_D = logD, theta = y)
  })
  bind_rows(lines)
}

# ---- Panel A: scaling diagram --------------------------------------
panel_scaling_diagram <- function(kid_traj, llm_lines, P) {

  # Population mean line for kids: slope = 1 + delta over the kid D range
  age_grid <- seq(12, 30, by = 0.5)
  age_eff <- pmax(age_grid - median(P$s), 0.01)
  log_age_rel <- log(age_eff / P$a0)
  pop_log10_D <- (P$mu_r + P$log_H + log(age_eff)) / log(10)
  pop_theta   <- median(P$delta + 1) * log_age_rel
  pop_line <- tibble::tibble(log10_D = pop_log10_D, theta = pop_theta)

  # Slope labels — placed to the right of each line's right endpoint, with
  # explicit nudge (in data units) so they clear the line cluster.
  # All three LLM lines pass close to y=0 at right end (slopes are tiny),
  # so we stagger their label y-positions explicitly.
  llm_labels <- llm_lines |>
    group_by(name, label, colour, beta) |>
    slice_max(log10_D, n = 1) |>
    ungroup() |>
    mutate(text = sprintf("β = %.3f", beta)) |>
    arrange(desc(beta)) |>
    mutate(y_offset = c(2.0, 0.5, -1.5)[seq_len(dplyr::n())],
           theta_lbl = y_offset)

  pop_label <- tibble::tibble(
    log10_D = max(pop_line$log10_D),
    theta   = max(pop_line$theta),
    text    = sprintf("1 + δ ≈ %.1f", median(P$delta) + 1)
  )

  ggplot() +
    # Vertical guide separating the two D regimes (purely visual)
    geom_vline(xintercept = 8.5, colour = "grey85",
               linetype = "dotted", linewidth = 0.4) +
    # Kid trajectories — light grey fan
    geom_line(data = kid_traj,
              aes(log10_D, theta, group = interaction(draw, kid)),
              alpha = 0.05, colour = "grey25") +
    # Population mean kid line — bold red
    geom_line(data = pop_line,
              aes(log10_D, theta), colour = "firebrick", linewidth = 1.3) +
    # LLM lines
    geom_line(data = llm_lines,
              aes(log10_D, theta, group = name, colour = label),
              linewidth = 0.95) +
    # Slope labels — pushed clear of lines
    geom_text(data = pop_label,
              aes(log10_D, theta, label = text),
              nudge_x = 0.25, hjust = 0, vjust = 0.5,
              size = 3.5, fontface = "bold", colour = "firebrick") +
    geom_segment(data = llm_labels,
                 aes(x = log10_D, xend = log10_D + 0.15,
                     y = theta, yend = theta_lbl, colour = label),
                 linewidth = 0.4, show.legend = FALSE) +
    geom_text(data = llm_labels,
              aes(log10_D + 0.18, theta_lbl,
                  label = text, colour = label),
              hjust = 0, vjust = 0.5,
              size = 3.2, fontface = "bold", show.legend = FALSE) +
    # Region labels at top
    annotate("text", x = 7,    y = 14, label = "Children\n(SM2 fit)",
             size = 3.5, fontface = "bold", colour = "firebrick") +
    annotate("text", x = 10.7, y = 14, label = "LLMs\n(scaling laws)",
             size = 3.5, fontface = "bold", colour = "#1f78b4") +
    scale_colour_manual(values = setNames(sapply(LLM_REF, `[[`, "colour"),
                                          sapply(LLM_REF, `[[`, "label"))) +
    coord_cartesian(xlim = c(5.8, 14), ylim = c(-15, 15), clip = "off") +
    labs(
      x = expression("Cumulative tokens "*log[10]*"(D)"),
      y = "Ability / evidence  (log-likelihood units, offset)",
      colour = NULL,
      title = "Scaling within and across children vs. LLM scaling laws"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position    = "bottom",
          legend.box.margin  = margin(t = -5),
          plot.margin        = margin(8, 14, 8, 8),
          plot.title         = element_text(face = "bold"))
}

# ---- Panel B: distribution of per-instance scaling exponents -------
panel_exponent_dist <- function(kid_traj, P) {

  # One slope per (draw, kid)
  slopes <- kid_traj |>
    distinct(draw, kid, slope)

  llm_pts <- tibble::tibble(
    label  = sapply(LLM_REF, `[[`, "label"),
    beta   = sapply(LLM_REF, `[[`, "beta"),
    colour = sapply(LLM_REF, `[[`, "colour")
  )

  # Range to show: include both LLM betas and kid slopes, with headroom
  x_max <- max(slopes$slope, na.rm = TRUE) + 1
  x_min <- -0.5

  # Density for kids; staircase placement of LLM labels to avoid overlap
  llm_pts <- llm_pts |>
    arrange(desc(beta)) |>
    mutate(label_y_frac = c(0.85, 0.65, 0.45)[seq_len(dplyr::n())])

  d <- density(slopes$slope, n = 256)
  ymax <- max(d$y) * 1.15

  ggplot() +
    geom_density(data = slopes, aes(slope),
                 fill = "grey60", colour = "grey30", alpha = 0.7) +
    geom_vline(data = llm_pts,
               aes(xintercept = beta, colour = label),
               linewidth = 0.9, linetype = "dashed") +
    geom_text(data = llm_pts,
              aes(x = beta, y = ymax * label_y_frac,
                  label = sprintf("%s  (β = %.3f)", label, beta),
                  colour = label),
              hjust = -0.05, vjust = 0.5, size = 3.0,
              show.legend = FALSE) +
    geom_segment(data = llm_pts,
                 aes(x = beta, xend = beta + 0.3,
                     y = ymax * label_y_frac, yend = ymax * label_y_frac,
                     colour = label),
                 linewidth = 0.4, show.legend = FALSE) +
    scale_colour_manual(values = setNames(llm_pts$colour, llm_pts$label)) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, ymax)) +
    labs(
      x = expression("Per-instance scaling exponent"),
      y = "Density (kids)",
      title = "Per-instance scaling exponents"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold"))
}

# ---- Compose, render & save ----------------------------------------
make_scaling_disanalogy_plot <- function(
  fit_name        = "long_slopes",
  output_basename = "D1_scaling_disanalogy",
  n_kids_per_draw = 30,
  n_draws_use     = 60,
  age_grid        = seq(12, 30, by = 0.5),
  width           = 13,
  height          = 6.5,
  dpi             = 150
) {
  P <- get_kid_scaling_params(fit_name)
  kid_traj <- simulate_kid_trajectories(P,
                                        n_kids_per_draw = n_kids_per_draw,
                                        n_draws_use     = n_draws_use,
                                        age_grid        = age_grid)
  llm_lines <- build_llm_lines()

  pA <- panel_scaling_diagram(kid_traj, llm_lines, P)
  pB <- panel_exponent_dist(kid_traj, P)

  composite <- pA + pB + plot_layout(widths = c(2, 1)) +
    plot_annotation(
      subtitle = "Slopes are the structurally meaningful comparison; y-axis offset is anchored for visual clarity. Kid D range comes from per-kid r_i.",
      caption = sprintf(
        "SM2 fit '%s': median 1+δ = %.2f [%.2f, %.2f]; σ_ζ = %.2f [%.2f, %.2f]; mu_r=%.2f, sigma_r=%.2f, a₀=%.0f.  Chinchilla and Kaplan exponents from Hoffmann et al. 2022 / Kaplan et al. 2020.",
        fit_name,
        median(P$delta) + 1,
        quantile(P$delta + 1, 0.025),
        quantile(P$delta + 1, 0.975),
        median(P$sigma_zeta),
        quantile(P$sigma_zeta, 0.025),
        quantile(P$sigma_zeta, 0.975),
        P$mu_r, P$sigma_r, P$a0
      ),
      theme = theme(
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        plot.caption  = element_text(size = 8, colour = "grey40", hjust = 0)
      )
    )

  out_png <- file.path(OUT_DIR, paste0(output_basename, ".png"))
  ggsave(out_png, composite, width = width, height = height, dpi = dpi)
  out_rds <- file.path(OUT_DIR, paste0(output_basename, "_data.rds"))
  saveRDS(list(kid_traj = kid_traj, llm_lines = llm_lines, params = P),
          out_rds)
  cat("Saved:", out_png, "\n      ", out_rds, "\n")

  invisible(list(plot = composite, kid_traj = kid_traj, llm_lines = llm_lines, P = P))
}

# ---- CLI ------------------------------------------------------------
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  fit_name <- if (length(args) >= 1) args[1] else "long_slopes"
  output_basename <- if (length(args) >= 2) args[2] else "D1_scaling_disanalogy"
  make_scaling_disanalogy_plot(fit_name = fit_name,
                               output_basename = output_basename)
}
