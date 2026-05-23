## Helper: build cross-sectional empirical reference for the quantile-fan
## plots. Fixes two artifacts of the longitudinal-as-cross-sectional
## approach used previously:
##   1. One admin per child (cross-sectional). Longitudinal data has
##      multiple admins per kid; treating them as separate dots biases
##      the empirical fan toward "stayer" kids.
##   2. Lower smoothing on the empirical quantile-regression fit
##      (lambda = 50 vs. previous 1000). Less aggressive flattening of
##      real curvature in the empirical quantiles.
##
## Used by: model/scripts/quantile_demo_*.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(quantregGrowth)
})

# Build cross-sectional empirical (one admin per child).
# Strategy: for each child, sample one admin at random from those with
# valid `produces` data within the age range. Restricts to bundle items.
# Build cross-sectional empirical from the broader wordbank corpus
# (long_items.rds), restricted to the bundle's item set and the
# specified language / form. Properly deduplicates by (child_id,
# age, item) -- long_items.rds repeats the same observation 2-3x for
# some kids (likely a wordbankr aggregation artifact). All repeated
# rows have identical `produces` values, so distinct() is safe.
#
# Use this when you want a much larger N empirical reference than
# the bundle alone provides -- the bundle is a stratified subset
# of ~500 kids; wordbank has ~thousands more matching the filters.
build_xsec_empirical_wordbank <- function(item_definitions,
                                           language_filter,
                                           form_filter = "WS",
                                           age_range = c(16, 30),
                                           long_items_path,
                                           cache_path = NULL,
                                           seed = 20260523) {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cat("Loading cached wordbank x-sec empirical:", cache_path, "\n")
    return(readRDS(cache_path))
  }
  if (!file.exists(long_items_path))
    stop("long_items.rds not found at ", long_items_path)
  d_long <- readRDS(long_items_path)

  admins <- d_long |>
    filter(language == language_filter, form == form_filter,
           item %in% item_definitions,
           age >= age_range[1], age <= age_range[2],
           !is.na(produces)) |>
    distinct(child_id, age, item, .keep_all = TRUE) |>
    group_by(child_id, age) |>
    summarise(vocab = sum(produces, na.rm = TRUE),
              n_items = n(), .groups = "drop")

  set.seed(seed)
  emp <- admins |>
    group_by(child_id) |>
    slice_sample(n = 1) |>
    ungroup()

  if (!is.null(cache_path)) saveRDS(emp, cache_path)
  emp
}

build_xsec_empirical <- function(bundle_df,
                                  age_range = c(16, 30),
                                  cache_path = NULL,
                                  seed = 20260523) {
  # Build cross-sectional empirical (one admin per child) from the
  # bundle's $df (which carries proper admin_key / aa identifying
  # individual admins). long_items.rds was previously used, but it
  # has no admin_id column -- multiple admins at the same age for the
  # same kid got collapsed into a single row, summing produces across
  # admins (max n_items hit 2040 = 3 admins x 680 items, inflating
  # vocab counts beyond J and breaking BEINF's [0,1] constraint).
  #
  # Arguments:
  #   bundle_df : the $df slot of a bundle (long format, one row per
  #               (admin, item)). Required cols: admin_key, child_id,
  #               age, produces.
  #   age_range : age window to keep, in months.
  #   cache_path: optional .rds path to memoize result.
  #   seed      : RNG seed for the random sub-sample per child.
  #
  # Returns a data.frame with columns: child_id, admin_key, age, vocab,
  # n_items.
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cat("Loading cached x-sec empirical:", cache_path, "\n")
    return(readRDS(cache_path))
  }

  # Different bundle pipelines use different child-ID columns:
  #   long_subset_data.rds (EN, NO):           child_id
  #   babyview_subset_data.rds, seedlings_*:   subject_id
  #   stanford_linked_subset_data.rds:         lab_subject_id
  # Detect which is present and aliased it to a stable name.
  id_col <- intersect(c("child_id", "subject_id", "lab_subject_id"),
                      names(bundle_df))
  if (length(id_col) == 0)
    stop("bundle_df has no recognized child-ID column (looked for: ",
         "child_id, subject_id, lab_subject_id)")
  id_col <- id_col[1]
  d <- bundle_df
  d$child_id <- d[[id_col]]

  # NOTE: bundle$df can contain duplicate (admin_key, item) rows --
  # ~3% of rows are duplicates in the EN bundle, where the same (kid,
  # age, item, produces) tuple appears 2-3 times. The model treats
  # these as independent observations; for the empirical reference,
  # deduplicate first so vocab counts max out at J (not 2J or 3J).
  admins <- d |>
    filter(age >= age_range[1], age <= age_range[2],
           !is.na(produces)) |>
    distinct(admin_key, item, .keep_all = TRUE) |>
    group_by(child_id, age, admin_key) |>
    summarise(vocab = sum(produces, na.rm = TRUE),
              n_items = n(), .groups = "drop")

  # One admin per child, randomly sampled. (Random, not earliest /
  # latest, so the age distribution isn't biased toward either end.)
  set.seed(seed)
  emp <- admins |>
    group_by(child_id) |>
    slice_sample(n = 1) |>
    ungroup()

  if (!is.null(cache_path)) saveRDS(emp, cache_path)
  emp
}

# Smooth monotone quantile-regression fan on the x-sec empirical.
# Lower lambda (50) vs. the historical 1000 -- preserves more real
# curvature without overfitting noise. Returns long-format predictions
# at the requested `age_grid` for `taus`.
fit_xsec_quantile_fan <- function(emp_df,
                                   taus = c(0.10, 0.25, 0.50, 0.75, 0.90),
                                   age_grid,
                                   J_total,
                                   lambda = 10000) {
  # GAMLSS quantile estimation, per the MB-CDI manual approach (Frank
  # et al., demo-vocab/gamlss_demo.R): fit prop_produced ~ pbm(age,
  # lambda) with a beta-family response, then read out centiles via
  # centiles.pred(). Handles the bounded [0, 1] support of vocab
  # proportions properly (unlike quantreg::rq on raw counts) and gives
  # smooth, monotone-by-construction quantile curves.
  #
  # Uses BEINF (beta inflated at 0 AND 1) so zero-vocab and ceiling
  # admins don't have to be filtered out. Returns the empirical
  # quantile fan in long format with columns: tau, age, vocab (where
  # vocab is on the raw-count scale, prop * J_total).
  if (!requireNamespace("gamlss", quietly = TRUE))
    stop("Install gamlss: install.packages('gamlss')")
  if (!requireNamespace("gamlss.dist", quietly = TRUE))
    stop("Install gamlss.dist: install.packages('gamlss.dist')")
  # gamlss formulas reference pbm/pb as bare symbols, looked up in the
  # search path -- gamlss::pbm() qualification doesn't survive formula
  # evaluation. So attach gamlss to the search path explicitly.
  suppressPackageStartupMessages(suppressWarnings(library(gamlss)))

  # gamlss is unhappy with tibble inputs (model.frame construction
  # ends up passing non-numeric subsets into pbm). Force a plain
  # data.frame and explicit numerics.
  d <- as.data.frame(emp_df)
  d$prop_produced <- as.numeric(d$vocab) / J_total
  d$prop_produced <- pmin(pmax(d$prop_produced, 0), 1)
  d$age <- as.numeric(d$age)
  d <- d[, c("age", "prop_produced"), drop = FALSE]

  # gamlss + centiles.pred use non-standard evaluation -- centiles.pred
  # looks up the data by NAME (model$call$data) in the calling env,
  # not the function env. So we have to expose the data in a place
  # centiles.pred can find. The least-invasive hack: write `d` into
  # the global env under a unique name, call gamlss with data = that
  # name, then clean up afterwards.
  gamlss_data_name <- paste0(".gamlss_xsec_", as.integer(Sys.time()))
  assign(gamlss_data_name, d, envir = globalenv())
  on.exit(rm(list = gamlss_data_name, envir = globalenv()), add = TRUE)
  data_sym <- as.name(gamlss_data_name)

  # mu = expected proportion; sigma = dispersion. pbm gives a monotone
  # P-spline on mu so the median (and other centiles via the link)
  # are non-decreasing in age. pb on sigma lets variance change with age.
  fml_mu    <- eval(substitute(prop_produced ~ pbm(age, lambda = LAM),
                                list(LAM = lambda)))
  fml_sigma <- ~ pb(age)
  # Build the gamlss call manually so model$call$data is the symbol,
  # not the local variable -- this is what centiles.pred wants.
  gamlss_call <- bquote(
    gamlss::gamlss(formula = .(fml_mu),
                   sigma.formula = .(fml_sigma),
                   family = gamlss.dist::BEINF(),
                   data = .(data_sym),
                   trace = FALSE)
  )
  m <- tryCatch(
    suppressWarnings(suppressMessages(eval(gamlss_call))),
    error = function(e) { message("gamlss error: ", conditionMessage(e)); NULL }
  )
  if (is.null(m)) return(NULL)

  # centiles.pred() returns a data frame with columns x, then one
  # column per centile. cent argument expects 0-100 percentile values.
  cents <- gamlss::centiles.pred(m,
                                  cent = taus * 100,
                                  xname = "age",
                                  xvalues = age_grid)
  # First column is x (age), remaining are the centile columns named
  # by the requested numbers.
  out <- as.data.frame(cents)
  age_col <- out[[1]]
  pred_cols <- out[, -1, drop = FALSE]
  names(pred_cols) <- as.character(taus)
  # Convert proportions back to vocab counts on the bundle's J_total scale.
  pred_cols <- pred_cols * J_total
  out <- data.frame(age = age_col, pred_cols, check.names = FALSE)
  tidyr::pivot_longer(out, cols = -age, names_to = "tau", values_to = "vocab")
}
