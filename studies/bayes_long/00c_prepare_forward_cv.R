## 00c_prepare_forward_cv.R -- build TRAIN/TEST bundles for forward (prospective) CV.
##
## WHY. The ladder's LOO holds out ONE item response at a time. A child contributes
## several hundred responses per administration, so the held-out response stays highly
## predictable from that child's remaining responses, that word's responses in other
## children, and often the same administration. Item-level LOO therefore measures
## local interpolation, and it structurally cannot penalise an over-flexible per-child
## slope: a kappa_i that merely interpolates a child's own sittings still scores well.
##
## The honest test of an ACCELERATION model is prospective: fit on a child's earlier
## administrations and predict their LATER one. That is what this builds.
##
## Split: for every child, the final administration (by age) goes to TEST; all earlier
## ones stay in TRAIN. Children need >=3 administrations so TRAIN retains >=2 -- the
## same identifiability requirement as the main analysis (one gap = one slope with zero
## residual df). Items are held fixed across the split so item difficulties transfer.
##
## Cost note: this needs only ONE refit per (dataset, model) -- every child's last
## administration is held out simultaneously -- not K refits. PSIS cannot approximate
## it because dropping a whole administration is far too large a perturbation for
## importance weights to survive.
##
## Usage:  Rscript studies/bayes_long/00c_prepare_forward_cv.R [suffix]     # default _a3
## Output: fits/bayes_long/bundle_<slug>_fcv.rds  = list(stan_data=TRAIN, test=..., meta=...)

suppressPackageStartupMessages({library(dplyr)})
SFX   <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else "_a3"
BL    <- file.path("fits", "bayes_long")
SLUGS <- c("thal", "smith", "marchman", "norwegian", "japanese")

one <- function(slug) {
  f <- file.path(BL, sprintf("bundle_%s%s.rds", slug, SFX))
  if (!file.exists(f)) { cat("skip", slug, "(no bundle)\n"); return(invisible(NULL)) }
  b <- readRDS(f); sd0 <- b$stan_data

  ## ---- choose the held-out administration per child: the LAST by age -------
  adm <- tibble(a = seq_len(sd0$A), child = sd0$admin_to_child, age = sd0$admin_age)
  n_by_child <- adm |> count(child, name = "n_adm")
  ## only children with >=3 admins can give up one and still identify a slope
  eligible <- n_by_child$child[n_by_child$n_adm >= 3]
  test_a <- adm |> filter(child %in% eligible) |>
    group_by(child) |> slice_max(age, n = 1, with_ties = FALSE) |> ungroup() |> pull(a)
  is_test <- seq_len(sd0$A) %in% test_a

  ## ---- reindex TRAIN administrations (children/items keep their global index) ----
  ## Child and item indices are deliberately NOT renumbered: the test set must score
  ## against the same xi_i / kappa_i / delta_j the training fit estimates.
  keep_a  <- which(!is_test)
  new_of  <- integer(sd0$A); new_of[keep_a] <- seq_along(keep_a)
  tr_obs  <- !is_test[sd0$aa]
  te_obs  <-  is_test[sd0$aa]

  train <- list(
    N = sum(tr_obs), A = length(keep_a), I = sd0$I, J = sd0$J,
    aa = new_of[sd0$aa[tr_obs]], jj = sd0$jj[tr_obs], y = sd0$y[tr_obs],
    admin_to_child = sd0$admin_to_child[keep_a], admin_age = sd0$admin_age[keep_a],
    log_H = sd0$log_H, a0 = sd0$a0)
  ## TEST keeps ORIGINAL child indexing so scoring lines up with the fitted per-child
  ## parameters; aa here indexes rows of the test administration table.
  te_a    <- which(is_test)
  new_te  <- integer(sd0$A); new_te[te_a] <- seq_along(te_a)
  test <- list(
    N = sum(te_obs), A = length(te_a),
    aa = new_te[sd0$aa[te_obs]], jj = sd0$jj[te_obs], y = sd0$y[te_obs],
    admin_to_child = sd0$admin_to_child[te_a], admin_age = sd0$admin_age[te_a])

  stopifnot(train$N + test$N == sd0$N,
            max(train$aa) == train$A, max(test$aa) == test$A,
            all(train$admin_to_child >= 1), all(train$admin_to_child <= sd0$I))

  meta <- list(slug = slug, sfx = SFX, n_kids = sd0$I, n_kids_tested = length(eligible),
               A_train = train$A, A_test = test$A, N_train = train$N, N_test = test$N,
               test_age_range = range(sd0$admin_age[te_a]),
               train_age_range = range(sd0$admin_age[keep_a]))
  out <- file.path(BL, sprintf("bundle_%s_fcv.rds", slug))
  saveRDS(list(stan_data = train, test = test, meta = meta), out)
  cat(sprintf("%-10s kids=%4d tested=%4d | train A=%5d N=%8d | test A=%4d N=%7d | test ages %.0f-%.0f\n",
              slug, sd0$I, length(eligible), train$A, train$N, test$A, test$N,
              meta$test_age_range[1], meta$test_age_range[2]))
}

cat(sprintf("Forward-CV bundles from '%s' bundles (hold out each child's LAST administration)\n", SFX))
invisible(lapply(SLUGS, one))
