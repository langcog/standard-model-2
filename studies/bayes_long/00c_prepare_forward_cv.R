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
## MIN_ELIG (2nd arg, default 3) -- IMPORTANT. At MIN_ELIG=3, holding out the last of
## three administrations leaves TWO in training: precisely the regime this project
## established cannot identify a per-child slope (one gap, zero residual df), and the
## reason the main analysis requires 3+. In the shallow datasets nearly every child sits
## at exactly three administrations (Thal 637/639, Marchman 149/162, Japanese 95/96), so
## MIN_ELIG=3 tests M3 in a configuration the paper explicitly excludes. Use MIN_ELIG=4
## so training retains >=3; only Norwegian (602 children) and Smith (77) have the depth.
##
## At MIN_ELIG>3 ineligible children are dropped from TRAIN as well, not just from the
## test set, so every child in the model has an identifiable slope. Child indices are
## therefore renumbered, and the test set carries the NEW indices.
##
## Usage:  Rscript studies/bayes_long/00c_prepare_forward_cv.R [suffix] [min_elig] [slug...]
##   e.g.  Rscript studies/bayes_long/00c_prepare_forward_cv.R _a3 4 norwegian
## Output: bundle_<slug>_fcv.rds (min_elig 3) or bundle_<slug>_fcv<N>.rds (min_elig N>3)

suppressPackageStartupMessages({library(dplyr)})
.args <- commandArgs(trailingOnly = TRUE)
SFX      <- if (length(.args) >= 1) .args[1] else "_a3"
MIN_ELIG <- if (length(.args) >= 2) as.integer(.args[2]) else 3L
BL    <- file.path("fits", "bayes_long")
SLUGS <- if (length(.args) >= 3) .args[-(1:2)] else
         c("thal", "smith", "marchman", "norwegian", "japanese")
TAG   <- if (MIN_ELIG > 3L) sprintf("_fcv%d", MIN_ELIG) else "_fcv"

one <- function(slug) {
  f <- file.path(BL, sprintf("bundle_%s%s.rds", slug, SFX))
  if (!file.exists(f)) { cat("skip", slug, "(no bundle)\n"); return(invisible(NULL)) }
  b <- readRDS(f); sd0 <- b$stan_data

  ## ---- choose the held-out administration per child: the LAST by age -------
  adm <- tibble(a = seq_len(sd0$A), child = sd0$admin_to_child, age = sd0$admin_age)
  n_by_child <- adm |> count(child, name = "n_adm")
  eligible <- sort(n_by_child$child[n_by_child$n_adm >= MIN_ELIG])
  if (!length(eligible)) { cat("skip", slug, "(no children with >=", MIN_ELIG, "admins)\n"); return(invisible(NULL)) }

  ## drop ineligible children entirely, then renumber children 1..I_keep
  keep_child <- seq_len(sd0$I) %in% eligible
  cmap <- integer(sd0$I); cmap[eligible] <- seq_along(eligible)
  adm_keep <- keep_child[sd0$admin_to_child]              # admins belonging to kept children
  test_a <- adm |> filter(child %in% eligible) |>
    group_by(child) |> slice_max(age, n = 1, with_ties = FALSE) |> ungroup() |> pull(a)
  is_test <- seq_len(sd0$A) %in% test_a
  ## an admin is TRAIN if its child is kept and it is not the held-out one
  is_train <- adm_keep & !is_test
  sd0$admin_to_child <- cmap[sd0$admin_to_child]          # renumbered child index per admin
  sd0$I <- length(eligible)

  ## ---- reindex administrations; child indices are the RENUMBERED ones above ----
  ## Both train and test carry the same (renumbered) child index, so the held-out
  ## administration scores against exactly the xi_i / kappa_i the training fit estimated.
  ## Item indices are untouched, so delta_j transfers directly.
  tr_a    <- which(is_train); te_a <- which(is_test)
  new_of  <- integer(sd0$A); new_of[tr_a] <- seq_along(tr_a)
  new_te  <- integer(sd0$A); new_te[te_a] <- seq_along(te_a)
  tr_obs  <- is_train[sd0$aa]
  te_obs  <- is_test[sd0$aa]

  train <- list(
    N = sum(tr_obs), A = length(tr_a), I = sd0$I, J = sd0$J,
    aa = new_of[sd0$aa[tr_obs]], jj = sd0$jj[tr_obs], y = sd0$y[tr_obs],
    admin_to_child = sd0$admin_to_child[tr_a], admin_age = sd0$admin_age[tr_a],
    log_H = sd0$log_H, a0 = sd0$a0)
  test <- list(
    N = sum(te_obs), A = length(te_a),
    aa = new_te[sd0$aa[te_obs]], jj = sd0$jj[te_obs], y = sd0$y[te_obs],
    admin_to_child = sd0$admin_to_child[te_a], admin_age = sd0$admin_age[te_a])

  ## every kept child must retain >= MIN_ELIG-1 training administrations
  tr_per_child <- tabulate(train$admin_to_child, nbins = sd0$I)
  stopifnot(max(train$aa) == train$A, max(test$aa) == test$A,
            all(train$admin_to_child >= 1), all(train$admin_to_child <= sd0$I),
            all(tr_per_child >= MIN_ELIG - 1L),
            length(te_a) == sd0$I,                       # exactly one held-out admin per child
            !any(is_train & is_test))

  meta <- list(slug = slug, sfx = SFX, min_elig = MIN_ELIG,
               n_kids = sd0$I, n_kids_tested = length(eligible),
               A_train = train$A, A_test = test$A, N_train = train$N, N_test = test$N,
               train_adm_per_child = mean(tr_per_child),
               test_age_range = range(sd0$admin_age[te_a]),
               train_age_range = range(sd0$admin_age[tr_a]))
  out <- file.path(BL, sprintf("bundle_%s%s.rds", slug, TAG))
  saveRDS(list(stan_data = train, test = test, meta = meta), out)
  cat(sprintf("%-10s kids=%4d | train A=%5d (%.2f/kid) N=%8d | test A=%4d N=%7d | train %.0f-%.0f, test %.0f-%.0f\n",
              slug, sd0$I, train$A, mean(tr_per_child), train$N, test$A, test$N,
              meta$train_age_range[1], meta$train_age_range[2],
              meta$test_age_range[1], meta$test_age_range[2]))
}

cat(sprintf("Forward-CV bundles from '%s' bundles: hold out each child's LAST administration,\n", SFX))
cat(sprintf("  eligibility >= %d administrations (training retains >= %d) -> tag '%s'\n",
            MIN_ELIG, MIN_ELIG - 1L, TAG))
invisible(lapply(SLUGS, one))
