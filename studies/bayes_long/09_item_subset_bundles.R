## 09_item_subset_bundles.R -- does kappa track the learner, or the ruler?
##
## THE OBJECTION. kappa ~ 11 might be nothing but the spread of CDI item difficulties
## divided by the width of the observation window. If so it is a property of the
## checklist, not of the child, and the 11-versus-1 comparison against language models is
## not a comparison of learners at all.
##
## THE TEST. Refit the same model on the same children and the same administrations, with
## only the ITEM SET narrowed. If kappa is a rescaling of the instrument, halving the
## retained difficulty spread should roughly halve kappa. If kappa is anchored to the
## children, it should barely move.
##
## WHY A SIZE-MATCHED RANDOM CONTROL IS ESSENTIAL. Narrowing by difficulty does two things
## at once: it shrinks the retained spread (the effect of interest) and it throws away half
## the data (a pure precision loss that can move kappa on its own). A random subset of the
## SAME SIZE loses the same information while preserving the spread, so the contrast
## between them isolates the part that matters. The original proposal had one random half;
## every narrowed condition needs its own size-matched control, hence rand25 as well.
##
## Conditions (per dataset):
##   all      every item (the existing fit; no bundle written)
##   mid50    middle 50% by fitted difficulty -- narrows spread AND truncates the extremes
##   mid25    middle 25% -- the same, harder
##   easy50   easiest half        } same size as mid50/rand50, roughly preserved spread,
##   hard50   hardest half        } shifted location. The cleanest discriminator: if kappa
##                                  is real both halves recover it, each informative for a
##                                  different age band; if kappa is the ruler, each gives
##                                  about half.
##   rand50   random half   -- size-matched control for mid50/easy50/hard50
##   rand25   random 25%    -- size-matched control for mid25
##
## Difficulty ordering comes from the fitted delta_j of the baseline _a3 M3 fit, joined to
## the bundle by the item index jj that the psi export carries.
##
## Children and administrations are NEVER dropped. An older child answering only hard items
## still contributes an all-zero pattern, which is informative under the model; removing
## them would change the sample and confound the very comparison being made. Only
## administrations left with literally zero observations are dropped, since a zero-length
## group breaks the Stan indexing.
##
## Usage:  Rscript studies/bayes_long/09_item_subset_bundles.R <slug_a3> [...]
## Output: fits/bayes_long/bundle_<slug>_<cond>.rds  + a console manifest

suppressPackageStartupMessages({library(dplyr)})
BL   <- file.path("fits", "bayes_long")
SUMM <- file.path(BL, "summaries")
SEED <- 20260811L

slugs <- commandArgs(trailingOnly = TRUE)
if (!length(slugs)) slugs <- c("thal_a3", "marchman_a3")

subset_bundle <- function(b, keep_jj, cond) {
  sd0 <- b$stan_data
  keep_jj <- sort(unique(keep_jj))
  ## reindex items to 1..J_new; map[old] -> new
  jmap <- integer(sd0$J); jmap[keep_jj] <- seq_along(keep_jj)
  ix <- which(sd0$jj %in% keep_jj)

  aa <- sd0$aa[ix]; jj <- jmap[sd0$jj[ix]]; y <- sd0$y[ix]
  ## drop administrations with no surviving observations, then reindex
  live_a <- sort(unique(aa))
  amap <- integer(sd0$A); amap[live_a] <- seq_along(live_a)
  aa <- amap[aa]

  new <- list(N = length(ix), A = length(live_a), I = sd0$I, J = length(keep_jj),
              aa = aa, jj = jj, y = y,
              admin_to_child = sd0$admin_to_child[live_a],
              admin_age = sd0$admin_age[live_a],
              log_H = sd0$log_H, a0 = sd0$a0)
  stopifnot(max(new$aa) == new$A, max(new$jj) == new$J,
            length(new$aa) == new$N, length(new$y) == new$N,
            max(new$admin_to_child) <= new$I)
  meta <- b$meta
  meta$cond <- cond; meta$n_items <- new$J
  meta$n_admins_dropped <- sd0$A - new$A
  list(stan_data = new, child_ix = b$child_ix, item_ix = b$item_ix,
       admin_ix = b$admin_ix, meta = meta)
}

for (slug in slugs) {
  bf <- file.path(BL, sprintf("bundle_%s.rds", slug))
  pf <- file.path(SUMM, sprintf("%s_m3_psi.csv", slug))
  if (!file.exists(bf) || !file.exists(pf)) { cat("skip", slug, "(missing bundle or psi)\n"); next }
  b   <- readRDS(bf)
  psi <- read.csv(pf)
  stopifnot(nrow(psi) == b$stan_data$J, !anyNA(psi$delta_j))
  d <- psi$delta_j[order(psi$jj)]          # delta_j indexed by bundle item index
  J <- length(d); o <- order(d)            # item indices, easiest -> hardest
  n50 <- floor(J / 2); n25 <- floor(J / 4)

  set.seed(SEED + sum(utf8ToInt(slug)))
  CONDS <- list(
    mid50  = o[(floor(J * 0.25) + 1):(floor(J * 0.75))],
    mid25  = o[(floor(J * 0.375) + 1):(floor(J * 0.625))],
    easy50 = o[seq_len(n50)],
    hard50 = o[(J - n50 + 1):J],
    rand50 = sample(J, n50),
    rand25 = sample(J, n25))

  cat(sprintf("\n=== %s (J=%d, delta sd=%.2f, range=%.2f) ===\n", slug, J, sd(d), diff(range(d))))
  cat(sprintf("  %-8s %5s %7s %8s %9s %8s\n", "cond", "items", "d_sd", "d_range", "obs", "adm_drop"))
  for (cond in names(CONDS)) {
    keep <- CONDS[[cond]]
    nb <- subset_bundle(b, keep, cond)
    out <- file.path(BL, sprintf("bundle_%s_%s.rds", slug, cond))
    saveRDS(nb, out)
    cat(sprintf("  %-8s %5d %7.2f %8.2f %9d %8d\n", cond, length(keep),
                sd(d[keep]), diff(range(d[keep])), nb$stan_data$N,
                nb$meta$n_admins_dropped))
  }
}
cat("\nwrote bundles into", BL, "\n")
