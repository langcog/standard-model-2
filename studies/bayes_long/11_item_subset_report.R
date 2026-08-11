## 11_item_subset_report.R -- does kappa track the learner, or the ruler?
##
## Reads the item-narrowed refits built by 09_item_subset_bundles.R and asks whether the
## recovered kappa follows the retained item-difficulty spread.
##
## THE TWO HYPOTHESES, MADE QUANTITATIVE IN ADVANCE.
##   ruler   kappa is a rescaling of the instrument, so it is proportional to the retained
##           difficulty spread: kappa_cond = kappa_all * sd_cond / sd_all.
##   learner kappa is a property of the children, so it is invariant: kappa_cond = kappa_all.
## Every condition is scored against both, and the size-matched random controls say how
## much of any movement is simply the loss of half the data.
##
## Usage:  Rscript studies/bayes_long/11_item_subset_report.R
## Output: paper/cache/si_item_subset.rds + console report

suppressPackageStartupMessages({library(dplyr)})
BL <- file.path("fits","bayes_long"); SUMM <- file.path(BL, "summaries")
CACHE <- file.path("paper","cache")
DS <- c(thal_a3 = "English (Thal)", marchman_a3 = "English (Marchman)")
CONDS <- c("mid50","mid25","easy50","hard50","rand50","rand25")
NARROWED <- c("mid50","mid25","easy50","hard50")

get1 <- function(slug, cond) {
  sf <- if (is.na(cond)) file.path(SUMM, sprintf("%s_m3.summary.rds", slug))
        else             file.path(SUMM, sprintf("%s_%s_m3.summary.rds", slug, cond))
  bf <- if (is.na(cond)) file.path(BL, sprintf("bundle_%s.rds", slug))
        else             file.path(BL, sprintf("bundle_%s_%s.rds", slug, cond))
  if (!file.exists(sf) || !file.exists(bf)) return(NULL)
  s <- as.data.frame(readRDS(sf))
  data.frame(slug = slug, cond = if (is.na(cond)) "all" else cond,
             kappa = s$median[s$variable == "kappa_pop"],
             kappa_lo = s$q5[s$variable == "kappa_pop"],
             kappa_hi = s$q95[s$variable == "kappa_pop"],
             rhat = s$rhat[s$variable == "kappa_pop"],
             ess = s$ess_bulk[s$variable == "kappa_pop"],
             row.names = NULL)
}

## The retained item set is not recoverable from the subset bundle's reindexed jj, so
## rebuild the cuts exactly as 09_item_subset_bundles.R did, from the same baseline
## delta_j and the same seed.
retained_sd <- function(slug) {
  psi <- read.csv(file.path(SUMM, sprintf("%s_m3_psi.csv", slug)))
  d <- psi$delta_j[order(psi$jj)]; J <- length(d); o <- order(d)
  n50 <- floor(J/2); n25 <- floor(J/4)
  set.seed(20260811L + sum(utf8ToInt(slug)))
  cuts <- list(mid50  = o[(floor(J*0.25)+1):(floor(J*0.75))],
               mid25  = o[(floor(J*0.375)+1):(floor(J*0.625))],
               easy50 = o[seq_len(n50)],
               hard50 = o[(J-n50+1):J],
               rand50 = sample(J, n50),
               rand25 = sample(J, n25))
  c(list(all = data.frame(cond="all", d_sd=sd(d), n_items=J)),
    lapply(names(cuts), function(cn)
      data.frame(cond=cn, d_sd=sd(d[cuts[[cn]]]), n_items=length(cuts[[cn]])))) |>
    bind_rows()
}

## `sl_` rather than `slug`: inside mutate() the bare name `slug` resolves to the data
## frame's own slug COLUMN, not to the loop variable, and DS[[<vector>]] then errors.
res <- bind_rows(lapply(names(DS), function(sl_) {
  k <- bind_rows(c(list(get1(sl_, NA)), lapply(CONDS, function(c) get1(sl_, c))))
  if (is.null(k) || !nrow(k)) return(NULL)
  left_join(k, retained_sd(sl_), by = "cond") |>
    mutate(lang = DS[[sl_]],
           sd_ratio   = d_sd / d_sd[cond == "all"],
           kappa_all  = kappa[cond == "all"],
           pred_ruler = kappa_all * sd_ratio,          # proportional to retained spread
           pred_learn = kappa_all,                      # invariant
           kappa_ratio = kappa / kappa_all)
}))
if (!nrow(res)) { cat("no completed subset fits\n"); quit(save="no") }

dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(CACHE, "si_item_subset.rds"))

for (sl in names(DS)) {
  r <- filter(res, slug == sl)
  cat(sprintf("\n=== %s (baseline kappa %.2f, delta sd %.2f) ===\n",
              DS[[sl]], r$kappa_all[1], r$d_sd[r$cond=="all"]))
  cat(sprintf("  %-8s %6s %7s %8s %9s %9s %7s\n",
              "cond","items","d_sd","kappa","ruler pred","learner","kappa/all"))
  for (i in seq_len(nrow(r))) with(r[i,],
    cat(sprintf("  %-8s %6d %7.2f %8.2f %9.2f %9.2f %7.2f\n",
                cond, n_items, d_sd, kappa, pred_ruler, pred_learn, kappa_ratio)))
}

## Scoring: across the narrowed conditions, is kappa closer to proportional or invariant?
cat("\n=== which hypothesis fits better across the narrowed conditions? ===\n")
for (sl in names(DS)) {
  r <- filter(res, slug == sl, cond %in% NARROWED)
  e_ruler <- sqrt(mean((r$kappa - r$pred_ruler)^2))
  e_learn <- sqrt(mean((r$kappa - r$pred_learn)^2))
  cat(sprintf("  %-20s RMSE vs ruler %6.2f | vs learner %6.2f  -> %s\n",
              DS[[sl]], e_ruler, e_learn,
              if (e_learn < e_ruler) "LEARNER" else "ruler"))
}
cat("\n=== the controls: how much is just losing half the data? ===\n")
for (sl in names(DS)) {
  r <- filter(res, slug == sl)
  g <- function(c) r$kappa_ratio[r$cond == c]
  cat(sprintf("  %-20s rand50 %.2f  rand25 %.2f   (vs mid50 %.2f, mid25 %.2f)\n",
              DS[[sl]], g("rand50"), g("rand25"), g("mid50"), g("mid25")))
}
cat("\nwrote", file.path(CACHE, "si_item_subset.rds"), "\n")
