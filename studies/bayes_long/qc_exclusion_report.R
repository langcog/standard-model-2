## QC exclusion report: spaghetti with unified-outlier-excluded kids/admins in red,
## + a per-dataset count/percentage table. Sources the pull/harmonize/clean logic from
## 00_prepare_bundles.R so the numbers match the bundle build exactly (same order:
## >=MIN_ADMINS pre-filter, item filter, /J proportion, clean_child, re-apply >=MIN_ADMINS).
##
## THRESHOLD. Computed at BOTH 3+ and 2+. The main text says "three or more observations
## per child were included" and its analysis sample is the _a3 bundles, so the exclusion
## percentages it quotes must be the 3+ ones -- previously this script was hard-coded to
## MIN_ADMINS=2, so the prose described 3+ while the numbers were 2+. The 3+ rate is the
## higher of the two: a child with exactly three administrations who loses one to the
## filter falls below the threshold, whereas at 2+ that child is retained.
## `sample` in the cache is the MAIN-TEXT (3+) table; `sample_a2` keeps the 2+ variant for
## the SI, and `qc_spaghetti_data.rds` follows the main-text threshold.
##
## Usage: Rscript studies/bayes_long/qc_exclusion_report.R
suppressPackageStartupMessages({library(wordbankr); library(dplyr); library(tidyr); library(ggplot2); library(here)})
MAIN_MIN_ADMINS <- 3L
## pull the function defs + constants from 00 without running its build loop
src <- readLines(here("studies","bayes_long","00_prepare_bundles.R"))
end <- grep("^## ---- run:", src)[1] - 1          # everything above the main run loop
eval(parse(text = paste(src[1:end], collapse="\n")))

N_SPAG <- 250
## Memoise the Wordbank pull: run_threshold is called once per threshold, and without
## this the script re-pulls every language each time -- twice the network round-trips
## for identical data, on a service that has already been down once this week.
.pull_cache <- new.env(parent = emptyenv())
pull_items_cached <- function(language) {
  key <- gsub("[^A-Za-z]", "_", language)
  if (is.null(.pull_cache[[key]]))
    .pull_cache[[key]] <- harmonize_items(pull_language(language)$items)
  .pull_cache[[key]]
}

run_threshold <- function(MIN_ADMINS) {
set.seed(1)
report <- list(); traj <- list()
for (i in seq_len(nrow(UNITS))) {
  u <- UNITS[i,]
  it <- pull_items_cached(u$language)
  it_u <- if (is.na(u$dataset)) it else filter(it, dataset_name==u$dataset)
  df <- it_u |> group_by(ckey, age, item) |> summarise(produces=max(produces), .groups="drop")
  keep <- df |> distinct(ckey, age) |> count(ckey) |> filter(n>=MIN_ADMINS) |> pull(ckey); df <- filter(df, ckey %in% keep)
  it_keep <- df |> count(item) |> filter(n>=MIN_ITEM_OBS) |> pull(item); df <- filter(df, item %in% it_keep)
  J_qc <- n_distinct(df$item)   # proportion over the full checklist (matches 00_prepare_bundles /J fix)
  prop <- df |> group_by(ckey, age) |> summarise(v=sum(produces)/J_qc, .groups="drop") |> arrange(ckey, age)
  ka <- prop |> group_by(ckey) |> group_modify(~ mutate(.x, keep=clean_child(.x$age, .x$v))) |> ungroup()
  ## a child is "excluded" if any admin removed OR it drops below MIN_ADMINS after cleaning
  surv <- ka |> filter(keep) |> count(ckey) |> filter(n>=MIN_ADMINS) |> pull(ckey)
  ka <- ka |> mutate(excluded = !(ckey %in% surv) | !keep)
  report[[u$slug]] <- tibble(dataset=u$slug,
    kids=n_distinct(ka$ckey), admins=nrow(ka),
    kids_excluded=n_distinct(ka$ckey[!(ka$ckey %in% surv)]),
    admins_removed=sum(!ka$keep))
  ## trajectories for the plot (sample survivors for grey; all excluded kids in red)
  exk <- unique(ka$ckey[!(ka$ckey %in% surv)])
  smp <- sample(surv, min(N_SPAG, length(surv)))
  traj[[u$slug]] <- ka |> filter(ckey %in% c(smp, exk)) |>
    mutate(dataset=u$slug, bad = ckey %in% exk)
}
tab <- bind_rows(report) |> mutate(pct_kids=100*kids_excluded/kids, pct_admins=100*admins_removed/admins)
cat(sprintf("=== QC exclusion table (MIN_ADMINS=%d, unified local-outlier filter) ===\n", MIN_ADMINS))
print(as.data.frame(tab), digits=3)
TR <- bind_rows(traj) |> mutate(dataset=factor(dataset, levels=UNITS$slug))
list(tab = tab, TR = TR)
}

R3 <- run_threshold(MAIN_MIN_ADMINS)      # main text
R2 <- run_threshold(2L)                   # SI / full longitudinal sample
tab <- R3$tab; TR <- R3$TR
saveRDS(tab, here("studies","bayes_long","qc_exclusion_table.rds"))

## ---- paper cache: sample/exclusion table (methods) + spaghetti data (SI figure) ----
## Self-contained (bundles only), so it renders without the still-running fits.
LAB <- c(thal="English (Thal)", smith="English (Smith)", marchman="English (Marchman)",
         norwegian="Norwegian", japanese="Japanese")
## Read the bundle MATCHING the threshold: _a3 for 3+, unsuffixed for 2+. Using the 2+
## bundle for a 3+ exclusion table was the original mismatch.
make_sample <- function(tabx, min_adm) {
  sfx <- if (min_adm >= 3L) "_a3" else ""
  samp <- bind_rows(lapply(UNITS$slug, function(s){
    m <- readRDS(here("fits","bayes_long", paste0("bundle_",s,sfx,".rds")))$meta
    e <- tabx[tabx$dataset==s,]
    tibble(dataset=s, label=LAB[[s]], n_kids_raw=e$kids, n_kids=m$n_kids, n_admins=m$n_admins,
           n_obs=m$n_obs, age_lo=m$age_range[1], age_hi=m$age_range[2],
           kids_excluded=e$kids_excluded, admins_removed=e$admins_removed,
           pct_kids_excluded=100*e$kids_excluded/e$kids)
  }))
  bind_rows(samp, summarise(samp, dataset="total", label="All", n_kids_raw=sum(n_kids_raw),
    n_kids=sum(n_kids), n_admins=sum(n_admins), n_obs=sum(n_obs), age_lo=min(age_lo), age_hi=max(age_hi),
    kids_excluded=sum(kids_excluded), admins_removed=sum(admins_removed),
    pct_kids_excluded=100*sum(kids_excluded)/sum(n_kids_raw)))
}
samp    <- make_sample(R3$tab, MAIN_MIN_ADMINS)   # main text (3+)
samp_a2 <- make_sample(R2$tab, 2L)                # SI (full longitudinal sample)
## consistency check: the retained counts must equal the fitted _a3 analysis sample
stopifnot(all(samp$n_kids[samp$dataset != "total"] ==
              vapply(UNITS$slug, function(s)
                readRDS(here("fits","bayes_long", paste0("bundle_",s,"_a3.rds")))$stan_data$I, integer(1))))
saveRDS(list(sample=samp, sample_a2=samp_a2, min_admins=MAIN_MIN_ADMINS,
             qc_rule=list(rel_tol=0.25, peak_floor=0.10, drop_floor=0.05, rate_max=0.40, jump_base=0.10)),
        here("paper","cache","bayes_long_sample.rds"))
saveRDS(TR, here("paper","cache","qc_spaghetti_data.rds"))
cat("wrote paper/cache/{bayes_long_sample,qc_spaghetti_data}.rds\n")
p <- ggplot() +
  geom_line(data=filter(TR,!bad), aes(age, v, group=ckey), color="grey70", alpha=0.20, linewidth=0.25) +
  geom_line(data=filter(TR, bad), aes(age, v, group=ckey), color="#d7191c", alpha=0.70, linewidth=0.4) +
  geom_point(data=filter(TR, bad, !keep), aes(age, v), color="#d7191c", size=0.7) +
  facet_wrap(~dataset, nrow=1) + scale_y_continuous(limits=c(0,1)) +
  labs(x="Age (months)", y="proportion of items produced",
       title="Unified local-outlier QC: red = excluded children (dots = removed administrations)",
       subtitle=sprintf("MIN_ADMINS=%d (main-text analysis sample); crater(>25%% below peak) or jump(>0.40/mo from base<0.10)", MAIN_MIN_ADMINS)) +
  theme_minimal(base_size=10) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank(), plot.title=element_text(size=11,face="bold"))
out <- here("studies","bayes_long","qc_exclusion_spaghetti.png")
ggsave(out, p, width=3.1*nrow(UNITS)+0.5, height=3.6, dpi=150); cat("wrote", out, "\n")
