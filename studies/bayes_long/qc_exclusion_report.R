## QC exclusion report: spaghetti with unified-outlier-excluded kids/admins in red,
## + a per-dataset count/percentage table. Run at MIN_ADMINS=2 (the full longitudinal
## sample) so nothing is hidden. Sources the pull/harmonize/clean logic from
## 00_prepare_bundles.R so the numbers match the bundle build exactly.
##
## Usage: Rscript studies/bayes_long/qc_exclusion_report.R
suppressPackageStartupMessages({library(wordbankr); library(dplyr); library(tidyr); library(ggplot2); library(here)})
Sys.setenv(MIN_ADMINS="2")
## pull the function defs + constants from 00 without running its build loop
src <- readLines(here("studies","bayes_long","00_prepare_bundles.R"))
end <- grep("^## ---- run:", src)[1] - 1          # everything above the main run loop
eval(parse(text = paste(src[1:end], collapse="\n")))

set.seed(1); N_SPAG <- 250
report <- list(); traj <- list()
for (i in seq_len(nrow(UNITS))) {
  u <- UNITS[i,]
  pl <- pull_language(u$language); it <- harmonize_items(pl$items)
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
cat("=== QC exclusion table (MIN_ADMINS=2, unified local-outlier filter) ===\n"); print(as.data.frame(tab), digits=3)
saveRDS(tab, here("studies","bayes_long","qc_exclusion_table.rds"))

TR <- bind_rows(traj) |> mutate(dataset=factor(dataset, levels=UNITS$slug))

## ---- paper cache: sample/exclusion table (methods) + spaghetti data (SI figure) ----
## Self-contained (bundles only), so it renders without the still-running fits.
LAB <- c(thal="English (Thal)", smith="English (Smith)", marchman="English (Marchman)",
         norwegian="Norwegian", japanese="Japanese")
samp <- bind_rows(lapply(UNITS$slug, function(s){
  m <- readRDS(here("fits","bayes_long", paste0("bundle_",s,".rds")))$meta
  e <- tab[tab$dataset==s,]
  tibble(dataset=s, label=LAB[[s]], n_kids_raw=e$kids, n_kids=m$n_kids, n_admins=m$n_admins,
         n_obs=m$n_obs, age_lo=m$age_range[1], age_hi=m$age_range[2],
         kids_excluded=e$kids_excluded, admins_removed=e$admins_removed,
         pct_kids_excluded=100*e$kids_excluded/e$kids)
}))
samp <- bind_rows(samp, summarise(samp, dataset="total", label="All", n_kids_raw=sum(n_kids_raw),
  n_kids=sum(n_kids), n_admins=sum(n_admins), n_obs=sum(n_obs), age_lo=min(age_lo), age_hi=max(age_hi),
  kids_excluded=sum(kids_excluded), admins_removed=sum(admins_removed),
  pct_kids_excluded=100*sum(kids_excluded)/sum(n_kids_raw)))
saveRDS(list(sample=samp, min_admins=2,
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
       subtitle="MIN_ADMINS=2 (full longitudinal sample); crater(>25% below peak) or jump(>0.40/mo from base<0.10)") +
  theme_minimal(base_size=10) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank(), plot.title=element_text(size=11,face="bold"))
out <- here("studies","bayes_long","qc_exclusion_spaghetti.png")
ggsave(out, p, width=3.1*nrow(UNITS)+0.5, height=3.6, dpi=150); cat("wrote", out, "\n")
