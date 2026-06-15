## Reproducible glmer ladder: simple per-child measurement summaries as an
## "unbiased" benchmark for the Bayesian io-proc model. Fits 5 models and writes
## a tidy coefficient table (Wald 95% CIs). RUN LOCALLY (lme4, ~40 min; the
## processing-full fit is 689k obs). Caches models so re-runs are cheap.
suppressPackageStartupMessages({library(dplyr); library(lme4); library(here)})
b  <- readRDS(here("fits/joint_io_proc_mm_subset_data.rds")); sd <- b$stan_data; a0 <- sd$a0
zstudy <- function(x, s) ave(x, s, FUN = function(v) as.numeric(scale(v)))   # within-study z
ctrl   <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))

## ---- processing_z: simple RT measurement model, per-kid mean residual ----
rt   <- b$lwl %>% mutate(log_age = log(lwl_age / a0), dataset = dataset_name)
mrt  <- lmer(lwl_log_rt ~ log_age + (1 | dataset), rt, REML = FALSE)
rt$resid <- residuals(mrt)
proc <- rt %>% group_by(ii) %>% summarise(presid = mean(resid), .groups = "drop") %>%
  mutate(study = sd$study_of_child[ii], proc_z = zstudy(-presid, study))   # higher = FASTER
## ---- input_z: mean observed log input per kid, within-study z ----
inp <- tibble(ii = sd$rec_to_child, li = sd$log_input_obs) %>%
  group_by(ii) %>% summarise(mli = mean(li), .groups = "drop") %>%
  mutate(study = sd$study_of_child[ii], input_z = zstudy(mli, study))
## ---- full vocab obs (all kids) ----
obs <- tibble(child = sd$admin_to_child[sd$aa], item = sd$jj, produces = sd$y,
              age = sd$admin_age[sd$aa]) %>%
  mutate(log_age = log(age / a0), dataset = factor(sd$study_of_child[child]),
         child = factor(child), item = factor(item)) %>%
  left_join(inp  %>% transmute(child = factor(ii), input_z), by = "child") %>%
  left_join(proc %>% transmute(child = factor(ii), proc_z),  by = "child")
both <- intersect(unique(inp$ii), unique(proc$ii))
common <- obs %>% filter(child %in% factor(both))

fit_cache <- function(tag, form, data) {
  f <- here("studies/io_proc_glmer/cache", paste0(tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  m <- glmer(form, data, binomial, control = ctrl); saveRDS(m, f); m
}
F_log <- produces ~ input_z * log_age + (log_age | child) + (1 | item) + (1 | dataset)
models <- list(
  input_common = list(fit_cache("input_common", produces ~ input_z*log_age + (log_age|child)+(1|item)+(1|dataset), common), 142),
  proc_common  = list(fit_cache("proc_common",  produces ~ proc_z*log_age  + (log_age|child)+(1|item)+(1|dataset), common), 142),
  both_common  = list(fit_cache("both_common",  produces ~ input_z*log_age + proc_z*log_age + (log_age|child)+(1|item)+(1|dataset), common), 142),
  input_full   = list(fit_cache("input_full",   produces ~ input_z*log_age + (log_age|child)+(1|item)+(1|dataset), filter(obs, !is.na(input_z))), 193),
  proc_full    = list(fit_cache("proc_full",    produces ~ proc_z*log_age  + (log_age|child)+(1|item)+(1|dataset), filter(obs, !is.na(proc_z))), 326)
)

## ---- tidy coefficient table (Wald 95% CIs) ----
rows <- lapply(names(models), function(tag) {
  m <- models[[tag]][[1]]; n <- models[[tag]][[2]]; co <- summary(m)$coefficients
  keep <- grep("input_z|proc_z", rownames(co), value = TRUE)   # catches input_z:log_age AND log_age:proc_z
  do.call(rbind, lapply(keep, function(tm) {
    est <- co[tm,"Estimate"]; se <- co[tm,"Std. Error"]
    chan <- if (grepl("input", tm)) "input" else "processing"
    term <- if (grepl("log_age", tm)) "acceleration" else "level"
    data.frame(spec = tag, n_kids = n, channel = chan, term = term,
               est = est, lo = est - 1.96*se, hi = est + 1.96*se, p = co[tm,"Pr(>|z|)"],
               model_type = "glmer")
  }))
})
coefs <- do.call(rbind, rows)
write.csv(coefs, here("studies/io_proc_glmer/results/glmer_ladder_coefs.csv"), row.names = FALSE)
cat("wrote results/glmer_ladder_coefs.csv\n"); print(coefs, digits = 3)
