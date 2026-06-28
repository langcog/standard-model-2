## Cheap per-cohort glmer fits to localize the input->acceleration destabilization.
## Same channel grammar as the English ladder (input_z within-cohort, produces ~ input_z*log_age),
## fit separately on the NEW cohorts so we can put them on one coefficient plot:
##   SLENA-item  (29 Spanish kids: item-level CDI + LENA input)  -- the CLEAN test
##   HABLA-sum   (103 Spanish kids: sumscore CDI + input)        -- count-only
##   ELENA-sum   (34 English kids: sumscore CDI + input)         -- count-only
## If clean Spanish item-level recovers a healthy input->accel but the sumscore cohorts don't,
## the destabilizer is the COUNT data (κ unpinned -> σ_ζ inflated), not "Spanish" per se.
## RUN LOCALLY. -> results/cohort_glmer_coefs.csv
suppressPackageStartupMessages({library(dplyr); library(lme4); library(here)})
b <- readRDS(here("fits/joint_io_proc_bilingual_subset_data.rds")); d <- b$stan_data; a0 <- d$a0
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
z <- function(x) as.numeric(scale(x))

input_z_for <- function(kids) tibble(ii = d$rec_to_child, li = d$log_input_obs) %>%
  filter(ii %in% kids) %>% group_by(ii) %>% summarise(mli = mean(li), .groups = "drop") %>%
  mutate(input_z = z(mli)) %>% transmute(child = ii, input_z)

grab <- function(m, n, tag) {
  co <- summary(m)$coefficients; out <- list()
  for (tm in grep("input_z", rownames(co), value = TRUE)) {
    est <- co[tm, "Estimate"]; se <- co[tm, "Std. Error"]
    out[[tm]] <- data.frame(spec = tag, n_kids = n, channel = "input",
      term = if (grepl("log_age", tm)) "acceleration" else "level",
      est = est, lo = est - 1.96*se, hi = est + 1.96*se,
      p = co[tm, "Pr(>|z|)"], model_type = "glmer")
  }
  do.call(rbind, out)
}

res <- list()

## ---- SLENA Spanish: item-level CDI + input (the clean test) ----
sl_kids <- which(d$study_of_child == 7)
iz <- input_z_for(sl_kids)
obs_sl <- tibble(child = d$admin_to_child[d$aa], item = d$jj, produces = d$y, age = d$admin_age[d$aa]) %>%
  filter(child %in% sl_kids, child %in% iz$child) %>%
  mutate(log_age = log(age / a0)) %>% inner_join(iz, by = "child") %>%
  mutate(child = factor(child), item = factor(item))
cat(sprintf("SLENA item-level: %d kids, %d responses\n", n_distinct(obs_sl$child), nrow(obs_sl)))
m_sl <- glmer(produces ~ input_z * log_age + (log_age | child) + (1 | item), obs_sl, binomial, control = ctrl)
res$slena <- grab(m_sl, n_distinct(obs_sl$child), "SLENA-item (ES, clean)")

## ---- sumscore cohorts: cbind(k, n-k) ~ input_z*log_age (frequentist count analog) ----
sum_glmer <- function(study, tag) {
  kids <- which(d$study_of_child == study)
  iz <- input_z_for(kids)
  rows <- tibble(child = d$sum_child, k = d$sum_k, form = d$sum_form, log_age = d$sum_log_age) %>%
    filter(child %in% kids, child %in% iz$child) %>%
    mutate(n = d$form_len[form], fail = n - k) %>% inner_join(iz, by = "child") %>%
    mutate(child = factor(child))
  cat(sprintf("%s: %d kids, %d sumscore admins\n", tag, n_distinct(rows$child), nrow(rows)))
  m <- glmer(cbind(k, fail) ~ input_z * log_age + (log_age | child), rows, binomial, control = ctrl)
  grab(m, n_distinct(rows$child), tag)
}
res$habla <- sum_glmer(8, "HABLA-sum (ES, count)")
res$elena <- sum_glmer(3, "ELENA-sum (EN, count)")   # FMW2013 sumscore kids

coefs <- bind_rows(res)
write.csv(coefs, here("studies/io_proc_glmer/results/cohort_glmer_coefs.csv"), row.names = FALSE)
cat("\nwrote results/cohort_glmer_coefs.csv\n")
coefs %>% filter(term == "acceleration") %>%
  transmute(spec, n_kids, `input->accel` = round(est, 3), ci = sprintf("[%.2f, %.2f]", lo, hi), p = round(p, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
