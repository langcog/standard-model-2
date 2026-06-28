## SLENA Spanish processing->efficiency (the missing Spanish point for the comparison forest).
## Same proc_z construction as the English ladder: per-child mean RT residual, flipped, z-scored.
## If Spanish proc->eff lands normal (~0.6, positive) while Spanish input->accel is negative,
## that points to a channel-SPECIFIC language difference, not generic SLENA noise.
## RUN LOCALLY. -> results/slena_proc.csv
suppressPackageStartupMessages({library(dplyr); library(lme4); library(here)})
b <- readRDS(here("fits/joint_io_proc_bilingual_subset_data.rds")); d <- b$stan_data; a0 <- d$a0
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
sl <- which(d$study_of_child == 7)

## proc_z: per-child mean residual from a simple RT measurement model, flipped (higher=faster)
rt <- tibble(ii = d$lwl_to_child, log_age = d$lwl_log_age, log_rt = d$lwl_log_rt) %>% filter(ii %in% sl)
mrt <- lm(log_rt ~ log_age, rt); rt$res <- residuals(mrt)
proc <- rt %>% group_by(ii) %>% summarise(pr = mean(res), .groups = "drop") %>%
  mutate(proc_z = as.numeric(scale(-pr))) %>% transmute(child = ii, proc_z)
obs <- tibble(child = d$admin_to_child[d$aa], item = d$jj, produces = d$y, age = d$admin_age[d$aa]) %>%
  filter(child %in% proc$child) %>% mutate(log_age = log(age / a0)) %>%
  inner_join(proc, by = "child") %>% mutate(child = factor(child), item = factor(item))
cat(sprintf("SLENA proc: %d kids, %d responses\n", n_distinct(obs$child), nrow(obs)))
m <- glmer(produces ~ proc_z * log_age + (log_age | child) + (1 | item), obs, binomial, control = ctrl)
co <- summary(m)$coefficients
out <- do.call(rbind, lapply(c("proc_z", "proc_z:log_age"), function(tm) {
  est <- co[tm, "Estimate"]; se <- co[tm, "Std. Error"]
  data.frame(spec = "SLENA-proc", n_kids = n_distinct(obs$child), channel = "processing",
             term = if (grepl("log_age", tm)) "acceleration" else "level",
             est = est, lo = est - 1.96*se, hi = est + 1.96*se, p = co[tm, "Pr(>|z|)"])
}))
write.csv(out, here("studies/io_proc_glmer/results/slena_proc.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
