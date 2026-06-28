## Is input->acceleration a robust English phenomenon, or concentrated / a language thing?
## Recompute the CURRENT both-channel overlap, then fit input->accel SEPARATELY per dataset that
## has item-level CDI + input (RT not needed for this channel). Each: produces ~ input_z*log_age
## + (log_age|child) + (1|item), input_z within-dataset. -> results/input_accel_by_dataset.csv
suppressPackageStartupMessages({library(dplyr); library(lme4); library(here)})
b <- readRDS(here("fits/joint_io_proc_bilingual_subset_data.rds")); d <- b$stan_data; a0 <- d$a0
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
nm <- c("AM2018","FM2012","FMW2013","totlot","babyview","seedlings","SLENA(ES)","HABLA(ES)")
lang <- c("English","English","English","English","English","English","Spanish","Spanish")

## ---- current overlap: input AND RT AND CDI (the historical "both"=142?) ----
item <- unique(d$admin_to_child); rt <- unique(d$lwl_to_child); inp <- unique(d$rec_to_child)
en <- which(d$study_of_child <= 6)
cat(sprintf("CURRENT English overlap: item+input = %d | item+RT = %d | item+input+RT = %d\n",
            sum(en %in% inp & en %in% item), sum(en %in% rt & en %in% item),
            sum(en %in% inp & en %in% rt & en %in% item)))
cat(sprintf("  (glmer 'both_common' was 142 = input∩RT; input→accel only needs item+input)\n"))

fit_one <- function(s) {
  kids <- which(d$study_of_child == s)
  iz <- tibble(ii = d$rec_to_child, li = d$log_input_obs) %>% filter(ii %in% kids) %>%
    group_by(ii) %>% summarise(mli = mean(li), .groups = "drop")
  if (nrow(iz) < 8 || sd(iz$mli) == 0) return(NULL)         # need input variation
  iz <- iz %>% mutate(input_z = as.numeric(scale(mli))) %>% transmute(child = ii, input_z)
  obs <- tibble(child = d$admin_to_child[d$aa], item = d$jj, produces = d$y, age = d$admin_age[d$aa]) %>%
    filter(child %in% iz$child) %>% mutate(log_age = log(age / a0)) %>%
    inner_join(iz, by = "child") %>% mutate(child = factor(child), item = factor(item))
  m <- tryCatch(glmer(produces ~ input_z * log_age + (log_age | child) + (1 | item), obs, binomial, control = ctrl),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients; tm <- "input_z:log_age"
  est <- co[tm, "Estimate"]; se <- co[tm, "Std. Error"]
  data.frame(dataset = nm[s], language = lang[s], n_kids = n_distinct(obs$child),
             est = est, lo = est - 1.96*se, hi = est + 1.96*se, p = co[tm, "Pr(>|z|)"])
}
rows <- lapply(c(1, 3, 5, 6, 7), fit_one)   # datasets w/ item+input: AM2018, FMW2013, babyview, seedlings, SLENA
res <- bind_rows(rows)
write.csv(res, here("studies/io_proc_glmer/results/input_accel_by_dataset.csv"), row.names = FALSE)
cat("\nwrote results/input_accel_by_dataset.csv\n")
res %>% transmute(dataset, language, n_kids, `input->accel` = round(est, 2),
                  ci = sprintf("[%.1f, %.1f]", lo, hi), p = round(p, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
