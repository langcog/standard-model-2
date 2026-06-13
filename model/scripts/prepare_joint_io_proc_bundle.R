## Build the JOINT input + processing bundle (log_irt_long_proc_dp_joint.stan).
##
## Extends the validated proc_dp bundle (AM2018 + FM2012 + FMW2013: CDI + LWL-RT
## + observed LENA input) with the two input-only datasets BabyView + SEEDLingS
## (CDI + observed input, NO RT) pulled from io_pooled_subset_data.rds. The two
## are disjoint from proc_dp's datasets, so children are simply appended (new ii)
## -- no cross-source ID reconciliation. Items are restricted to the proc_dp
## chosen set so J is shared. sigma_lena becomes per-study (headcam != LENA).
##
## RUN LOCALLY.  Output: fits/joint_io_proc_subset_data.rds
source("model/R/config.R"); source("model/R/helpers.R")
suppressPackageStartupMessages({ library(dplyr) })

pd  <- readRDS(file.path(PATHS$fits_dir, "proc_dp_all_subset_data.rds"))
iob <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))

## ---- 1. proc_dp pieces (studies 1-3) ---- ##
cdi3   <- pd$df
wi     <- pd$word_info                       # jj -> item, prob, cc
child3 <- pd$child_info                       # ii -> lab_subject_id, dataset_name, study
lwl    <- pd$lwl                              # ii, lwl_age, lwl_log_rt
lena3  <- pd$lena                             # ii, z_lena (one row per input child)
sdp    <- pd$stan_data
I3 <- sdp$I; A3 <- sdp$A; J <- sdp$J; C <- sdp$C
DS3 <- pd$datasets                            # AM2018, FM2012, FMW2013  (studies 1,2,3)
a0  <- sdp$a0
item2jj <- setNames(wi$jj, wi$item)
chosen_items <- wi$item
cat(sprintf("proc_dp base: I=%d A=%d J=%d S=%d (datasets: %s)\n",
            I3, A3, J, length(DS3), paste(DS3, collapse=", ")))

## ---- 2. BabyView + SEEDLingS CDI from io_pooled, restricted to chosen items ---- ##
io_cdi <- iob$df %>%
  filter(study %in% c("BabyView", "SEEDLingS"), item %in% chosen_items) %>%
  ## drop io_pooled's own indices (ii/aa/jj/cc/s) to avoid join collisions
  select(study, ckey, age, item, produces) %>%
  mutate(cid       = paste(study, ckey, sep = "::"),
         study_idx = ifelse(study == "BabyView", length(DS3) + 1L, length(DS3) + 2L),
         jj        = unname(item2jj[item]))
cids <- io_cdi %>% distinct(cid, study, study_idx)
cids$ii <- I3 + seq_len(nrow(cids))
io_cdi <- io_cdi %>%
  left_join(cids %>% select(cid, ii), by = "cid") %>%
  mutate(admin_key = paste(cid, age, sep = "_"))
adm_new <- io_cdi %>% distinct(admin_key) %>% mutate(aa = A3 + row_number())
io_cdi <- io_cdi %>% left_join(adm_new, by = "admin_key")
cat(sprintf("added CDI: BabyView+SEEDLingS kids=%d admins=%d rows=%d (items kept=%d/%d chosen)\n",
            nrow(cids), nrow(adm_new), nrow(io_cdi), n_distinct(io_cdi$item), length(chosen_items)))

## ---- 3. Combined CDI long frame (proc_dp rows + new rows) ---- ##
y_all  <- c(sdp$y, io_cdi$produces)
aa_all <- c(sdp$aa, io_cdi$aa)
jj_all <- c(sdp$jj, io_cdi$jj)

## admin_to_child + admin_age for the new admins
adm_child_new <- io_cdi %>% distinct(aa, ii, age) %>% arrange(aa)
admin_to_child <- c(sdp$admin_to_child, adm_child_new$ii)
admin_age      <- c(sdp$admin_age,      adm_child_new$age)

## study_of_child for all I = I3 + new
study_of_child <- c(sdp$study_of_child, cids$study_idx[order(cids$ii)])
I <- I3 + nrow(cids); A <- A3 + nrow(adm_new); S <- length(DS3) + 2L

## ---- 4. BabyView + SEEDLingS observed input (one standardized z per child) ---- ##
## recover per-recording log input from io_pooled stan_data, map io-ii -> ckey
ii2ckey <- iob$df %>% distinct(ii, ckey, study)
rec <- tibble::tibble(logr = iob$stan_data$log_r_obs, io_ii = iob$stan_data$recording_to_child) %>%
  left_join(ii2ckey, by = c("io_ii" = "ii")) %>%
  filter(study %in% c("BabyView", "SEEDLingS")) %>%
  mutate(cid = paste(study, ckey, sep = "::"))
## per-child mean log input + within-child noise var of that mean
perchild <- rec %>% group_by(study, cid) %>%
  summarise(lmean = mean(logr), n_rec = n(),
            wvar  = ifelse(n() >= 2, var(logr), NA_real_), .groups = "drop")
## standardize within study (noise-corrected), per-study sigma_lena
add_in <- perchild %>% group_by(study) %>%
  group_modify(function(g, key) {
    noise_var_mean <- mean(g$wvar / g$n_rec, na.rm = TRUE)         # avg var of the child-mean read
    sd_obs  <- sd(g$lmean)
    sd_true <- sqrt(pmax(sd_obs^2 - noise_var_mean, 1e-4))
    g$z_lena      <- (g$lmean - mean(g$lmean)) / sd_true
    g$sigma_lena  <- sqrt(noise_var_mean) / sd_true
    g
  }) %>% ungroup() %>%
  left_join(cids %>% select(cid, ii, study_idx), by = "cid")

## per-study sigma_lena vector (study order 1..5)
sigma_lena_vec <- numeric(S)
## proc_dp studies (1..length(DS3)): the pooled LENA sd for the input-having
## ones (those represented among proc_dp's input children, lena3), placeholder
## for the RT-only ones (FM2012, fernald_totlot -- no observed input).
pd_input_studies <- unique(child3$study[match(lena3$ii, child3$ii)])
for (s in seq_len(length(DS3)))
  sigma_lena_vec[s] <- if (s %in% pd_input_studies) sdp$sigma_lena else 1.0
sigma_lena_vec[length(DS3) + 1] <- add_in$sigma_lena[add_in$study == "BabyView"][1]
sigma_lena_vec[length(DS3) + 2] <- add_in$sigma_lena[add_in$study == "SEEDLingS"][1]

## combined input arrays: proc_dp input children (z_lena3) + new
rec_to_child <- c(lena3$ii, add_in$ii)
z_lena       <- c(lena3$z_lena, add_in$z_lena)
V <- length(z_lena)
cat(sprintf("input channel: V=%d (proc_dp %d + BV/SE %d). sigma_lena per study: %s\n",
            V, nrow(lena3), nrow(add_in), paste(sprintf("%.3f", sigma_lena_vec), collapse=", ")))

## ---- 5. stan_data (same shape as proc_dp; sigma_lena now vector[S]) ---- ##
## Pin sigma_r to the apples-to-apples value (channel-matched US English
## between-child log-input SD; validated by the GCP sigma_r pins). The input
## share scales as sigma_r^2 (var_input_xi = sigma_r^2), so use 0.44 here too
## rather than proc_dp's inherited 0.53 Sperry-CDS value.
SIGMA_R_PIN <- 0.44
stan_data <- modifyList(sdp, list(
  N = length(y_all), A = A, I = I, S = S,
  y = y_all, aa = aa_all, jj = jj_all,
  admin_to_child = admin_to_child, admin_age = admin_age,
  study_of_child = study_of_child,
  sigma_r = SIGMA_R_PIN,
  V = V, rec_to_child = rec_to_child, z_lena = z_lena,
  sigma_lena = sigma_lena_vec
))

## child_info for all I
child_new <- cids %>% transmute(ii, lab_subject_id = cid,
                                dataset_name = study, study = study_idx)
child_info <- bind_rows(child3, child_new[order(child_new$ii), ]) %>% arrange(ii)

bundle <- list(stan_data = stan_data, word_info = wi, child_info = child_info,
               lwl = lwl, datasets = c(DS3, "babyview", "seedlings"),
               sigma_lena_by_study = sigma_lena_vec,
               language = "English (joint io+proc: AM2018, FM2012, FMW2013, BabyView, SEEDLingS)")
out <- file.path(PATHS$fits_dir, "joint_io_proc_subset_data.rds")
saveRDS(bundle, out)

cat(sprintf("\n== joint bundle ==\n  I=%d A=%d J=%d C=%d S=%d N=%d N_lwl=%d V=%d a0=%d\n",
            I, A, J, C, S, length(y_all), sdp$N_lwl, V, a0))
cat(sprintf("  RT children: %d (studies 1-%d) ; input children: %d ; both: %d\n",
            length(unique(sdp$lwl_to_child)), length(DS3), V,
            length(intersect(unique(sdp$lwl_to_child), rec_to_child))))
cat(sprintf("Saved %s\n", out))
