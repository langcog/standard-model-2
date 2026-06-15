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

## ---- 2. Input-only CDI from io_pooled, restricted to chosen items ---- ##
## BabyView + SEEDLingS = new studies (5,6). PLUS the FMW2013 input-only SALVAGE:
## TLO has 51 LENA kids but proc_dp keeps only the ~31 with RT; the other 20 have
## CDI+input but no RT. io_pooled keeps all 51 (CDI+input). Recover the 20 not in
## proc_dp's fmw_2013 set as input-only kids of the EXISTING study 3 (fmw_2013).
proc_fmw_subj <- as.character(child3$lab_subject_id[child3$dataset_name == "fmw_2013"])
salv_ckeys <- iob$child_info %>%
  filter(s == which(iob$studies$study == "FMW2013"), !subject_id %in% proc_fmw_subj) %>% pull(ckey)
FMW_IDX <- match("fmw_2013", DS3)                       # existing study index (3)
io_cdi <- iob$df %>%
  filter((study %in% c("BabyView", "SEEDLingS")) | (study == "FMW2013" & ckey %in% salv_ckeys),
         item %in% chosen_items) %>%
  ## drop io_pooled's own indices (ii/aa/jj/cc/s) to avoid join collisions
  select(study, ckey, age, item, produces) %>%
  mutate(cid       = paste(study, ckey, sep = "::"),
         study_idx = case_when(study == "BabyView"  ~ length(DS3) + 1L,
                               study == "SEEDLingS" ~ length(DS3) + 2L,
                               study == "FMW2013"   ~ FMW_IDX),   # salvage -> existing study 3
         jj        = unname(item2jj[item]))
cids <- io_cdi %>% distinct(cid, ckey, study, study_idx)
cids$ii <- I3 + seq_len(nrow(cids))
## carry the real subject_id (retained in io_pooled child_info) for robust joins
cids <- cids %>% left_join(iob$child_info %>% distinct(ckey, subject_id), by = "ckey")
io_cdi <- io_cdi %>%
  left_join(cids %>% select(cid, ii), by = "cid") %>%
  mutate(admin_key = paste(cid, age, sep = "_"))
adm_new <- io_cdi %>% distinct(admin_key) %>% mutate(aa = A3 + row_number())
io_cdi <- io_cdi %>% left_join(adm_new, by = "admin_key")
cat(sprintf("added CDI (BabyView+SEEDLingS new + FMW2013 input-only salvage): kids=%d (incl %d salvage) admins=%d rows=%d (items kept=%d/%d chosen)\n",
            nrow(cids), length(salv_ckeys), nrow(adm_new), nrow(io_cdi), n_distinct(io_cdi$item), length(chosen_items)))

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

## ---- 4. Input MEASUREMENT MODEL: RAW per-recording log-input ---- ##
## (mm variant) Source ALL observed input from io_pooled (raw log_r_obs, NOT
## standardized). Map each recording's io child to its JOINT ii via subject_id:
## proc kids match by lab_subject_id; io-added (BV/SE/salvage) by the retained
## subject_id. instr = io measurement type (1=head-cam BabyView, 2=LENA).
## Match on (subject_id AND dataset) -- bare subject_id collides because TotLot
## cohorts (AM2018/totlot) share lab-id ranges. Canonicalize both sides to a
## dataset key, then join on (key, ds).
io_ds <- c(BabyView = "BabyView", SEEDLingS = "SEEDLingS",
           AM2018 = "adams_marchman_2018", FMW2013 = "fmw_2013")
io_ci <- iob$child_info %>%
  mutate(ds = io_ds[iob$studies$study[s]], key = as.character(subject_id)) %>%
  filter(!is.na(key)) %>% transmute(io_ii = ii, key, ds)
joint_keys <- bind_rows(
  child3 %>% transmute(jii = ii, key = as.character(lab_subject_id), ds = dataset_name),
  cids   %>% filter(!is.na(subject_id)) %>%
    transmute(jii = ii, key = as.character(subject_id),
              ds = recode(study, FMW2013 = "fmw_2013")))   # BabyView/SEEDLingS unchanged
io2j <- io_ci %>% inner_join(joint_keys, by = c("key", "ds")) %>% distinct(io_ii, jii)
recs <- tibble::tibble(logr  = iob$stan_data$log_r_obs,
                       io_ii = iob$stan_data$recording_to_child,
                       instr = iob$stan_data$meas_of_study[iob$stan_data$study_of_recording]) %>%
  inner_join(io2j, by = "io_ii")
rec_to_child  <- recs$jii
log_input_obs <- recs$logr
instr_of_rec  <- recs$instr
V_obs   <- nrow(recs)
n_instr <- max(instr_of_rec)
## per-study observed-input level (anchor for mu_r_s); input studies, else overall mean
rec_study <- study_of_child[rec_to_child]
mu_r_s_prior_mean <- rep(mean(log_input_obs), S)
for (s in unique(rec_study)) mu_r_s_prior_mean[s] <- mean(log_input_obs[rec_study == s])
ikc <- tibble::tibble(jii = rec_to_child, s = rec_study) %>% distinct(jii, s) %>% count(s)
cat(sprintf("input MM: V_obs=%d recordings / %d kids; n_instr=%d (1=headcam,2=LENA); kids/study: %s\n",
            V_obs, length(unique(rec_to_child)), n_instr,
            paste(sprintf("s%d=%d", ikc$s, ikc$n), collapse = " ")))

## ---- 4b. SEEDLingS LWL-RT channel (Zhu et al., derived) -> processing study 6 ---- ##
## Per-trial RT from model/scripts/prepare_seedlings_lwl_rt.R, joined by the REAL
## subject_id (NOT the dense io ii). The model reads RT-study via
## study_of_child[lwl_to_child], so these obs auto-route to tau_s[6]/psi_s[6];
## no Stan change. Moves SEEDLingS kids input-only -> both-channel.
se_rt  <- readr::read_csv("data/seedlings/seedlings_lwl_rt.csv", show_col_types = FALSE)
se_map <- cids %>% filter(study == "SEEDLingS") %>% select(subject_id, ii)
se_lwl <- se_rt %>% inner_join(se_map, by = c("lab_subject_id" = "subject_id")) %>%
  transmute(ii, lwl_age, lwl_log_rt)
lwl_to_child_all <- c(sdp$lwl_to_child, se_lwl$ii)
lwl_log_age_all  <- c(sdp$lwl_log_age,  log(se_lwl$lwl_age / a0))
lwl_log_rt_all   <- c(sdp$lwl_log_rt,   se_lwl$lwl_log_rt)
n_se_unmatched <- nrow(se_rt) - nrow(se_lwl)
cat(sprintf("SEEDLingS RT: joined %d obs / %d kids by subject_id (%d RT rows unmatched); N_lwl %d -> %d\n",
            nrow(se_lwl), n_distinct(se_lwl$ii), n_se_unmatched, sdp$N_lwl, length(lwl_log_rt_all)))

## ---- 5. stan_data (mm: input is a measurement model; sigma_r estimated) ---- ##
stan_data <- modifyList(sdp, list(
  N = length(y_all), A = A, I = I, S = S,
  y = y_all, aa = aa_all, jj = jj_all,
  admin_to_child = admin_to_child, admin_age = admin_age,
  study_of_child = study_of_child,
  ## drop the pinned-input fields (sigma_r now a parameter; raw input replaces z)
  sigma_r = NULL, V = NULL, z_lena = NULL, sigma_lena = NULL,
  ## input measurement model: raw per-recording log input + indices + priors
  V_obs = V_obs, rec_to_child = rec_to_child, log_input_obs = log_input_obs,
  n_instr = n_instr, instr_of_rec = instr_of_rec,
  sigma_r_prior_mean = 0.44, sigma_r_prior_sd = 0.10,            # apples-to-apples anchor
  mu_r_s_prior_mean = mu_r_s_prior_mean, mu_r_s_prior_sd = 1.0,  # weak; data identifies
  sigma_meas_prior_sd = 1.0,
  ## RT channel: SEEDLingS (study 6) appended to proc_dp studies 1-4
  grainsize = 1L,
  N_lwl = length(lwl_log_rt_all), lwl_to_child = lwl_to_child_all,
  lwl_log_age = lwl_log_age_all, lwl_log_rt = lwl_log_rt_all,
  ## frank_etal_2026 informative RT priors (full Fig-2 d_sub, centered @ a0=21mo)
  mu_rt_prior_mean = 6.84, mu_rt_prior_sd = 0.20,
  psi_prior_mean = -0.35, mu_rtslope_prior_sd = 0.10,
  sigma_rt0_prior_mean = 0.143, sigma_rt0_prior_sd = 0.04,
  sigma_rt1_prior_mean = 0.26,  sigma_rt1_prior_sd = 0.08,
  sigma_lwl_prior_mean = 0.24,  sigma_lwl_prior_sd = 0.05
))

## child_info for all I (carry subject_id for SEEDLingS/BabyView)
child_new <- cids %>% transmute(ii, lab_subject_id = cid, subject_id,
                                dataset_name = study, study = study_idx)
child_info <- bind_rows(child3, child_new[order(child_new$ii), ]) %>% arrange(ii)

## bundle lwl frame: proc_dp RT + SEEDLingS RT (for inspection/plots)
lwl_all <- bind_rows(lwl, se_lwl %>% mutate(dataset_name = "seedlings_zhu"))

bundle <- list(stan_data = stan_data, word_info = wi, child_info = child_info,
               lwl = lwl_all, datasets = c(DS3, "babyview", "seedlings"),
               input_recs = tibble::tibble(ii = rec_to_child, log_input = log_input_obs, instr = instr_of_rec),
               language = "English (joint io+proc MEASUREMENT-MODEL: +RT, raw input)")
out <- file.path(PATHS$fits_dir, "joint_io_proc_mm_subset_data.rds")
saveRDS(bundle, out)

cat(sprintf("\n== joint io-proc MM bundle ==\n  I=%d A=%d J=%d C=%d S=%d N=%d N_lwl=%d V_obs=%d a0=%d\n",
            I, A, J, C, S, length(y_all), length(lwl_log_rt_all), V_obs, a0))
cat(sprintf("  RT children: %d ; input children: %d (recordings %d) ; both-channel: %d\n",
            length(unique(lwl_to_child_all)), length(unique(rec_to_child)), V_obs,
            length(intersect(unique(lwl_to_child_all), unique(rec_to_child)))))
cat(sprintf("Saved %s\n", out))
