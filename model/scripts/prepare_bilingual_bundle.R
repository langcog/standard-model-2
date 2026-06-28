## Build the BILINGUAL ("bi-lean") bundle: the full English io-proc MM bundle
## (fits/joint_io_proc_mm_subset_data.rds) AUGMENTED with the Spanish + sumscore extension.
##
## Adds, on top of the shared population geometry:
##   - Spanish item-level CDI (SLENA / WF2013), canonicalized to Wordbank-ES via the
##     cdi_short_code_map_es_{wg,ws}.csv maps -> a NON-OVERLAPPING Spanish item bank
##     (one new lexical class cc = C_en+1; shared mu_mu_c hyperprior softly ties scales).
##   - Spanish LWL (SLENA n29 + HABLA DRT) and Spanish input (SLENA LENA + HABLA AWC),
##     routed to two new Spanish studies (WF2013, HABLA) via study_of_child.
##   - The SUMSCORE count branch: ELENA-WS (English, attaches to existing fmw kids as a
##     2nd timepoint where ids match) + HABLA-WG/WS (Spanish, count-only). Each admin is
##     k of n produced over its form's item set (English-WS / Spanish-WG / Spanish-WS).
##
## RUN LOCALLY.  Output: fits/joint_io_proc_bilingual_subset_data.rds
suppressPackageStartupMessages({library(here); library(dplyr); library(readr); library(readxl); library(tidyr)})
RAW <- function(...) here("data/raw", ...)

b  <- readRDS(here("fits/joint_io_proc_mm_subset_data.rds"))
sd <- b$stan_data
I0 <- sd$I; A0 <- sd$A; J0 <- sd$J; C0 <- sd$C; S0 <- sd$S; a0 <- sd$a0
ci <- b$child_info
cat(sprintf("English base: I=%d A=%d J=%d C=%d S=%d N=%d N_lwl=%d V_obs=%d a0=%g\n",
            I0, A0, J0, C0, S0, sd$N, sd$N_lwl, sd$V_obs, a0))

## ===== 1. Spanish item bank (canonical, non-overlapping with English) =====
mwg <- read_csv(here("data/intermediates/cdi_short_code_map_es_wg.csv"), show_col_types = FALSE) %>%
  filter(!is.na(item_definition)) %>% transmute(form = "WG", short, def = item_definition)
mws <- read_csv(here("data/intermediates/cdi_short_code_map_es_ws.csv"), show_col_types = FALSE) %>%
  filter(!is.na(item_definition)) %>% transmute(form = "WS", short, def = item_definition)
codemap <- bind_rows(mwg, mws)                                  # (form, short) -> canonical def
sp_items <- tibble(def = sort(unique(codemap$def))) %>% mutate(jj = J0 + row_number())
J_es <- nrow(sp_items)
SP_CC <- C0 + 1L                                               # single new Spanish class
cat(sprintf("Spanish item bank: %d canonical items (jj %d..%d), class cc=%d\n",
            J_es, J0 + 1L, J0 + J_es, SP_CC))

## ===== 2. SLENA Spanish item-level CDI -> long (canonical jj) =====
sl <- read_csv(here("data/intermediates/slena_cdi_items_long.csv"), show_col_types = FALSE) %>%
  mutate(subject_id = as.character(subject_id)) %>%
  inner_join(codemap, by = c("form", "item" = "short")) %>%
  left_join(sp_items, by = "def") %>%
  group_by(subject_id, age, form, jj) %>% summarise(produces = max(produces), .groups = "drop")  # dedup canonical collisions
slena_kids <- tibble(subject_id = sort(unique(sl$subject_id)))  # 42 Spanish item-level kids
cat(sprintf("SLENA item-level: %d kids, %d (kid,age,form) admins, %d responses\n",
            nrow(slena_kids), n_distinct(paste(sl$subject_id, sl$age, sl$form)), nrow(sl)))

## ===== 3. HABLA Spanish kids (count + LWL + input, no item-level) =====
hb <- read_csv(RAW("Bang2025/Habla1.0_LENALWLCDISumScores.csv"), show_col_types = FALSE) %>%
  mutate(ID = as.character(ID))
habla_kids <- tibble(subject_id = sort(unique(hb$ID)))

## ===== 4. Assign new ii: 10 ELENA-new English (study 3) + SLENA (study 7) + HABLA (study 8) =====
el <- read_excel(RAW("FMW2013/elena/ELENA_WS_SumScores.xlsx")) %>%
  transmute(subject_id = as.character(ParticipantId), age = as.numeric(CDIAge), k = as.numeric(VOCAB)) %>%
  filter(!is.na(k), !is.na(age))
eng_fmw_ids <- ci %>% filter(dataset_name == "fmw_2013") %>% transmute(subject_id = as.character(lab_subject_id), ii)
el <- el %>% left_join(eng_fmw_ids, by = "subject_id")          # ii where matched, NA -> new English kid
elena_new_ids <- sort(unique(el$subject_id[is.na(el$ii)]))      # 10 count-only English kids

WF_STUDY <- S0 + 1L; HB_STUDY <- S0 + 2L                       # 7 = WF2013, 8 = HABLA
new_kids <- bind_rows(
  tibble(subject_id = elena_new_ids,       study = 3L,        dataset_name = "fmw_2013"),
  tibble(subject_id = slena_kids$subject_id, study = WF_STUDY, dataset_name = "WF2013"),
  tibble(subject_id = habla_kids$subject_id, study = HB_STUDY, dataset_name = "Bang2025")
) %>% mutate(ii = I0 + row_number())
I <- I0 + nrow(new_kids); S <- HB_STUDY
study_of_child <- c(sd$study_of_child, new_kids$study[order(new_kids$ii)])
## resolve ELENA ii (existing-or-new) now that new English kids have ii
el <- el %>% left_join(new_kids %>% filter(dataset_name == "fmw_2013") %>% select(subject_id, ii_new = ii),
                       by = "subject_id") %>%
  mutate(ii = coalesce(ii, ii_new)) %>% select(-ii_new)
sl_ii <- new_kids %>% filter(dataset_name == "WF2013")  %>% transmute(subject_id, ii)
hb_ii <- new_kids %>% filter(dataset_name == "Bang2025") %>% transmute(subject_id, ii)
cat(sprintf("New children: %d ELENA-new(EN) + %d SLENA(ES) + %d HABLA(ES) = %d; I %d -> %d; S=%d\n",
            length(elena_new_ids), nrow(slena_kids), nrow(habla_kids), nrow(new_kids), I0, I, S))

## ===== 5. Append Spanish item-level CDI responses (new admins) =====
sl <- sl %>% left_join(sl_ii, by = "subject_id") %>% mutate(admin_key = paste(ii, age, form, sep = "_"))
adm_sl <- sl %>% distinct(admin_key, ii, age) %>% mutate(aa = A0 + row_number())
sl <- sl %>% left_join(adm_sl %>% select(admin_key, aa), by = "admin_key")
A <- A0 + nrow(adm_sl)
y_all  <- c(sd$y,  sl$produces)
aa_all <- c(sd$aa, sl$aa)
jj_all <- c(sd$jj, sl$jj)
admin_to_child <- c(sd$admin_to_child, adm_sl$ii[order(adm_sl$aa)])
admin_age      <- c(sd$admin_age,      adm_sl$age[order(adm_sl$aa)])
cat(sprintf("Spanish CDI: +%d admins, +%d responses; A %d -> %d, N %d -> %d\n",
            nrow(adm_sl), nrow(sl), A0, A, sd$N, length(y_all)))

## ===== 6. Spanish LWL (SLENA n29 + HABLA DRT) =====
slwl <- read_csv(RAW("WF2013/SLENA_LWL_LENA_n29.csv"), show_col_types = FALSE) %>%
  transmute(subject_id = as.character(SubNum), r18 = rtmsec18known_3001800, r24 = RT24_D300) %>%
  pivot_longer(c(r18, r24), names_to = "tp", values_to = "rt") %>%
  transmute(subject_id, age = ifelse(tp == "r18", 18, 24), rt) %>% filter(!is.na(rt), rt > 0) %>%
  inner_join(sl_ii, by = "subject_id")
hlwl <- hb %>% transmute(subject_id = ID, r18 = DRT18mKnown, r21 = DRT21mKnown, r25 = DRT25mKnown) %>%
  pivot_longer(c(r18, r21, r25), names_to = "tp", values_to = "rt") %>%
  transmute(subject_id, age = as.numeric(sub("r", "", tp)), rt) %>% filter(!is.na(rt), rt > 0) %>%
  inner_join(hb_ii, by = "subject_id")
sp_lwl <- bind_rows(slwl, hlwl)
lwl_to_child_all <- c(sd$lwl_to_child, sp_lwl$ii)
lwl_log_age_all  <- c(sd$lwl_log_age,  log(sp_lwl$age / a0))
lwl_log_rt_all   <- c(sd$lwl_log_rt,   log(sp_lwl$rt))
cat(sprintf("Spanish LWL: SLENA %d + HABLA %d obs; N_lwl %d -> %d\n",
            nrow(slwl), nrow(hlwl), sd$N_lwl, length(lwl_log_rt_all)))

## ===== 7. Spanish input (SLENA LENA @18 + HABLA AWC @18,25) =====
sin <- read_csv(RAW("WF2013/SLENA_LWL_LENA_n29.csv"), show_col_types = FALSE) %>%
  transmute(subject_id = as.character(SubNum), log_input = log(FreqAWChr)) %>%
  filter(is.finite(log_input)) %>% inner_join(sl_ii, by = "subject_id")
hin <- hb %>% transmute(subject_id = ID, a18 = AD18AWChr, a25 = AD25AWChr) %>%
  pivot_longer(c(a18, a25), names_to = "tp", values_to = "awc") %>%
  transmute(subject_id, log_input = log(awc)) %>% filter(is.finite(log_input)) %>%
  inner_join(hb_ii, by = "subject_id")
sp_in <- bind_rows(sin, hin)
rec_to_child  <- c(sd$rec_to_child,  sp_in$ii)
log_input_obs <- c(sd$log_input_obs, sp_in$log_input)
instr_of_rec  <- c(sd$instr_of_rec,  rep(2L, nrow(sp_in)))      # LENA = instr 2
V_obs <- length(rec_to_child)
## per-study input level prior (extend to S): Spanish studies get their own observed means
mu_r_s_prior_mean <- c(sd$mu_r_s_prior_mean, rep(mean(log_input_obs), 2))
mu_r_s_prior_mean[WF_STUDY] <- mean(sin$log_input); mu_r_s_prior_mean[HB_STUDY] <- mean(hin$log_input)
cat(sprintf("Spanish input: SLENA %d + HABLA %d recs; V_obs %d -> %d\n",
            nrow(sin), nrow(hin), sd$V_obs, V_obs))

## ===== 8. Sumscore count branch (ELENA-WS EN + HABLA-WG/WS ES) =====
## Forms: 1 = English-WS (all English items), 2 = Spanish-WG, 3 = Spanish-WS.
f1 <- seq_len(J0)                                              # English bank
f2 <- sp_items$jj[sp_items$def %in% (mwg$def)]                 # Spanish WG-form items
f3 <- sp_items$jj[sp_items$def %in% (mws$def)]                 # Spanish WS-form items
form_item  <- c(f1, f2, f3)
form_len   <- c(length(f1), length(f2), length(f3))
form_start <- c(1L, 1L + cumsum(form_len)[-3])
n_forms <- 3L
## admins
el_adm <- el %>% transmute(ii, age, k, form = 1L)
hwg <- hb %>% transmute(subject_id = ID, ii = hb_ii$ii[match(ID, hb_ii$subject_id)],
                        age = as.numeric(CDIAge18m), k = as.numeric(WordsProd18m), form = 2L)
hws <- bind_rows(
  hb %>% transmute(ID, age = as.numeric(CDIWSAge),     k = as.numeric(CDIVocPost21)),
  hb %>% transmute(ID, age = as.numeric(CDIAgePost25), k = as.numeric(CDIVocPost25))) %>%
  transmute(ii = hb_ii$ii[match(ID, hb_ii$subject_id)], age, k, form = 3L)
sum_adm <- bind_rows(el_adm, hwg, hws) %>% filter(!is.na(ii), !is.na(age), !is.na(k), k >= 0) %>%
  mutate(k = pmin(k, form_len[form]))                          # guard k <= form size
n_sum <- nrow(sum_adm)
cat(sprintf("Count branch: %d admins (ELENA-WS %d, HABLA-WG %d, HABLA-WS %d); forms n=%s\n",
            n_sum, sum(sum_adm$form == 1), sum(sum_adm$form == 2), sum(sum_adm$form == 3),
            paste(form_len, collapse = "/")))

## ===== 9. Assemble stan_data =====
cc_all <- c(sd$cc, rep(SP_CC, J_es))
stan_data <- modifyList(sd, list(
  N = length(y_all), A = A, I = I, J = J0 + J_es, C = SP_CC, S = S,
  y = y_all, aa = aa_all, jj = jj_all, cc = cc_all,
  admin_to_child = admin_to_child, admin_age = admin_age, study_of_child = study_of_child,
  N_lwl = length(lwl_log_rt_all), lwl_to_child = lwl_to_child_all,
  lwl_log_age = lwl_log_age_all, lwl_log_rt = lwl_log_rt_all,
  V_obs = V_obs, rec_to_child = rec_to_child, log_input_obs = log_input_obs, instr_of_rec = instr_of_rec,
  mu_r_s_prior_mean = mu_r_s_prior_mean,
  n_sum = n_sum, sum_child = sum_adm$ii, sum_log_age = log(sum_adm$age / a0),
  sum_k = as.integer(sum_adm$k), sum_form = as.integer(sum_adm$form),
  n_forms = n_forms, n_form_items = length(form_item),
  form_item = as.integer(form_item), form_start = as.integer(form_start), form_len = as.integer(form_len)
))

## word_info + child_info (aux, extended)
word_info  <- bind_rows(b$word_info, sp_items %>% transmute(jj, item = def, prob = NA_real_, cc = SP_CC))
child_info <- bind_rows(ci %>% mutate(subject_id = as.character(subject_id)),
                        new_kids %>% transmute(ii, lab_subject_id = subject_id,
                                               dataset_name, study, subject_id)) %>% arrange(ii)

bundle <- list(stan_data = stan_data, word_info = word_info, child_info = child_info,
               lwl = b$lwl, datasets = c(b$datasets, "WF2013", "Bang2025"),
               input_recs = tibble(ii = rec_to_child, log_input = log_input_obs, instr = instr_of_rec),
               language = "Bilingual (English + Spanish; bi-lean: item-level + sumscore count branch)")
out <- here("fits/joint_io_proc_bilingual_subset_data.rds")
saveRDS(bundle, out)
cat(sprintf("\n== bilingual bundle ==\n  I=%d A=%d J=%d C=%d S=%d N=%d N_lwl=%d V_obs=%d n_sum=%d\n  Saved %s\n",
            I, A, J0 + J_es, SP_CC, S, length(y_all), length(lwl_log_rt_all), V_obs, n_sum, out))
