## Build the bundle for the D'0-D'3 processing regression ladder
## (log_irt_long_proc_dp.stan).
##
## RUN LOCALLY. Combines, for AM2018 + FM2012 + FMW2013:
##   - item-level CDI (stanford_cdi_items_long.csv; WG+WS, <=30mo)
##   - LWL RT (d_sub joined to 2026.1 lab ids; <=30mo, n_trials_rt>=5, winsorized)
##   - observed LENA input (TL3TLO_LENA.csv; AM2018 + FMW2013), standardized,
##     with sigma_lena pinned empirically from the 16/18mo replicates.
##
## Usage: Rscript model/scripts/prepare_proc_dp_bundle.R [n_items] [datasets]
##   datasets: comma-sep peekbank names (default all 3); e.g. adams_marchman_2018
##
## Output: fits/proc_dp_<tag>_subset_data.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

args      <- commandArgs(trailingOnly = TRUE)
n_items   <- as.integer(if (length(args) >= 1) args[1] else 200)
DATASETS  <- if (length(args) >= 2) {
  strsplit(args[2], ",")[[1]]
} else {
  c("adams_marchman_2018", "fernald_marchman_2012", "fmw_2013", "fernald_totlot")
}
TAG       <- if (length(DATASETS) == 4) "all" else paste(gsub("_.*","",DATASETS), collapse="_")

PB_DIR <- file.path(PROJECT_ROOT, "data/peekbank")
SEED   <- 20260609
N_DIFF_BINS <- 4
AGE_CAP <- 30

# lab study <-> peekbank dataset, and study index
STUDY_MAP <- c(totlot3 = "adams_marchman_2018",
               totlot2 = "fernald_marchman_2012",
               tlo     = "fmw_2013",
               totlot  = "fernald_totlot")     # original TotLot: RT + CDI, no LENA input
DS_LEVELS <- c("adams_marchman_2018", "fernald_marchman_2012", "fmw_2013", "fernald_totlot")
DS_LEVELS <- DS_LEVELS[DS_LEVELS %in% DATASETS]

# ---- 1. Item-level CDI (studies x WG/WS, <=30mo). fernald_totlot CDI comes
#         from the separately-parsed TL_{18,21,25}m_WS.xlsx (parse_totlot_cdi.R). ---- #
.read_cdi <- function(f) read_csv(file.path(PB_DIR, f), show_col_types = FALSE, progress = FALSE) %>%
  transmute(lab_subject_id = as.character(lab_subject_id), study = as.character(study),
            age = as.integer(age), form = as.character(form), item = as.character(item),
            produces = as.integer(produces))
cdi <- bind_rows(.read_cdi("stanford_cdi_items_long.csv"), .read_cdi("totlot_cdi_items_long.csv")) %>%
  mutate(dataset_name = STUDY_MAP[study]) %>%
  filter(dataset_name %in% DATASETS, age <= AGE_CAP, !is.na(produces))
cat(sprintf("CDI rows: %d (kids %d, admins %d) across %s\n",
            nrow(cdi), n_distinct(cdi$lab_subject_id),
            n_distinct(paste(cdi$lab_subject_id, cdi$age)),
            paste(unique(cdi$dataset_name), collapse=",")))

# ---- 2. CHILDES p_j + lexical class (io_pooled = canonical English freq;
#         long_items.rds lost its prob column) ---- #
iop <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))$df %>%
  distinct(item, prob, lexical_category) %>% filter(!is.na(prob), prob > 0)
prob_lookup  <- iop %>% distinct(item, prob)
class_lookup <- iop %>% distinct(item, lexical_category) %>%
  group_by(item) %>% slice(1) %>% ungroup()
cdi <- cdi %>%
  left_join(prob_lookup, by = "item") %>%
  left_join(class_lookup, by = "item") %>%
  filter(!is.na(prob), !is.na(lexical_category))
cat(sprintf("After CHILDES match: %d rows, %d items\n", nrow(cdi), n_distinct(cdi$item)))

# ---- 3. Stratified item subsample (class x difficulty) ---- #
set.seed(SEED)
item_summary <- cdi %>% distinct(item, lexical_category, prob) %>%
  mutate(log_p = log(prob)) %>%
  group_by(lexical_category) %>% mutate(diff_bin = ntile(log_p, N_DIFF_BINS)) %>% ungroup()
n_classes <- n_distinct(item_summary$lexical_category)
per_cell  <- max(1, floor(n_items / (n_classes * N_DIFF_BINS)))
chosen <- item_summary %>% group_by(lexical_category, diff_bin) %>%
  slice_sample(n = per_cell) %>% ungroup() %>% pull(item)
if (length(chosen) > n_items) chosen <- sample(chosen, n_items)
if (length(chosen) < n_items) {
  extras <- setdiff(item_summary$item, chosen)
  chosen <- c(chosen, sample(extras, min(n_items - length(chosen), length(extras))))
}
cdi <- cdi %>% filter(item %in% chosen)
cat(sprintf("Items kept: %d\n", n_distinct(cdi$item)))

# ---- 4. LWL RT: d_sub joined to 2026.1 lab ids, QC'd ---- #
dsub <- readRDS(file.path(PB_DIR, "1_d_sub.Rds")) %>% ungroup()
amp  <- readRDS(file.path(PB_DIR, "_pb2026_admins.rds")) %>%
  mutate(lab_subject_id = as.character(lab_subject_id))
lwl <- dsub %>% filter(dataset_name %in% DATASETS, !is.na(log_rt)) %>%
  inner_join(amp %>% select(administration_id, lab_subject_id), by = "administration_id") %>%
  filter(age <= AGE_CAP, n_trials_rt >= 5) %>%
  group_by(dataset_name) %>%
  mutate(log_rt_w = pmin(log_rt, quantile(log_rt, 0.99, na.rm = TRUE))) %>%
  ungroup() %>%
  select(dataset_name, lab_subject_id, lwl_age = age, lwl_log_rt = log_rt_w)
cat(sprintf("LWL admins (QC'd): %d (kids %d)\n", nrow(lwl), n_distinct(lwl$lab_subject_id)))

# ---- 5. Child set = CDI ∩ RT; build indices ---- #
kids <- intersect(unique(cdi$lab_subject_id), unique(lwl$lab_subject_id))
cdi <- cdi %>% filter(lab_subject_id %in% kids)
lwl <- lwl %>% filter(lab_subject_id %in% kids)

cdi <- cdi %>%
  mutate(admin_key = paste(lab_subject_id, age, sep = "_"),
         aa = as.integer(factor(admin_key)),
         ii = as.integer(factor(lab_subject_id, levels = kids)),
         jj = as.integer(factor(item)),
         cc = as.integer(factor(lexical_category)))
admin_info <- cdi %>% distinct(aa, ii, lab_subject_id, age, dataset_name) %>% arrange(aa)
word_info  <- cdi %>% group_by(jj) %>%
  summarise(item = first(item), prob = first(prob), cc = first(cc), .groups = "drop") %>% arrange(jj)
child_info <- cdi %>% distinct(ii, lab_subject_id, dataset_name) %>% arrange(ii) %>%
  mutate(study = as.integer(factor(dataset_name, levels = DS_LEVELS)))

I <- max(cdi$ii); A <- max(cdi$aa); J <- max(cdi$jj); C <- max(cdi$cc); S <- length(DS_LEVELS)
a0_dataset <- round(median(admin_info$age))

# LWL -> child ii
lwl <- lwl %>% inner_join(child_info %>% select(lab_subject_id, ii), by = "lab_subject_id")
N_lwl <- nrow(lwl)
cat(sprintf("Bundle: I=%d A=%d J=%d C=%d S=%d N=%d N_lwl=%d a0=%d\n",
            I, A, J, C, S, nrow(cdi), N_lwl, a0_dataset))

# ---- 6. Observed LENA input + empirical sigma_lena ---- #
LENA_STUDY <- c(TL3 = "adams_marchman_2018", TLO = "fmw_2013")
lena_raw <- read_csv(file.path(PB_DIR, "TL3TLO_LENA.csv"), show_col_types = FALSE, progress = FALSE) %>%
  mutate(lab_subject_id = as.character(SubjectID1), dataset_name = LENA_STUDY[Study]) %>%
  filter(dataset_name %in% DATASETS, lab_subject_id %in% kids)
# log AWC at 16 + 18 mo (replicates)
lena <- lena_raw %>%
  transmute(lab_subject_id, dataset_name,
            l16 = log(AWCHr16M), l18 = log(AWCHr18M)) %>%
  rowwise() %>% mutate(lmean = mean(c(l16, l18), na.rm = TRUE)) %>% ungroup() %>%
  filter(is.finite(lmean))
# reliability from the 16-18mo replicate (kids with both)
both <- lena %>% filter(is.finite(l16), is.finite(l18))
sigma_w2  <- mean((both$l16 - both$l18)^2, na.rm = TRUE) / 2          # single-read noise var
# standardize within study by the noise-corrected between-child SD
lena <- lena %>% group_by(dataset_name) %>%
  mutate(sd_obs = sd(lmean, na.rm = TRUE),
         sd_true = sqrt(pmax(sd_obs^2 - sigma_w2/2, 1e-4)),
         z_lena = (lmean - mean(lmean, na.rm = TRUE)) / sd_true) %>% ungroup()
sigma_lena <- sqrt(sigma_w2/2) / mean(lena$sd_true)                   # noise of the child-mean read, std units
lena <- lena %>% inner_join(child_info %>% select(lab_subject_id, ii), by = "lab_subject_id")
V <- nrow(lena)
cat(sprintf("LENA: %d kids with input, sigma_lena=%.3f (reliability=%.2f)\n",
            V, sigma_lena, 1 - sigma_lena^2/(1 + sigma_lena^2)))

# ---- 7. stan_data ---- #
prior_r <- load_input_rate_prior()
stan_data <- c(
  list(
    N = nrow(cdi), A = A, I = I, J = J, C = C, S = S,
    aa = cdi$aa, jj = cdi$jj,
    admin_to_child = admin_info$ii, cc = word_info$cc, y = cdi$produces,
    study_of_child = child_info$study,
    admin_age = admin_info$age, log_p = log(word_info$prob),
    log_H = MODEL_CONSTANTS$log_H, a0 = a0_dataset,
    mu_r = prior_r$mu_r, sigma_r = prior_r$sigma_r,
    # LWL
    N_lwl = N_lwl, lwl_to_child = lwl$ii,
    lwl_log_age = log(lwl$lwl_age / a0_dataset), lwl_log_rt = lwl$lwl_log_rt,
    mu_rt_prior_mean = 6.5, mu_rt_prior_sd = 1,
    mu_rtslope_prior_sd = 1, sigma_rt0_prior_sd = 1, sigma_rt1_prior_sd = 1,
    sigma_lwl_prior_sd = 1,
    # LENA
    V = V, rec_to_child = lena$ii, z_lena = lena$z_lena, sigma_lena = sigma_lena,
    # ladder defaults = D'0 (gamma_in on, betas pinned near 0); fit driver overrides
    beta_xi_prior_sd = 0.001,
    beta_k0_prior_sd = 0.001, beta_k1_prior_sd = 0.001
  ),
  # slopes free (longitudinal); gamma_in on (D'0 baseline). gamma_in_prior_sd
  # already lives in DEFAULT_PRIORS, so set it here to avoid a duplicate name.
  modifyList(DEFAULT_PRIORS, list(sigma_zeta_prior_sd = 1, gamma_in_prior_sd = 1))
)

bundle <- list(stan_data = stan_data, admin_info = admin_info, word_info = word_info,
               child_info = child_info, lwl = lwl, lena = lena, df = cdi,
               datasets = DS_LEVELS,
               language = sprintf("English (proc_dp: %s)", paste(DS_LEVELS, collapse=", ")))
out <- file.path(PATHS$fits_dir, sprintf("proc_dp_%s_subset_data.rds", TAG))
saveRDS(bundle, out)
cat(sprintf("\nSaved %s\n", out))
