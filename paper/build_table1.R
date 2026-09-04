## build_table1.R — dataset-characteristics table (Table 1).
## Writes paper/cache/table1_datasets.csv (COMMITTED) with one row per
## dataset used anywhere in the paper, in three groups:
##   * Wordbank longitudinal — the five by-study glmer/IRT datasets
##     (fits/glmer_ladder/data_<study>.rds; EN D pools Thal+Smith+Marchman)
##   * Input / processing — the io_pooled + proc_dp datasets; for AM2018 and
##     FMW2013 (in both analyses, slightly different QC subsets) N / admins /
##     ages are the UNION of the two subsamples
##   * Wordbank cross-sectional — the 31 languages in the fig-demographics
##     cross-sectional arm (studies/cross_sectional_demographics/00_build.R; one
##     admin per child, full archive N per language as of the 2026-06-11
##     uncapped refit — no per-language subsample cap)
##
## Stats are computed from the analysis bundles/frames, so the table shows
## the data the models actually saw (post-QC), not raw archive sizes.
## Run: Rscript paper/build_table1.R   (needs local frames + fit bundles)

suppressPackageStartupMessages({ library(here); library(dplyr) })

## ---- Wordbank longitudinal (the five by-study datasets in the paper) ----
LONG <- tibble::tribble(
  ~slug,        ~citation,                ~language,
  "thal",       "Thal et al. (20XX)",     "English (American)",
  "smith",      "Smith et al. (20XX)",    "English (American)",
  "marchman",   "Marchman et al. (20XX)", "English (American)",
  "norwegian",  "Simonsen et al. (2014)", "Norwegian",
  "japanese",   "Hagihara et al. (2023)", "Japanese"
)
long_rows <- bind_rows(lapply(seq_len(nrow(LONG)), function(i) {
  b  <- readRDS(here("fits", "glmer_ladder", paste0("data_", LONG$slug[i], ".rds")))
  ad <- b$df |> distinct(child_id, age)
  tibble::tibble(group = "Wordbank longitudinal",
                 citation = LONG$citation[i], language = LONG$language[i],
                 n = b$n_kids, n_admins = nrow(ad),
                 mean_age = mean(ad$age), min_age = min(ad$age), max_age = max(ad$age),
                 longitudinal = "yes", input = "", n_recordings = NA, n_lwl = NA)
}))

## ---- Input / processing datasets ----
## 2026-09: read straight from the bundle the wired io-proc fit trained on
## (joint_io_proc_lean_d2_enct_2k <- joint_io_proc_english_count bundle), so
## Table 1 is the analysis sample by construction. Replaces the old stitch of
## the io_pooled + proc_dp bundles, which predates fernald_totlot, the
## SEEDLingS RT wiring, the all-items switch and the ELENA count kids.
## Study index -> dataset (b$datasets order); the 95 ELENA-WS sumscore
## admins belong to fmw_2013 and are counted as administrations there.
jb <- readRDS(here("fits", "joint_io_proc_english_count_subset_data.rds")); jsd <- jb$stan_data
ci <- jb$child_info
adm <- bind_rows(
  tibble::tibble(ii = jsd$admin_to_child, age = jsd$admin_age),                       # item-level CDI admins
  tibble::tibble(ii = jsd$sum_child,      age = exp(jsd$sum_log_age) * jsd$a0)) |>    # sumscore (count) admins
  left_join(ci |> select(ii, study), by = "ii")
n_rec <- tibble::tibble(ii = jsd$rec_to_child) |> left_join(ci |> select(ii, study), by = "ii") |> count(study, name = "n_recordings")
n_lwl <- tibble::tibble(ii = jsd$lwl_to_child) |> left_join(ci |> select(ii, study), by = "ii") |> count(study, name = "n_lwl")
IOPROC <- tibble::tribble(
  ~study, ~citation,                           ~language,            ~input,
  5L,     "Long et al. (2024)",                "English (American)", "headcam",
  6L,     "Egan-Dailey & Bergelson (2025)",    "English (American)", "LENA",
  1L,     "Adams et al. (2018)",               "English (American)", "LENA",
  3L,     "Fernald et al. (2013)",             "English (American)", "LENA",
  2L,     "Fernald & Marchman (2012)",         "English (American)", "",
  4L,     "Fernald, Perfors & Marchman (2006)","English (American)", ""
)
io_rows <- IOPROC |>
  rowwise() |>
  mutate(n = sum(ci$study == study),
         n_admins = sum(adm$study == study),
         mean_age = mean(adm$age[adm$study == study]),
         min_age  = min(adm$age[adm$study == study]),
         max_age  = max(adm$age[adm$study == study])) |>
  ungroup() |>
  left_join(n_rec, by = "study") |>
  left_join(n_lwl, by = "study") |>
  mutate(group = "Input / processing", longitudinal = "yes") |>
  select(group, citation, language, n, n_admins, mean_age, min_age, max_age,
         longitudinal, input, n_recordings, n_lwl)

## ---- Wordbank cross-sectional (fig-demographics arm) ----
xs_rows <- bind_rows(lapply(
  sort(list.files(here("studies", "cross_sectional_demographics", "cache", "frames"), full.names = TRUE)),
  function(f) {
    fr <- readRDS(f); ch <- fr |> distinct(child_id, age)
    tibble::tibble(group = "Wordbank cross-sectional",
                   citation = "Frank et al. (2017)", language = fr$language[1],
                   n = nrow(ch), n_admins = nrow(ch),
                   mean_age = mean(ch$age), min_age = min(ch$age), max_age = max(ch$age),
                   longitudinal = "", input = "", n_recordings = NA, n_lwl = NA)
  })) |> arrange(language)

tab <- bind_rows(long_rows, io_rows, xs_rows) |> mutate(mean_age = round(mean_age, 1))
out <- here("paper", "cache", "table1_datasets.csv")
write.csv(tab, out, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows: %s)\n", out, nrow(tab),
            paste(table(factor(tab$group, levels = unique(tab$group))), collapse = "/")))
print(as.data.frame(tab |> filter(group != "Wordbank cross-sectional")), row.names = FALSE)
