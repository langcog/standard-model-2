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
iob <- readRDS(here("fits", "io_pooled_subset_data.rds"))
pb  <- readRDS(here("fits", "proc_dp_all_subset_data.rds"))
io_adm <- iob$df |> distinct(study, ckey, age)
pr_adm <- pb$df  |> distinct(dataset_name, id = as.character(lab_subject_id), age)
n_rec  <- iob$recordings |> count(study)
n_lwl  <- pb$lwl |> count(dataset_name)

## AM2018/FMW2013: union of the io + proc subsamples, keyed on the original
## lab_subject_id (the pooled io ckey is an internal index; recover ids from
## the per-dataset io bundles).
io_ids <- function(bundle, study) {
  readRDS(here("fits", bundle))$df |>
    distinct(id = as.character(subject_id), age) |> mutate(study = study)
}
union_stats <- function(io_df, proc_ds) {
  u <- bind_rows(io_df |> select(id, age),
                 pr_adm |> filter(dataset_name == proc_ds) |> select(id, age)) |> distinct()
  tibble::tibble(n = n_distinct(u$id), n_admins = nrow(u),
                 mean_age = mean(u$age), min_age = min(u$age), max_age = max(u$age))
}
io_only_stats <- function(study) {
  a <- io_adm |> filter(study == !!study)
  tibble::tibble(n = n_distinct(a$ckey), n_admins = nrow(a),
                 mean_age = mean(a$age), min_age = min(a$age), max_age = max(a$age))
}
IOPROC <- tibble::tribble(
  ~citation,                        ~language,            ~input,     ~rec_study, ~lwl_ds,
  "Long et al. (2024)",             "English (American)", "headcam",  "BabyView",  NA,
  "Egan-Dailey & Bergelson (2025)", "English (American)", "LENA",     "SEEDLingS", NA,
  "Adams et al. (2018)",            "English (American)", "LENA",     "AM2018",    "adams_marchman_2018",
  "Fernald et al. (2013)",          "English (American)", "LENA",     "FMW2013",   "fmw_2013",
  "Fernald & Marchman (2012)",      "English (American)", "",         NA,          "fernald_marchman_2012"
)
io_rows <- bind_rows(
  io_only_stats("BabyView"), io_only_stats("SEEDLingS"),
  union_stats(io_ids("io_am2018_subset_data.rds",  "AM2018"),  "adams_marchman_2018"),
  union_stats(io_ids("io_fmw2013_subset_data.rds", "FMW2013"), "fmw_2013"),
  pr_adm |> filter(dataset_name == "fernald_marchman_2012") |>
    summarise(n = n_distinct(id), n_admins = n(), mean_age = mean(age),
              min_age = min(age), max_age = max(age))
) |>
  bind_cols(IOPROC |> select(citation, language, input, rec_study, lwl_ds)) |>
  left_join(n_rec, by = c(rec_study = "study")) |> rename(n_recordings = n.y, n = n.x) |>
  left_join(n_lwl, by = c(lwl_ds = "dataset_name")) |> rename(n_lwl = n.y, n = n.x) |>
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
