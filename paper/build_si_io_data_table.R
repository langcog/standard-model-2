## Writes paper/cache/si_io_data_table.csv (COMMITTED): one row per dataset in the
## joint io-proc bundle, with the channels (CDI vocab / LWL RT / observed input)
## each contributes. Mirrors the working data table; the SI .qmd renders this CSV.
## RUN LOCALLY.
suppressPackageStartupMessages({ library(dplyr); library(here) })

b  <- readRDS(here("fits", "joint_io_proc_mm_subset_data.rds")); sd <- b$stan_data
ci <- b$child_info
adm <- tibble(ii = sd$admin_to_child, age = sd$admin_age)
rt_ii  <- unique(sd$lwl_to_child)
inp_ii <- unique(sd$rec_to_child)

## study index -> dataset_name + publication (see data/peekbank/README.md key)
META <- tibble::tribble(
  ~study, ~dataset,                 ~publication,
  1L, "adams_marchman_2018",   "Adams, Marchman et al. (2018)",
  2L, "fernald_marchman_2012", "Fernald \\& Marchman (2012)",
  3L, "fmw_2013",              "Fernald, Marchman \\& Weisleder (2013)",
  4L, "fernald_totlot",        "Fernald, Perfors \\& Marchman (2006)",
  5L, "BabyView",              "BabyView (Long et al.)",
  6L, "SEEDLingS",             "SEEDLingS; LWL from Zhu et al. (submitted)")

tab <- META %>% rowwise() %>% mutate(
  kids       = sum(ci$study == study),
  cdi_admins = sum(adm$ii %in% ci$ii[ci$study == study]),
  rt_obs     = sum(sd$lwl_to_child %in% ci$ii[ci$study == study]),
  input_kids = length(intersect(inp_ii, ci$ii[ci$study == study])),
  age_lo     = min(adm$age[adm$ii %in% ci$ii[ci$study == study]]),
  age_hi     = max(adm$age[adm$ii %in% ci$ii[ci$study == study]]),
  ages       = sprintf("%d--%d", age_lo, age_hi)
) %>% ungroup() %>%
  mutate(vocab = "Y", rt = ifelse(rt_obs > 0, "Y", "--"), input = ifelse(input_kids > 0, "Y", "--")) %>%
  select(dataset, publication, kids, cdi_admins, rt_obs, input_kids, ages, vocab, rt, input)

total <- tibble(dataset = "Total", publication = "",
  kids = sum(tab$kids), cdi_admins = sum(tab$cdi_admins), rt_obs = sum(tab$rt_obs),
  input_kids = sum(tab$input_kids), ages = "", vocab = "", rt = "", input = "")
tab <- bind_rows(tab, total)

out <- here("paper", "cache", "si_io_data_table.csv")
write.csv(tab, out, row.names = FALSE)
cat(sprintf("Wrote %s (%d datasets; both-channel input&RT = %d kids)\n",
            out, nrow(tab) - 1L, length(intersect(rt_ii, inp_ii))))
print(as.data.frame(tab), row.names = FALSE)
