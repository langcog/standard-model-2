## Combined io-proc data-check figure: 3 panels, dataset = colour (constant across panels).
##   A. CDI vocabulary spaghetti (per-child lines, no points) + per-dataset smoother
##   B. LWL reaction-time spaghetti + per-dataset trend
##   C. Observed input spaghetti + per-dataset trend (BabyView averaged to monthly means
##      so its many head-cam videos don't swamp the LENA datasets; input is in study-
##      specific units, so cross-dataset levels are not comparable -- only the trends).
## Defines io_panels_fig() (returns the patchwork) so the SI can source + render it;
## writes figs/data_checks/io_panels.png when run directly. RUN LOCALLY.
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})

io_panels_fig <- function() {
  b <- readRDS(here("fits/joint_io_proc_mm_subset_data.rds"))
  DSORD <- c("BabyView","SEEDLingS","AM2018","FM2012","FMW2013","fernald_totlot")
  PAL <- setNames(c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#a6761d"), DSORD)
  fct <- function(x) factor(x, levels = DSORD)

  ## (1) CDI: per-child vocabulary count vs age (no points) + per-dataset smoother
  stanford <- read_csv(here("data/intermediates/stanford_cdi_items_long.csv"), show_col_types = FALSE) %>%
    transmute(ds = recode(study, totlot2 = "FM2012", totlot3 = "AM2018", tlo = "FMW2013"),
              id = as.character(lab_subject_id), age, form, item, produces)
  totlot <- read_csv(here("data/intermediates/totlot_cdi_items_long.csv"), show_col_types = FALSE) %>%
    transmute(ds = "fernald_totlot", id = as.character(lab_subject_id), age, form = "WS", item, produces)
  iop <- readRDS(here("fits/io_pooled_subset_data.rds"))$df %>% filter(study %in% c("BabyView","SEEDLingS")) %>%
    transmute(ds = study, id = as.character(ckey), age, form, item, produces)
  cdi <- bind_rows(stanford, totlot, iop) %>%
    group_by(ds, id, age, form, item) %>% summarise(produces = max(produces, na.rm = TRUE), .groups = "drop") %>%
    group_by(ds, id, age, form) %>% summarise(vocab = sum(produces == 1, na.rm = TRUE), .groups = "drop") %>%
    mutate(ds = fct(ds))
  p1 <- ggplot(cdi, aes(age, vocab, color = ds)) +
    geom_line(aes(group = interaction(id, form)), alpha = 0.28, linewidth = 0.25) +
    geom_smooth(aes(group = ds), method = "loess", se = FALSE, linewidth = 1, span = 1) +
    scale_color_manual(values = PAL, name = "dataset", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.4, alpha = 1), nrow = 1)) +
    labs(x = "age (months)", y = "words produced", title = "A. CDI vocabulary") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

  ## (2) LWL reaction time
  nm <- c(adams_marchman_2018 = "AM2018", fernald_marchman_2012 = "FM2012", fmw_2013 = "FMW2013",
          fernald_totlot = "fernald_totlot", seedlings_zhu = "SEEDLingS")
  rt <- b$lwl %>% mutate(ds = fct(nm[dataset_name]), gid = paste0(dataset_name, "::", ii)) %>%
    group_by(ds, gid, lwl_age) %>% summarise(rt = mean(exp(lwl_log_rt)), .groups = "drop")
  p2 <- ggplot(rt, aes(lwl_age, rt, color = ds)) +
    geom_line(aes(group = gid), alpha = 0.3, linewidth = 0.25) +
    geom_smooth(aes(group = ds), method = "loess", se = FALSE, linewidth = 1, span = 1) +
    scale_y_log10(breaks = c(400, 600, 1000, 1500)) +
    scale_color_manual(values = PAL, drop = FALSE, guide = "none") +
    labs(x = "age (months)", y = "RT (ms, log scale)", title = "B. LWL reaction time") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

  ## (3) Observed input (BabyView averaged to monthly means)
  lena <- read_csv(here("data/raw/AM2018/lena_am2018_fmw2013.csv"), show_col_types = FALSE) %>% filter(Study == "TL3") %>%
    transmute(ds = "AM2018", child = as.character(SubjectID1), a16 = AGE16M, r16 = AWCHr16M, a18 = AGE18M, r18 = AWCHr18M) %>%
    pivot_longer(c(a16, r16, a18, r18), names_to = c(".value","tp"), names_pattern = "([ar])(16|18)") %>%
    transmute(ds, child, age_mo = a, log_input = log(r))
  fmw <- read_csv(here("data/raw/FMW2013/TLOELENA_LENA_1824.csv"), show_col_types = FALSE) %>%
    transmute(ds = "FMW2013", child = as.character(SubjectID1), a18 = AGE18M, r18 = AWCHr18M, a24 = AGE24M, r24 = AWCHr24M) %>%
    pivot_longer(c(a18, r18, a24, r24), names_to = c(".value","tp"), names_pattern = "([ar])(18|24)") %>%
    transmute(ds, child, age_mo = a, log_input = log(r))
  seed <- read_csv(here("data/raw/seedlings/lena_data.csv"), show_col_types = FALSE) %>%
    transmute(ds = "SEEDLingS", child = as.character(subj), age_mo = month, log_input = log(awc_perhr)) %>% filter(age_mo <= 30)
  bv <- readRDS(here("fits/babyview_subset_data.rds"))$videos %>%
    transmute(child = as.character(subject_id), age_mo, rate = exp(log_r_obs)) %>%
    mutate(month = round(age_mo)) %>% group_by(child, month) %>%
    summarise(age_mo = mean(age_mo), log_input = log(mean(rate)), .groups = "drop") %>%
    transmute(ds = "BabyView", child, age_mo, log_input)
  inp <- bind_rows(lena, fmw, seed, bv) %>% filter(is.finite(log_input), age_mo > 0) %>% mutate(ds = fct(ds))
  p3 <- ggplot(inp, aes(age_mo, log_input, color = ds)) +
    geom_line(aes(group = child), alpha = 0.3, linewidth = 0.25) +
    geom_smooth(aes(group = ds), method = "loess", se = FALSE, linewidth = 1, span = 1) +
    scale_color_manual(values = PAL, drop = FALSE, guide = "none") +
    labs(x = "age (months)", y = "log input rate (study-specific units)", title = "C. Observed input") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

  (p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

if (sys.nframe() == 0) {   # run directly -> write the standalone artifact
  dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
  ggsave(here("figs", "data_checks", "io_panels.png"), io_panels_fig(), width = 14, height = 5, dpi = 140)
  message("wrote figs/data_checks/io_panels.png")
}
