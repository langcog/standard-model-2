## Data figure for the io-proc analysis: rows = the three data streams, columns
## = the six datasets (one child = one line; thick line = per-dataset loess).
##   A. CDI vocabulary (all six datasets)
##   B. LWL reaction time (five; BabyView has no LWL task)
##   C. Observed input (four; FM2012 / FPM2006 have no recordings). BabyView is
##      averaged to monthly means so its many head-cam videos don't dominate;
##      input is in study-specific units, so only the age trends are comparable
##      across columns, not the levels.
## Colour = data stream (vocabulary neutral, processing, input), matching the
## graphical model and the io-proc figure; empty cells mark channels a dataset
## did not collect. Defines io_panels_fig() (returns the patchwork; the paper's
## fig_data.R passes its palette/theme and adds the panel letters); writes
## figs/data_checks/io_panels.png when run directly. RUN LOCALLY.
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})

io_panels_fig <- function(pal = c(vocab = "#404040", processing = "#AA3377", input = "#009E73"),
                          base_theme = theme_bw(base_size = 8), txt = 7 / .pt) {
  ## english_count = the bundle the wired fit (joint_io_proc_lean_d2_enct_2k)
  ## was trained on; only b$lwl is read here, but keep figure + fit on one bundle.
  b <- readRDS(here("fits/joint_io_proc_english_count_subset_data.rds"))
  ## column order: datasets with recordings first, LWL-only last, so the empty
  ## cells sit at the edges of rows B and C
  DSORD <- c("BabyView","SEEDLingS","AM2018","FMW2013","FM2012","fernald_totlot")
  DSLAB <- c(BabyView = "BabyView", SEEDLingS = "SEEDLingS", AM2018 = "AM2018",
             FMW2013 = "FMW2013", FM2012 = "FM2012", fernald_totlot = "FPM2006")
  fct  <- function(x) factor(DSLAB[x], levels = DSLAB)
  AGES <- c(5, 37)
  missing_cells <- function(present, y)          # "not collected" placeholder
    tibble(ds = fct(setdiff(DSORD, present)), x = mean(AGES), y = y)
  spaghetti <- function(col, group) list(
    geom_line(aes(group = .data[[group]]), colour = col, alpha = 0.35, linewidth = 0.2),
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = col, linewidth = 0.7, span = 1),
    facet_wrap(~ ds, nrow = 1, drop = FALSE),
    coord_cartesian(xlim = AGES),
    scale_x_continuous(breaks = c(12, 24, 36)),
    base_theme, theme(panel.grid.minor = element_blank()))
  no_strip <- theme(strip.text = element_blank(), strip.background = element_blank())

  ## (A) CDI: per-child vocabulary count vs age + per-dataset smoother
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
    mutate(ds = fct(ds), line = paste(id, form))
  pA <- ggplot(cdi, aes(age, vocab)) + spaghetti(pal[["vocab"]], "line") +
    labs(x = NULL, y = "words produced", title = "CDI vocabulary")

  ## (B) LWL reaction time (session means)
  nm <- c(adams_marchman_2018 = "AM2018", fernald_marchman_2012 = "FM2012", fmw_2013 = "FMW2013",
          fernald_totlot = "fernald_totlot", seedlings_zhu = "SEEDLingS")
  rt <- b$lwl %>% mutate(ds = fct(nm[dataset_name]), gid = paste0(dataset_name, "::", ii)) %>%
    group_by(ds, gid, lwl_age) %>% summarise(rt = mean(exp(lwl_log_rt)), .groups = "drop")
  pB <- ggplot(rt, aes(lwl_age, rt)) + spaghetti(pal[["processing"]], "gid") +
    geom_text(data = missing_cells(unname(nm), 800), aes(x, y), label = "not collected",
              colour = "grey55", size = txt) +
    scale_y_log10(breaks = c(400, 600, 1000, 1500)) +
    labs(x = NULL, y = "RT (ms, log scale)", title = "Looking-while-listening reaction time") + no_strip

  ## (C) Observed input (BabyView averaged to monthly means)
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
  inp <- bind_rows(lena, fmw, seed, bv) %>% filter(is.finite(log_input), age_mo > 0) %>%
    mutate(child = paste(ds, child), ds = fct(ds))
  pC <- ggplot(inp, aes(age_mo, log_input)) + spaghetti(pal[["input"]], "child") +
    geom_text(data = missing_cells(c("AM2018","FMW2013","SEEDLingS","BabyView"), 7), aes(x, y),
              label = "not collected", colour = "grey55", size = txt) +
    labs(x = "age (months)", y = "log input rate", title = "Observed input") + no_strip

  pA / pB / pC
}

if (sys.nframe() == 0) {   # run directly -> write the standalone artifact
  dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
  ggsave(here("figs", "data_checks", "io_panels.png"), io_panels_fig() + plot_annotation(tag_levels = "A"),
         width = 7, height = 3.6, dpi = 200)
  message("wrote figs/data_checks/io_panels.png")
}
