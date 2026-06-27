## Three by-dataset diagnostic figures (facet_wrap per dataset) over ALL io-proc data,
## including the new Spanish + sumscore cohorts. English item-level comes from the
## harmonized CSVs; Spanish/sumscore cohorts are read from their raw sources.
##   figs/data_checks/io_diag_cdi.png    vocabulary (produced) by age
##   figs/data_checks/io_diag_rt.png     LWL reaction time by age
##   figs/data_checks/io_diag_input.png  input rate by age
## RUN LOCALLY.
suppressPackageStartupMessages({library(here); library(dplyr); library(tidyr); library(readr); library(readxl); library(ggplot2)})
DSORD <- c("babyview","SEEDLingS","AM2018","FM2012","FMW2013","FPM2006","WF2013","Bang2025")
fct <- function(x) factor(x, levels = DSORD)
LANG <- c(WF2013="Spanish", Bang2025="Spanish")   # default English
lang_of <- function(pc) ifelse(pc %in% names(LANG), "Spanish", "English")

## ===== CDI: vocabulary produced per admin =====
eng <- read_csv(here("data/harmonized/cdi_item_level.csv"), show_col_types = FALSE) %>%
  filter(item_canonical) %>% mutate(child_id = as.character(child_id)) %>%
  group_by(paper_code, child_id, age, form) %>% summarise(vocab = sum(produces == 1), .groups = "drop") %>%
  mutate(grain = "item")
# ELENA WS sumscores (FMW2013)
ew <- read_excel(here("data/raw/FMW2013/elena/ELENA_WS_SumScores.xlsx")) %>%
  transmute(paper_code = "FMW2013", child_id = as.character(ParticipantId), age = as.numeric(CDIAge),
            form = "WS", vocab = as.numeric(VOCAB), grain = "sumscore")
# WF2013 SLENA (VOCAB totals, WG + WS)
slena <- bind_rows(
  read_excel(here("data/raw/WF2013/SLENA_WGCDIS_toWordbank.xlsx"), .name_repair="unique_quiet") %>% mutate(form="WG"),
  read_excel(here("data/raw/WF2013/SLENA_CDIS_toWordbank.xlsx"),  .name_repair="unique_quiet") %>% mutate(form="WS")) %>%
  transmute(paper_code = "WF2013", child_id = as.character(ParticipantId), age = as.numeric(CDIAge),
            form, vocab = as.numeric(VOCAB), grain = "item") %>% filter(!is.na(vocab))
# HABLA (Bang2025) production sumscores at 18 (WG) / 21,25 (WS)
hb <- read_csv(here("data/raw/Bang2025/Habla1.0_LENALWLCDISumScores.csv"), show_col_types = FALSE)
hbcdi <- bind_rows(
  hb %>% transmute(child_id=as.character(ID), age=as.numeric(CDIAge18m),  form="WG", vocab=as.numeric(WordsProd18m)),
  hb %>% transmute(child_id=as.character(ID), age=as.numeric(CDIWSAge),    form="WS", vocab=as.numeric(CDIVocPost21)),
  hb %>% transmute(child_id=as.character(ID), age=as.numeric(CDIAgePost25),form="WS", vocab=as.numeric(CDIVocPost25))) %>%
  mutate(paper_code="Bang2025", grain="sumscore") %>% filter(!is.na(vocab), !is.na(age))
cdi <- bind_rows(eng, ew, slena, hbcdi) %>% filter(!is.na(vocab), !is.na(age)) %>%
  mutate(paper_code = fct(paper_code), form = factor(form, c("WG","WS")))

pc <- ggplot(cdi, aes(age, vocab, color = form)) +
  geom_line(aes(group = interaction(child_id, form)), color = "grey75", alpha = 0.4, linewidth = 0.25) +
  geom_point(aes(shape = grain), alpha = 0.6, size = 0.8) +
  facet_wrap(~ paper_code, scales = "free_y") +
  scale_color_manual(values = c(WG = "#1f78b4", WS = "#e31a1c"), name = "form") +
  scale_shape_manual(values = c(item = 16, sumscore = 1), name = "grain") +
  labs(x = "age (months)", y = "words produced",
       title = "CDI vocabulary by dataset (open = sumscore-only; English item-level from harmonized table)") +
  theme_bw(base_size = 9) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(here("figs/data_checks/io_diag_cdi.png"), pc, width = 11, height = 6, dpi = 130)

## ===== RT: LWL reaction time per session =====
rt_eng <- read_csv(here("data/harmonized/lwl_session_level.csv"), show_col_types = FALSE) %>%
  transmute(paper_code, child_id = as.character(child_id), age, rt_ms)
sl_lwl <- read_csv(here("data/raw/WF2013/SLENA_LWL_LENA_n29.csv"), show_col_types = FALSE) %>%
  transmute(child_id = as.character(SubNum), r18 = rtmsec18known_3001800, r24 = RT24_D300) %>%
  pivot_longer(c(r18, r24), names_to = "tp", values_to = "rt_ms") %>%
  transmute(paper_code = "WF2013", child_id, age = ifelse(tp == "r18", 18, 24), rt_ms) %>% filter(!is.na(rt_ms))
hb_rt <- hb %>% transmute(child_id = as.character(ID), r18 = DRT18mKnown, r21 = DRT21mKnown, r25 = DRT25mKnown) %>%
  pivot_longer(c(r18, r21, r25), names_to = "tp", values_to = "rt_ms") %>%
  transmute(paper_code = "Bang2025", child_id, age = as.numeric(sub("r","",tp)), rt_ms) %>% filter(!is.na(rt_ms))
rt <- bind_rows(rt_eng, sl_lwl, hb_rt) %>% filter(!is.na(rt_ms), rt_ms > 0) %>% mutate(paper_code = fct(paper_code))
pr <- ggplot(rt, aes(age, rt_ms)) +
  geom_line(aes(group = child_id), color = "grey75", alpha = 0.3, linewidth = 0.25) +
  geom_point(alpha = 0.5, size = 0.7, color = "#33a02c") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.7, span = 1) +
  facet_wrap(~ paper_code, scales = "free_y") + scale_y_log10() +
  labs(x = "age (months)", y = "RT (ms, log)", title = "LWL reaction time by dataset") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave(here("figs/data_checks/io_diag_rt.png"), pr, width = 11, height = 6, dpi = 130)

## ===== INPUT: rate per recording =====
in_eng <- read_csv(here("data/harmonized/input_level.csv"), show_col_types = FALSE) %>%
  transmute(paper_code, child_id = as.character(child_id), age, log_input)
sl_in <- read_csv(here("data/raw/WF2013/SLENA_LWL_LENA_n29.csv"), show_col_types = FALSE) %>%
  transmute(paper_code = "WF2013", child_id = as.character(SubNum), age = 18, log_input = log(FreqAWChr)) %>%
  filter(is.finite(log_input))
hb_in <- hb %>% transmute(child_id = as.character(ID), a18 = AD18AWChr, a25 = AD25AWChr) %>%
  pivot_longer(c(a18, a25), names_to = "tp", values_to = "awc") %>%
  transmute(paper_code = "Bang2025", child_id, age = ifelse(tp == "a18", 18, 25), log_input = log(awc)) %>%
  filter(is.finite(log_input))
inp <- bind_rows(in_eng, sl_in, hb_in) %>% filter(is.finite(log_input), age > 0) %>% mutate(paper_code = fct(paper_code))
pin <- ggplot(inp, aes(age, log_input)) +
  geom_line(aes(group = child_id), color = "grey75", alpha = 0.3, linewidth = 0.25) +
  geom_point(alpha = 0.5, size = 0.7, color = "#6a3d9a") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.7, span = 1) +
  facet_wrap(~ paper_code, scales = "free_y") +
  labs(x = "age (months)", y = "log input rate (study-specific units)", title = "Observed input by dataset") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave(here("figs/data_checks/io_diag_input.png"), pin, width = 11, height = 6, dpi = 130)

cat("wrote figs/data_checks/io_diag_{cdi,rt,input}.png\n")
cat("\ndatasets present per channel:\n")
cat(" CDI:  ", paste(levels(droplevels(cdi$paper_code)), collapse=", "), "\n")
cat(" RT:   ", paste(levels(droplevels(rt$paper_code)),  collapse=", "), "\n")
cat(" input:", paste(levels(droplevels(inp$paper_code)), collapse=", "), "\n")
