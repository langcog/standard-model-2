## IO diagnostic: per-child vocabulary spaghetti by CDI form, faceted by dataset.
## Built from the FULL item sets (not the model subsample): raw sumscore (top)
## and proportion of the form's items (bottom; comparable across WG/WS, which
## have different denominators). RUN LOCALLY.
## Output: figs/data_checks/io_diagnostic_spaghetti.png
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(ggplot2); library(patchwork) })

## full per-admin CDI from the source long files (all items, item-level so we can dedup)
stanford <- read_csv(here("data/peekbank/stanford_cdi_items_long.csv"), show_col_types = FALSE) %>%
  transmute(ds = recode(study, totlot2 = "FM2012", totlot3 = "AM2018", tlo = "FMW2013"),
            id = as.character(lab_subject_id), age, form, item, produces)
totlot <- read_csv(here("data/peekbank/totlot_cdi_items_long.csv"), show_col_types = FALSE) %>%
  transmute(ds = "fernald_totlot", id = as.character(lab_subject_id), age, form = "WS", item, produces)
iop <- readRDS(here("fits/io_pooled_subset_data.rds"))$df %>%
  filter(study %in% c("BabyView", "SEEDLingS")) %>%
  transmute(ds = study, id = as.character(ckey), age, form = form, item, produces)

## dedup to one row per (child, age, form, item) -- collapse same-age duplicate
## admins / duplicate item rows (produced-if-produced-in-any), matching the
## bundle's dedup. Without this, retested kids show 2x the item denominator.
adm <- bind_rows(stanford, totlot, iop) %>%
  group_by(ds, id, age, form, item) %>% summarise(produces = max(produces, na.rm = TRUE), .groups = "drop") %>%
  group_by(ds, id, age, form) %>%
  summarise(vocab = sum(produces == 1, na.rm = TRUE), n_items = n(), .groups = "drop") %>%
  mutate(prop = vocab / n_items)
cat("form sizes (items) by dataset:\n")
print(adm %>% distinct(ds, form, n_items) %>% arrange(ds) %>% as.data.frame(), row.names = FALSE)

DSORD <- c("BabyView", "SEEDLingS", "AM2018", "FM2012", "FMW2013", "fernald_totlot")
adm <- adm %>% mutate(ds = factor(ds, levels = DSORD))
FPAL <- c("WG" = "#1f78b4", "WS" = "#e31a1c")
## raw vocab count only; common (shared) x-axis across facets for comparison
p <- ggplot(adm, aes(age, vocab, color = form)) +
  geom_line(aes(group = id), color = "grey70", alpha = 0.4, linewidth = 0.25) +
  geom_point(alpha = 0.6, size = 0.7) +
  facet_wrap(~ ds, nrow = 1, scales = "free_y") +          # shared x; free y (WG/WS denominators differ)
  scale_x_continuous(limits = range(adm$age, na.rm = TRUE)) +
  scale_color_manual(values = FPAL, name = "Form") +
  labs(x = "age (months)", y = "vocabulary count (raw)",
       title = "IO diagnostic: per-kid vocabulary by form (lines connect a child's admins across age)") +
  theme_bw(base_size = 9) + theme(legend.position = "top", panel.grid.minor = element_blank())
dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figs", "data_checks", "io_diagnostic_spaghetti.png"), p, width = 14, height = 3.6, dpi = 130)
cat("wrote figs/data_checks/io_diagnostic_spaghetti.png\n")
