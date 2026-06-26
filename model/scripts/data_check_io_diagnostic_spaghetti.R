## IO diagnostic: per-child vocabulary spaghetti by CDI form, faceted by dataset.
## Now driven by the harmonized item-level table (data/harmonized/cdi_item_level.csv,
## built by build_harmonized_cdi.R) so this doubles as a CHECK on that file. Uses the
## canonical-vocab set (item_canonical == TRUE) so WG/WS denominators are comparable
## across datasets. RUN LOCALLY.  Output: figs/data_checks/io_diagnostic_spaghetti.png
##
## babyview items are mapped through the webCDI dictionary in parse_babyview_cdi.R, so it
## now lands at WG=393 / WS=676 (vs canonical 396/680); the ~7 residual are homonym
## disambiguations (dress/toy/telephone) flagged for the crosswalk hand-fix.
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(ggplot2) })

harm <- read_csv(here("data/harmonized/cdi_item_level.csv"), show_col_types = FALSE) %>%
  filter(item_canonical)

## dedup to one row per (dataset, child, age, form, item); collapse same-age duplicate
## admins (produced-if-produced-in-any), then count vocabulary per admin.
adm <- harm %>%
  group_by(paper_code, child_id, age, form, item) %>%
  summarise(produces = max(produces, na.rm = TRUE), .groups = "drop") %>%
  group_by(paper_code, child_id, age, form) %>%
  summarise(vocab = sum(produces == 1, na.rm = TRUE), n_items = n(), .groups = "drop")
cat("form sizes (canonical items) by dataset:\n")
print(adm %>% distinct(paper_code, form, n_items) %>% arrange(paper_code) %>% as.data.frame(), row.names = FALSE)

DSORD <- c("babyview", "SEEDLingS", "AM2018", "FM2012", "FMW2013", "FPM2006")
adm <- adm %>% mutate(paper_code = factor(paper_code, levels = DSORD))
FPAL <- c("WG" = "#1f78b4", "WS" = "#e31a1c")
p <- ggplot(adm, aes(age, vocab, color = form)) +
  geom_line(aes(group = child_id), color = "grey70", alpha = 0.4, linewidth = 0.25) +
  geom_point(alpha = 0.6, size = 0.7) +
  facet_wrap(~ paper_code, nrow = 1) +
  scale_color_manual(values = FPAL, name = "Form") +
  labs(x = "age (months)", y = "vocabulary count (canonical items)",
       title = "IO diagnostic: per-kid vocabulary by form (lines connect a child's admins across age)",
       subtitle = "harmonized cdi_item_level.csv, item_canonical (production-vocab) only — all cohorts at WG=396/WS=680 canonical") +
  theme_bw(base_size = 9) + theme(legend.position = "top", panel.grid.minor = element_blank())
dir.create(here("figs", "data_checks"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figs", "data_checks", "io_diagnostic_spaghetti.png"), p, width = 14, height = 3.6, dpi = 130)
cat("wrote figs/data_checks/io_diagnostic_spaghetti.png\n")
