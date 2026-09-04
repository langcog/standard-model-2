## Data figure (two column): the observed data entering the io-proc analysis.
## Rows = data streams (A CDI vocabulary, B LWL reaction time, C observed
## input), columns = the six datasets; one line per child, thick line = the
## per-dataset loess; empty cells = channels a dataset did not collect. Row
## colours are the stream hues shared with the graphical model and the
## io-proc figure (model/scripts/data_check_io_panels.R builds it).
##
## Reads: fits/joint_io_proc_english_count_subset_data.rds (for b$lwl) and the
##        raw CDI / LENA / BabyView sources listed in data_check_io_panels.R.
##        NOTE: this is the one figure that needs data/ + the fitted bundle
##        locally; everything else reads paper/cache/ only.
## Run:   Rscript paper/figures/fig_data.R
source(here::here("paper", "figures", "_common.R"))
source(here("model", "scripts", "data_check_io_panels.R"))

p <- io_panels_fig(pal = c(vocab = VOCAB_COL, processing = CHANNEL_PAL[["Processing"]],
                           input = CHANNEL_PAL[["Input"]]),
                   base_theme = theme_fig(), txt = TXT)
save_fig(tag_panels(p), "data", width = PNAS_2COL, height = 3.7)
