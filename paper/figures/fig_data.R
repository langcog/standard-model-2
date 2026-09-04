## Data figure (two column): the observed data entering the io-proc analysis,
## one colour per dataset -- (A) CDI vocabulary trajectories, (B) LWL reaction
## time, (C) observed input rate. Content unchanged from the manuscript's
## fig-io-data-panels chunk (model/scripts/data_check_io_panels.R).
##
## Reads: fits/joint_io_proc_english_count_subset_data.rds (for b$lwl) and the
##        raw CDI / LENA / BabyView sources listed in data_check_io_panels.R.
##        NOTE: this is the one figure that needs data/ + the fitted bundle
##        locally; everything else reads paper/cache/ only.
## Run:   Rscript paper/figures/fig_data.R
source(here::here("paper", "figures", "_common.R"))
source(here("model", "scripts", "data_check_io_panels.R"))

p <- io_panels_fig()
save_fig(p, "data", width = PNAS_2COL, height = PNAS_2COL * 0.4)
