## Standard-model-2 build targets.
##
## Usage (from project root):
##   make recovery         Parameter-recovery sim + fit.
##   make data             Prepare Wordbank subsample Stan data.
##   make baseline         Fit baseline on subsample (no diagnostic sweep).
##   make diagnostics      Fit the 2x2 diagnostic sweep + LOO comparison.
##   make analyze          Produce RQ-aligned plots from fits.
##   make analyze-all      Run analyze over all 4 variants.
##   make explainer        Recompile outputs/model_explainer.pdf from .tex.
##   make all              recovery -> data -> diagnostics -> analyze-all.
##   make clean-fits       Delete fits/ (force refit next run).
##   make clean-figs       Delete outputs/figs/.
##   make clean            Both.
##
## Data-size knobs for prepare_data (defaults 500 children x 200 items):
##   make data N_CHILDREN=1000 N_ITEMS=300

RSCRIPT   := Rscript
R_FLAGS   :=
SCRIPTS   := model/scripts

N_CHILDREN ?= 500
N_ITEMS    ?= 200

.PHONY: all smoke recovery data baseline diagnostics variant 2pl \
        analyze analyze-all explainer migrate clean clean-fits clean-figs

all: recovery data diagnostics analyze-all

smoke:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/00_smoke.R

recovery:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/01_recovery.R

data:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/02_prepare_data.R $(N_CHILDREN) $(N_ITEMS)

# Just the baseline variant of the 2x2 (reuses the diagnostic script).
baseline: data
	@echo "The baseline variant is fit as part of 'make diagnostics'."
	@echo "Use 'make diagnostics' to fit all four. For just baseline:"
	$(RSCRIPT) $(R_FLAGS) -e 'source("model/R/config.R"); \
	  source("model/R/helpers.R"); \
	  b <- readRDS(file.path(PATHS$$fits_dir, "subset_data.rds")); \
	  fit_variant(b$$stan_data, "wordbank_baseline")'

diagnostics: data
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/03_diagnostics.R

# Fit a single named variant: make variant NAME=2pl
variant: data
	@if [ -z "$(NAME)" ]; then echo "Usage: make variant NAME=<variant>"; exit 1; fi
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/fit_variant.R $(NAME)

# Convenience: fit 2PL (the next planned variant)
2pl: data
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/fit_variant.R 2pl

# ---- Longitudinal pipeline ----
# Pull item-level longitudinal data from Wordbank (needs network).
longitudinal-pull:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/pull_longitudinal.R

# Prepare a stratified subset for Stan: make longitudinal-data [LANG=English] [N_CHILDREN=600] [N_ITEMS=200]
longitudinal-data:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/prepare_longitudinal_data.R \
	    $${LANG:-English (American)} $${N_CHILDREN:-600} $${N_ITEMS:-200}

# Fit: make longitudinal-fit VARIANT=long_slopes [DATASET=norwegian]
longitudinal-fit:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/fit_longitudinal.R \
	    $${VARIANT:-long_slopes} $${DATASET:-english}

# Analyze: make longitudinal-analyze VARIANT=long_slopes [DATASET=norwegian]
longitudinal-analyze:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/analyze_longitudinal.R \
	    $${VARIANT:-long_slopes} $${DATASET:-english}

# ---- Input-observed pipeline (BabyView, eventually Seedlings) ----
babyview-data:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/prepare_babyview.R $${N_ITEMS:-200}

# Fit: make io-fit VARIANT=io_2pl_slopes [DATASET=babyview]
io-fit:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/fit_io.R \
	    $${VARIANT:-io_2pl_slopes} $${DATASET:-babyview}

# Sensitivity on sigma_r: make sensitivity VARIANT=2pl [SIGMA_R=0.3,0.53,0.8,1.0]
sensitivity: data
	@if [ -z "$(VARIANT)" ]; then echo "Usage: make sensitivity VARIANT=<name> [SIGMA_R=a,b,c]"; exit 1; fi
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/sensitivity_sigma_r.R $(VARIANT) $(SIGMA_R)

analyze:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/04_analyze.R baseline

analyze-all:
	$(RSCRIPT) $(R_FLAGS) $(SCRIPTS)/04_analyze.R all

explainer:
	cd notes && pdflatex -interaction=nonstopmode model_explainer.tex > /dev/null \
	    && pdflatex -interaction=nonstopmode model_explainer.tex > /dev/null \
	    ; rm -f model_explainer.aux model_explainer.log \
	            model_explainer.out model_explainer.toc

# One-off migration: move outputs from the pre-refactor flat layout into
# fits + outputs/figs. Safe to run any number of times (no-op if the
# source files are gone).
migrate:
	@mkdir -p fits outputs/figs
	@for v in baseline fix_delta fix_s both_fixed; do \
	  [ -f model/wordbank_fit_$$v.rds ] && \
	    mv model/wordbank_fit_$$v.rds fits/wordbank_$$v.rds && \
	    echo "moved wordbank_fit_$$v.rds -> fits/wordbank_$$v.rds" || true; \
	done
	@[ -f model/wordbank_loos.rds ] && \
	  mv model/wordbank_loos.rds fits/loo_compare.rds && \
	  echo "moved wordbank_loos.rds -> fits/loo_compare.rds" || true
	@[ -f model/wordbank_data_subset.rds ] && \
	  cp model/wordbank_data_subset.rds fits/subset_data.rds && \
	  echo "copied wordbank_data_subset.rds -> fits/subset_data.rds" || true
	@[ -f model/recovery_fit.rds ] && \
	  mv model/recovery_fit.rds fits/recovery.rds && \
	  echo "moved recovery_fit.rds -> fits/recovery.rds" || true
	@for f in model/*.png; do \
	  [ -f "$$f" ] && mv "$$f" outputs/figs/ && echo "moved $$f -> figs/" || true; \
	done

clean-fits:
	rm -rf fits

clean-figs:
	rm -rf outputs/figs

clean: clean-fits clean-figs
