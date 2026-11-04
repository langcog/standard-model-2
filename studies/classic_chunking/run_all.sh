#!/bin/sh
# Reproduce the CLASSIC <-> v2 item-difficulty comparison end to end.
# usage (from the repo root): sh studies/classic_chunking/run_all.sh /path/to/cmudict.dict
# cmudict.dict: https://github.com/cmusphinx/cmudict (file `cmudict.dict`)
# Needs R (childesr, dplyr, readr; network for CHILDES) and python3 (numpy, pandas).
set -e
CMU=${1:?path to cmudict.dict}
D=studies/classic_chunking
W=$D/work
mkdir -p $W $D/results
Rscript $D/00_export_v2_items.R $W
[ -f $W/childes_eng_na_utts.csv.gz ] || OUTDIR=$W Rscript $D/01_pull_childes.R
python3 $D/02_prep_data.py $W "$CMU"
cd $D
# tag  utterances  runs  age_lo  age_hi  [learning rate]
python3 03_run_experiment.py work 60000 40 10 36 main
python3 03_run_experiment.py work 15000 40 10 36 n15k
python3 03_run_experiment.py work 120000 40 10 36 n120k
python3 03_run_experiment.py work 60000 40 10 60 age60
python3 03_run_experiment.py work 60000 40 10 36 lr20 0.2
for t in main n15k n120k age60 lr20; do
  python3 analyze.py work $t > results/${t}_correlations.txt
  python3 05_cv_baselines.py work $t > results/${t}_cv_r2.txt
  if [ $t = age60 ]; then A='10 60'; else A='10 36'; fi
  python3 06_features.py work $t $A > results/${t}_features.txt
done
cp work/item_table_main.csv work/growth_main.csv results/
