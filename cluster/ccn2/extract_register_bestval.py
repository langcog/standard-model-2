#!/usr/bin/env python3
"""Build register_bestval.csv from the 576 composition-control surprisal CSVs.
Per (corpus,subset,seed,budget,word): best-val surprisal = mean_nll at max logged step.
Dict-keyed per word so duplicate evals at the same step (on_train_end + last scheduled
step) collapse to one row -- same convention as extract_ladder_bestval.py."""
import csv, glob, os, re
OUT="/data2/mcfrank/ladder/register_bestval.csv"
RUNGW={"0.5M":500000,"0.75M":750000,"1M":1000000,"1.5M":1500000,"2M":2000000,"3M":3000000,
       "4M":4000000,"6M":6000000,"8M":8000000,"12M":12000000,"16M":16000000,"24M":24000000}
cells={}
for f in sorted(glob.glob("/data2/mcfrank/ladder/register_dev/surprisal_*.csv")):
    m=re.match(r"surprisal_(babylm_S\d|climbmix_C\d)_seed(\d+)_(\S+)\.csv$", os.path.basename(f))
    if not m: continue
    ds,seed,rung=m.group(1),int(m.group(2)),m.group(3)
    data=list(csv.DictReader(open(f)))
    if not data: continue
    mx=max(int(x["step"]) for x in data)
    cells[(ds,seed,rung)]={x["word"]: float(x["mean_nll"]) for x in data if int(x["step"])==mx}
with open(OUT,"w",newline="") as o:
    w=csv.writer(o, lineterminator="\n"); w.writerow(["corpus","subset","seed","words","rung","word","surprisal"])
    for (ds,seed,rung) in sorted(cells, key=lambda k:(k[0],k[1],RUNGW.get(k[2],0))):
        corpus="babylm" if ds.startswith("babylm") else "climbmix"
        for word in sorted(cells[(ds,seed,rung)]):
            w.writerow([corpus,ds,seed,RUNGW.get(rung,0),rung,word,"%.6f"%cells[(ds,seed,rung)][word]])
print("cells:",len(cells)," rows:",sum(len(v) for v in cells.values()))
