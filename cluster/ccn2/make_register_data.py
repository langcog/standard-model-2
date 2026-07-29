#!/usr/bin/env python3
"""Register/composition control data. Build 6 disjoint base corpora
(3 BabyLM non-CHILDES + 3 ClimbMix), then per (dataset,seed) NESTED-cumulative
developmental ladders (whitespace-word budgets) -- drop-in for train_gpt2_childes.py.
Datasets = fixed disjoint subsets (composition/sample); seed = init+order within a subset."""
import os, glob, random, json, sys
import pyarrow.parquet as pq
OUT_BASE="/data2/mcfrank/ladder/register_data"
OUT_CHUNK="/data2/mcfrank/ladder/register_chunks"
BABYLM="/data2/mcfrank/babylm/train_100M"
KEEP=["gutenberg","open_subtitles","simple_wiki","bnc_spoken","switchboard"]   # drop childes
CLIMB=os.path.expanduser("~/.cache/autoresearch/data")
RUNGS=[0.5,0.75,1,1.5,2,3,4,6,8,12,16,24]
SEEDS=[0,1,2,3,4,5,7,42]                 # 8 seeds
SPLIT_SEED=12345; N_SUB=3; TARGET=24_000_000
def wc(s): return len(s.split())
def babylm_docs():
    d=[]
    for s in KEEP:
        for ln in open(os.path.join(BABYLM,s+".train"),errors="ignore"):
            ln=ln.strip()
            if ln: d.append(ln)
    return d
def climb_docs(maxw):
    d=[]; tot=0
    for f in sorted(glob.glob(os.path.join(CLIMB,"*.parquet"))):
        for b in pq.ParquetFile(f).iter_batches(columns=["text"],batch_size=2000):
            for t in b.column("text").to_pylist():
                t=" ".join(t.split())
                if not t: continue
                d.append(t); tot+=wc(t)
                if tot>=maxw: return d
    return d
def split_disjoint(docs,n,target,seed):
    order=list(range(len(docs))); random.Random(seed).shuffle(order)
    subs=[[] for _ in range(n)]; subw=[0]*n
    for idx in order:
        cand=[i for i in range(n) if subw[i]<target]
        if not cand: break
        i=min(cand,key=lambda i: subw[i])
        subs[i].append(docs[idx]); subw[i]+=wc(docs[idx])
    return subs,subw
os.makedirs(OUT_BASE,exist_ok=True); os.makedirs(OUT_CHUNK,exist_ok=True)
datasets={}
bl=babylm_docs(); bs,bw=split_disjoint(bl,N_SUB,TARGET,SPLIT_SEED)
for i in range(N_SUB):
    n=f"babylm_S{i+1}"; datasets[n]=bs[i]; print(f"{n}: {bw[i]:,} words {len(bs[i]):,} docs")
cm=climb_docs(N_SUB*TARGET+6_000_000); cs,cw=split_disjoint(cm,N_SUB,TARGET,SPLIT_SEED)
for i in range(N_SUB):
    n=f"climbmix_C{i+1}"; datasets[n]=cs[i]; print(f"{n}: {cw[i]:,} words {len(cs[i]):,} docs")
man={}
for n,docs in datasets.items():
    p=os.path.join(OUT_BASE,n+".txt")
    with open(p,"w") as f:
        for x in docs: f.write(x+"\n")
    man[n]={"words":sum(wc(x) for x in docs),"docs":len(docs),"file":p}
targets=sorted(int(t*1e6) for t in RUNGS)
for n,docs in datasets.items():
    wcs=[wc(x) for x in docs]
    for seed in SEEDS:
        order=list(range(len(docs))); random.Random(seed).shuffle(order)
        dd=os.path.join(OUT_CHUNK,f"{n}_seed{seed}"); os.makedirs(dd,exist_ok=True)
        cum=0; ci=[]; ti=0
        for pos,di in enumerate(order):
            ci.append(di); cum+=wcs[di]
            while ti<len(targets) and (cum>=targets[ti] or pos==len(order)-1):
                mtag=f"{targets[ti]/1e6:g}M"
                with open(os.path.join(dd,f"{n}_seed{seed}_{mtag}.txt"),"w") as wf:
                    for i in ci: wf.write(docs[i]+"\n")
                ti+=1
            if ti>=len(targets): break
        print(f"[ladder] {n} seed{seed}: {ti} rungs top {cum:,}w")
json.dump(man,open(os.path.join(OUT_BASE,"register_manifest.json"),"w"),indent=2)
print("REGISTER_DATA_DONE")
