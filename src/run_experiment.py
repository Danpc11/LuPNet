#!/usr/bin/env python3
"""One progression run -> results/runs/<tag>.tsv and <tag>_snap.npz.

    python3 src/run_experiment.py --beta 3 --seed 1 --tag b3_s1
    python3 src/run_experiment.py --beta 3 --seed 1 --no-remodel --tag b3_s1_norem
Any Params field can be set with --set name=value.
"""
import argparse, os, sys, time, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lupnet import Params
from lupnet.sim import run

ap = argparse.ArgumentParser()
ap.add_argument("--beta", type=float, default=3.0)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--no-remodel", action="store_true")
ap.add_argument("--no-hpv", action="store_true")
ap.add_argument("--set", nargs="*", default=[])
ap.add_argument("--tag", required=True)
ap.add_argument("--out", default="results/runs")
a = ap.parse_args()
p = Params(beta=a.beta, seed=a.seed, remodel=not a.no_remodel, hpv=not a.no_hpv)
for kv in a.set:
    k, v = kv.split("=")
    cur = getattr(p, k)
    setattr(p, k, v == "True" if isinstance(cur, bool) else v if isinstance(cur, str) else type(cur)(float(v)))
os.makedirs(a.out, exist_ok=True)
t0 = time.time()
rows, snaps, L = run(p)
df = pd.DataFrame(rows)
df.to_csv(os.path.join(a.out, f"{a.tag}.tsv"), sep="\t", index=False)
np.savez_compressed(os.path.join(a.out, f"{a.tag}_snap.npz"), ix=L.ix, iy=L.iy,
                    **{f"{k}_{q}": v for q, s in snaps.items() for k, v in s.items()})
json.dump(p.dict(), open(os.path.join(a.out, f"{a.tag}_params.json"), "w"), indent=1)
print(f"{a.tag}: {len(df)} steps, final healthy {df.healthy.iloc[-1]:.3f}, {time.time()-t0:.0f} s")
