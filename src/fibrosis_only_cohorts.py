#!/usr/bin/env python3
"""Cohorts without vasculopathy (supplement S5; main text, "Fibrosis alone ...").

hom: 60 simulated cases in which vessels respond only to fibrosis and hypoxia (uniform fibrosis speed).
het: the same cases, but with the gain of pressure-driven maladaptive remodelling (eta_p) and the maximum strength
     of hypoxic vasoconstriction (hpv_max) drawn per case (lognormal, medians 0.3 and 3.0, SD of log 1.0 and 0.4).
     The manuscript uses the first 38 cases of this cohort.
Stage of the single visit: healthy tissue U(0.30, 0.85). Tree: 256 units; shared law with b = 0.925.

    python3 src/fibrosis_only_cohorts.py --arm hom --n 60
    python3 src/fibrosis_only_cohorts.py --arm het --n 38
Runs are resumable; results go to results/ph_cohorts/{arm}_NN.tsv (pack with cohort_io.pack).
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lupnet import Params
from lupnet.sim import run

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results", "ph_cohorts")
ap = argparse.ArgumentParser(); ap.add_argument("--arm", choices=["hom", "het"], required=True); ap.add_argument("--n", type=int, default=60)
a = ap.parse_args()
N = 60
rng = np.random.default_rng(2026)
eta = np.exp(rng.normal(np.log(0.3), 1.0, N)); hpv = np.exp(rng.normal(np.log(3.0), 0.4, N)); stop = rng.uniform(0.3, 0.85, N)
os.makedirs(OUT, exist_ok=True)
for i in range(a.n):
    f = os.path.join(OUT, f"{a.arm}_{i:02d}.tsv")
    if os.path.exists(f):
        continue
    e, h = (eta[i], hpv[i]) if a.arm == "het" else (0.0, 3.0)
    d = json.load(open(os.path.join(ROOT, "data", "calibrated_ipf_params.json"))); p = Params()
    for k, v in d.items():
        if hasattr(p, k):
            setattr(p, k, tuple(v) if isinstance(getattr(p, k), tuple) else v)
    p.law, p.b, p.G, p.seed, p.eta_p, p.hpv_max, p.beta, p.stop_healthy = "shared", 0.925, 8, 100 + i, float(e), float(h), 0.0, float(stop[i])
    t = time.time(); rows, _, _ = run(p, do_exercise=False)
    df = pd.DataFrame(rows); df["eta_p"] = e; df["hpv_max"] = h; df["stop"] = stop[i]; df["vasc_A"] = 0.0
    df.to_csv(f, sep="\t", index=False); print(a.arm, i, round(time.time() - t), "s", flush=True)
