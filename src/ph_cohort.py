#!/usr/bin/env python3
"""Virtual IPF cohorts for pulmonary hypertension (PH), with a structural vasculopathy independent of fibrosis.

Each simulated case has
  - its own fibrosis course (seed) and speed (years per model unit scaled by a lognormal factor, sd 0.5),
  - with probability P_VASC, a bounded structural small-artery vasculopathy: a non-remodelable series lesion
    R_les = R_ta0 * A * (1 - exp(-(t_yr + t0) / tau)), A ~ lognormal(A_MED, A_SIG), t0 ~ U(0, 10) yr,
    with capillary loss 1 / (1 + cap_frac * A * (...)) in every unit,
  - one cross-sectional visit at a random stage (fraction of healthy tissue in [lo, hi]).
Model: shared target law (law="shared", b = 0.925), calibrated aged IPF parameters, 256-unit tree by default.

Calibrated values (200-patient cohort, all seven targets of ARTEMIS-IPF / Zisman 2007 inside the 95% CI):
  P_VASC = 0.5, A_MED = 1.5, A_SIG = 1.1, tau = 2 yr, cap_frac = 0.05.

    python3 src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85
    python3 src/ph_cohort.py run --tag adv2 --n 60 --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45
    python3 src/ph_cohort.py analyze --tag big2
    python3 src/ph_cohort.py analyze --table data/virtual_ph_cohort.tsv --tag calibrated   # shipped summary table
The run is resumable: finished patients are skipped.
"""
import argparse, glob, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from lupnet import Params
from lupnet.sim import run

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results", "ph_cohorts")
BASE = os.path.join(ROOT, "data", "calibrated_ipf_params.json")
TARGETS = dict(lt20=51, b20_25=30, ge25=19, ge35=4, rDLCO=-0.30, rFVC=0.0, ratio=0.73)


def make_params(seed, G, stop, A, t0, ypu, tau, cap, vmode="independent"):
    d = json.load(open(BASE)); p = Params()
    for k, v in d.items():
        if hasattr(p, k):
            setattr(p, k, tuple(v) if isinstance(getattr(p, k), tuple) else v)
    p.law, p.b, p.G, p.seed, p.stop_healthy = "shared", 0.925, G, seed, stop
    p.vasc_A, p.vasc_tau, p.t_vasc0, p.cap_frac, p.vasc_mode = A, tau, t0, cap, vmode
    p.years_per_unit *= ypu; p.stop_mPAP = 1e9
    return p


def draws(n, seed, lo, hi, p_vasc, a_med, a_sig):
    g = np.random.default_rng(seed)
    stop = g.uniform(lo, hi, n); has = g.uniform(0, 1, n) < p_vasc
    A = np.exp(g.normal(np.log(a_med), a_sig, n)) * has
    t0 = g.uniform(0, 10, n); ypu = np.exp(g.normal(0, 0.5, n))
    return stop, A, t0, ypu


def cmd_run(a):
    os.makedirs(OUT, exist_ok=True)
    stop, A, t0, ypu = draws(a.n, a.seed, a.lo, a.hi, a.p_vasc, a.a_med, a.a_sig)
    for i in range(a.n):
        f = os.path.join(OUT, f"{a.tag}_{i:03d}.tsv")
        if os.path.exists(f):
            continue
        p = make_params(a.pseed + i, a.G, float(stop[i]), float(A[i]), float(t0[i]), float(ypu[i]), a.tau, a.cap, a.vmode)
        t = time.time(); rows, _, _ = run(p, do_exercise=False)
        df = pd.DataFrame(rows); df["vasc_A"] = A[i]; df["t_vasc0"] = t0[i]; df["ypu"] = ypu[i]; df["stop"] = stop[i]; df["G"] = a.G
        df.to_csv(f, sep="\t", index=False)
        print(a.tag, i, round(time.time() - t), "s", flush=True)


def cmd_analyze(a):
    from scipy.stats import pearsonr
    if a.table:
        T = pd.read_csv(a.table, sep="\t")
        R = T[T.cohort == a.tag][["mPAP", "FVC", "DLCO"]].reset_index(drop=True) if a.tag else T[["mPAP", "FVC", "DLCO"]]
    else:
        from cohort_io import cases
        rows = []
        for _, d in cases(OUT, a.tag):
            r0, x = d.iloc[0], d.iloc[-1]
            rows.append(dict(mPAP=x.mPAP, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO))
        R = pd.DataFrame(rows)
    if len(R) == 0:
        raise SystemExit("no patients found: run the cohort first, or pass --table data/virtual_ph_cohort.tsv")
    Rm = R[(R.FVC >= 50) & (R.FVC <= 90)].reset_index(drop=True)

    from ph_stats import bootstrap, point, TARGETS as T_
    ci, pt = bootstrap(Rm), point(Rm)
    print(f"{a.tag}: {len(R)} simulated cases, {len(Rm)} with FVC 50-90%")
    for k, v in T_.items():
        lo, hi = ci[k]
        print(f"  {k:7s} {pt[k]:7.2f}  (95% CI {lo:.2f} to {hi:.2f})  target {v}  {'inside' if lo <= v <= hi else 'OUTSIDE'}")
    adv = R[R.FVC < 50]
    if len(adv):
        print(f"  FVC < 50%: n = {len(adv)}, PH {100 * (adv.mPAP >= 25).mean():.0f}%, mPAP >= 40 {100 * (adv.mPAP >= 40).mean():.0f}%"
              "  (transplant candidates: 46%, ~9%)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--tag", required=True); r.add_argument("--n", type=int, default=200); r.add_argument("--seed", type=int, default=5000)
    r.add_argument("--pseed", type=int, default=20000, help="fibrosis seed of patient i is pseed + i")
    r.add_argument("--lo", type=float, default=0.25); r.add_argument("--hi", type=float, default=0.85); r.add_argument("--G", type=int, default=8)
    r.add_argument("--p_vasc", type=float, default=0.5); r.add_argument("--a_med", type=float, default=1.5); r.add_argument("--a_sig", type=float, default=1.1)
    r.add_argument("--tau", type=float, default=2.0); r.add_argument("--cap", type=float, default=0.05)
    r.add_argument("--vmode", default="independent", choices=["independent", "fibrosis"])
    z = sub.add_parser("analyze"); z.add_argument("--tag", default="")
    z.add_argument("--table", default="", help="summary table (data/virtual_ph_cohort.tsv); --tag selects its cohort column")
    a = ap.parse_args()
    cmd_run(a) if a.cmd == "run" else cmd_analyze(a)
