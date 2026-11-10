#!/usr/bin/env python3
"""The same simulated cases with a tree of 256 and of 1,024 units (table S6).

The first N cases of the calibrated cohort (same parameters, lesion and stage) are simulated with G = 10
(1,024 units) and compared with the stored cohort, which uses G = 8 (256 units). A case with 1,024 units takes
several minutes of one core.

    python3 src/tree_size.py run --n 12 [--i0 0 --i1 12]      # one file per case, resumable, parallel-safe
    python3 src/tree_size.py analyze                          # tree_size_cases.tsv and table S6
On several cores: WHAT=tree bash scripts/cohorts_local.sh
"""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exercise_response import BIG2

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
WORK = os.path.join(ROOT, "results", "tree_size"); FIG = os.path.join(ROOT, "results", "figures")
CASES = os.path.join(FIG, "tree_size_cases.tsv")


def cmd_run(a):
    from ph_cohort import make_params, draws
    from lupnet.sim import run
    os.makedirs(WORK, exist_ok=True)
    stop, A, t0, ypu = draws(BIG2["n"], BIG2["seed"], BIG2["lo"], BIG2["hi"], BIG2["p_vasc"], BIG2["a_med"], BIG2["a_sig"])
    for i in range(max(a.i0, 0), min(a.n if a.i1 < 0 else a.i1, a.n)):
        f = os.path.join(WORK, f"case_{i:03d}.tsv")
        if os.path.exists(f):
            continue
        p = make_params(BIG2["pseed"] + i, a.G, float(stop[i]), float(A[i]), float(t0[i]), float(ypu[i]), BIG2["tau"], BIG2["cap"])
        rows, _, _ = run(p, do_exercise=False); x, r0 = rows[-1], rows[0]
        r = dict(case=i, G=a.G, vasc_A=A[i], healthy=x["healthy"], mPAP=x["mPAP"], FVC=100 * x["FVC"] / r0["FVC"], DLCO=100 * x["DLCO"] / r0["DLCO"])
        tmp = f + f".tmp{os.getpid()}"; pd.DataFrame([r]).to_csv(tmp, sep="\t", index=False); os.replace(tmp, f)
        print("tree-size case", i, flush=True)


def cmd_analyze(a):
    import glob
    from cohort_io import final_visits
    files = sorted(glob.glob(os.path.join(WORK, "case_*.tsv"))); os.makedirs(FIG, exist_ok=True)
    if files:
        pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True).sort_values("case").round(4).to_csv(CASES, sep="\t", index=False)
    if not os.path.exists(CASES):
        raise SystemExit("no cases: run `tree_size.py run` first")
    B = pd.read_csv(CASES, sep="\t").set_index("case")
    S = final_visits(os.path.join(ROOT, "results", "ph_cohorts"), "big2").set_index("case").loc[B.index]
    rows = []
    for k, unit in (("mPAP", "mmHg"), ("FVC", "points"), ("DLCO", "points")):
        d = B[k] - S[k]
        rows.append(dict(quantity=k, unit=unit, n=len(d), mean_difference=round(d.mean(), 2), sd=round(d.std(), 2), r=round(float(np.corrcoef(B[k], S[k])[0, 1]), 3),
                         largest=round(float(d.abs().max()), 2)))
    T = pd.DataFrame(rows); same = int(((B.mPAP >= 25) == (S.mPAP >= 25)).sum())
    T.to_csv(os.path.join(FIG, "TableS6_tree_size.tsv"), sep="\t", index=False)
    print(T.to_string(index=False)); print(f"differences are 1,024 - 256 units; {same} of {len(B)} cases keep the same PH status; per case:")
    print(pd.DataFrame({"mPAP_1024": B.mPAP, "mPAP_256": S.mPAP, "FVC_1024": B.FVC, "FVC_256": S.FVC, "DLCO_1024": B.DLCO, "DLCO_256": S.DLCO}).round(1).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--n", type=int, default=12); r.add_argument("--G", type=int, default=10)
    r.add_argument("--i0", type=int, default=0); r.add_argument("--i1", type=int, default=-1)
    sub.add_parser("analyze"); a = ap.parse_args()
    cmd_run(a) if a.cmd == "run" else cmd_analyze(a)
