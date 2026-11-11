#!/usr/bin/env python3
"""Re-simulate the shipped cohorts with the current code and compare them with the stored ones.

Every case is simulated again from the command recorded in data/cohorts_manifest.tsv and each numeric column of
every visit is compared with the stored cohort. One small file per case goes to results/verify/, so the run is
resumable and parallel-safe.

    python3 src/verify_cohorts.py run --jobs 96                    # every case of the 21 cohorts (about 35 core-hours)
    python3 src/verify_cohorts.py run --sample 2 --jobs 8          # two cases per cohort
    python3 src/verify_cohorts.py run --cohorts big2 hom --cases 0 5
    python3 src/verify_cohorts.py summary                          # table of the largest difference per cohort
"""
import argparse, os, shlex, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import cases, MANUSCRIPT, DESIGN

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); F = os.path.join(ROOT, "results", "ph_cohorts")
WORK = os.path.join(ROOT, "results", "verify"); MAN = os.path.join(ROOT, "data", "cohorts_manifest.tsv")
DEFAULT = dict(DESIGN["big2"], G=8, p_vasc=MANUSCRIPT["p_vasc"], a_med=MANUSCRIPT["a_med"], a_sig=MANUSCRIPT["a_sig"],
               tau=MANUSCRIPT["tau"], cap=MANUSCRIPT["cap"], vmode="independent", t_max=3.0, PB=760.0, pattern="calibrated", where="diffuse", mix=0.5, a_scale=1.0)


def arguments(command):
    """Options of a `ph_cohort.py run` command of the manifest, over the defaults."""
    a = dict(DEFAULT); tok = shlex.split(command)
    for k, v in zip(tok, tok[1:]):
        if k.startswith("--") and k[2:] in a:
            a[k[2:]] = type(DEFAULT[k[2:]])(v)
    return a


def simulate(tag, command, i):
    if command.startswith("fibrosis_only"):
        import fibrosis_only_cohorts
        return fibrosis_only_cohorts.simulate(tag, i)
    from ph_cohort import draws, make_params
    from lupnet.sim import run
    a = arguments(command); stop, A, t0, ypu = draws(a["n"], a["seed"], a["lo"], a["hi"], a["p_vasc"], a["a_med"], a["a_sig"])
    p = make_params(a["pseed"] + i, a["G"], float(stop[i]), float(A[i]), float(t0[i]), float(ypu[i]), a["tau"], a["cap"], a["vmode"], a["t_max"],
                    a["PB"], a["pattern"], a["where"], a["mix"], a["a_scale"])
    rows, _, _ = run(p, do_exercise=False); df = pd.DataFrame(rows); df["vasc_A"] = A[i]; df["t_vasc0"] = t0[i]
    return df


def one(task):
    tag, command, i = task; f = os.path.join(WORK, f"{tag}__{i:03d}.tsv")
    if os.path.exists(f):
        return
    stored = {k: d for k, d in cases(F, tag)}[i]; new = simulate(tag, command, i)
    cols = [c for c in stored.columns if c in new.columns and pd.api.types.is_numeric_dtype(stored[c]) and pd.api.types.is_numeric_dtype(new[c])]
    same_len = len(stored) == len(new); worst, wcol = (np.inf, "visits") if not same_len else (0.0, "")
    if same_len:
        for c in cols:
            x, y = stored[c].to_numpy(float), new[c].to_numpy(float); ok = ~(np.isnan(x) & np.isnan(y))
            d = float(np.nanmax(np.abs(x[ok] - y[ok]) / np.maximum(1.0, np.abs(x[ok])))) if ok.any() else 0.0
            if not d <= worst:
                worst, wcol = d, c
    r = dict(cohort=tag, case=i, visits_stored=len(stored), visits_new=len(new), columns=len(cols), max_rel_diff=worst, column=wcol,
             mPAP_stored=float(stored.mPAP.iloc[-1]), mPAP_new=float(new.mPAP.iloc[-1]))
    tmp = f + f".tmp{os.getpid()}"; pd.DataFrame([r]).to_csv(tmp, sep="\t", index=False); os.replace(tmp, f)
    print(f"{tag} case {i}: max relative difference {worst:.2e}", flush=True)


def cmd_run(a):
    os.makedirs(WORK, exist_ok=True); M = pd.read_csv(MAN, sep="\t"); rng = np.random.default_rng(0); tasks = []
    for r in M.itertuples():
        if a.cohorts and r.cohort not in a.cohorts:
            continue
        idx = a.cases if a.cases else (sorted(rng.choice(r.cases, min(a.sample, r.cases), replace=False).tolist()) if a.sample else range(r.cases))
        tasks += [(r.cohort, r.command, int(i)) for i in idx if int(i) < r.cases]
    print(len(tasks), "cases on", a.jobs, "workers", flush=True)
    if a.jobs > 1:
        from multiprocessing import Pool
        with Pool(a.jobs) as pool:
            list(pool.imap_unordered(one, tasks, chunksize=1))
    else:
        for t in tasks:
            one(t)
    cmd_summary(a)


def cmd_summary(a):
    import glob
    files = sorted(glob.glob(os.path.join(WORK, "*__*.tsv")))          # one file per case: cohort__case.tsv
    if not files:
        raise SystemExit("no verified cases: run `verify_cohorts.py run` first")
    D = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    S = D.groupby("cohort").agg(cases=("case", "count"), identical=("max_rel_diff", lambda x: int((x < 1e-6).sum())), max_rel_diff=("max_rel_diff", "max"),
                                max_mPAP_diff=("mPAP_new", lambda x: float(np.abs(x - D.loc[x.index, "mPAP_stored"]).max()))).reset_index()
    M = pd.read_csv(MAN, sep="\t").set_index("cohort").cases
    complete = set(S.cohort) == set(M.index) and all(S.set_index("cohort").cases.reindex(M.index) == M)
    # the record in data/ is the full verification; a partial run never replaces it
    out = os.path.join(ROOT, "data", "cohorts_verification.tsv") if complete else os.path.join(WORK, "summary.tsv")
    S.to_csv(out, sep="\t", index=False); print(S.to_string(index=False))
    print(f"{len(D)} of {int(M.sum())} cases; {int((D.max_rel_diff < 1e-6).sum())} identical within 1e-6; {os.path.relpath(out, ROOT)} written"
          + ("" if complete else " (partial run: data/cohorts_verification.tsv left as it is)"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--jobs", type=int, default=1); r.add_argument("--sample", type=int, default=0)
    r.add_argument("--cohorts", nargs="*", default=[]); r.add_argument("--cases", nargs="*", type=int, default=[])
    sub.add_parser("summary"); a = ap.parse_args()
    cmd_run(a) if a.cmd == "run" else cmd_summary(a)
