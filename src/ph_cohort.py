#!/usr/bin/env python3
"""Virtual IPF cohorts for pulmonary hypertension (PH), with a structural vasculopathy independent of fibrosis.

Each simulated case has
  - its own fibrosis course (seed) and speed (years per model unit scaled by a lognormal factor, sd 0.5),
  - with probability P_VASC, a bounded structural small-artery vasculopathy: a non-remodelable series lesion
    R_les = R_ta0 * A * (1 - exp(-(t_yr + t0) / tau)), A ~ lognormal(A_MED, A_SIG), t0 ~ U(0, 10) yr,
    with capillary loss 1 / (1 + cap_frac * A * (...)) in every unit,
  - one cross-sectional visit at a random stage (fraction of healthy tissue in [lo, hi]).
Model: shared target law (law="shared", b = 0.925), calibrated aged IPF parameters, 256-unit tree by default.

Calibrated values (v0.3.0, honeycomb compliance hc_C = 0.15; 200 simulated cases, all seven targets of mild-to-moderate IPF
inside the 95% intervals):
  P_VASC = 0.7, A_MED = 1.5, A_SIG = 1.1, tau = 2 yr, cap_frac = 0.05.
(v0.2.x used hc_C = 0.60 and P_VASC = 0.5.)

    python3 src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85
    python3 src/ph_cohort.py run --tag adv2 --n 60 --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45
    python3 src/ph_cohort.py analyze --tag big2
    python3 src/ph_cohort.py analyze --table data/virtual_ph_cohort.tsv --tag calibrated   # shipped summary table
The run is resumable: finished patients are skipped.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import lupnet
from lupnet import Params
from lupnet.sim import run

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results", "ph_cohorts")
BASE = os.path.join(ROOT, "data", "calibrated_ipf_params.json")


PATTERNS = {
    # calibrated: the manuscript cohorts (patchy, no subpleural predominance, 22% honeycomb share)
    "calibrated": {},
    # UIP: subpleural seed, spread by contiguity, heterogeneous vulnerability, honeycombing
    "uip": dict(sub_gain=3.0, kn=1.5, kn_s=3.0, at2_sigma=1.5, lam_hc=0.5),
    # fibrotic NSIP: homogeneous in space and time, little honeycombing, diffuse membrane thickening, less stiff scar
    "nsip": dict(sub_gain=0.0, kn=0.0, kn_s=0.0, at2_sigma=0.3, lam_hc=0.05, mem_field=0.0, fib_C=0.30, nbr_D=1.0),
}


def make_params(seed, G, stop, A, t0, ypu, tau, cap, vmode="independent", t_max=3.0, PB=760.0, pattern="calibrated",
                where="diffuse", mix=0.5, a_scale=1.0):
    d = json.load(open(BASE)); p = Params()
    for k, v in d.items():
        if hasattr(p, k):
            setattr(p, k, tuple(v) if isinstance(getattr(p, k), tuple) else v)
    p.law, p.b, p.G, p.seed, p.stop_healthy = "shared", 0.925, G, seed, stop
    p.vasc_A, p.vasc_tau, p.t_vasc0, p.cap_frac, p.vasc_mode = A, tau, t0, cap, vmode
    p.years_per_unit *= ypu; p.stop_mPAP = 1e9; p.t_max = t_max; p.PB = PB
    for k, v in PATTERNS[pattern].items():
        setattr(p, k, v)
    p.pattern = pattern; p.vasc_where, p.vasc_mix = where, mix
    p.vasc_A = p.vasc_A * a_scale             # a field lesion needs a larger nominal severity for the same pressure
    return p


def draws(n, seed, lo, hi, p_vasc, a_med, a_sig):
    g = np.random.default_rng(seed)
    stop = g.uniform(lo, hi, n); has = g.uniform(0, 1, n) < p_vasc
    A = np.exp(g.normal(np.log(a_med), a_sig, n)) * has
    t0 = g.uniform(0, 10, n); ypu = np.exp(g.normal(0, 0.5, n))
    return stop, A, t0, ypu


def cmd_run(a):
    os.makedirs(OUT, exist_ok=True)
    # The draws depend only on --n and --seed, so a shard computes exactly the cases the
    # whole run would have computed: several shards of one cohort can run side by side.
    stop, A, t0, ypu = draws(a.n, a.seed, a.lo, a.hi, a.p_vasc, a.a_med, a.a_sig)
    i1 = a.n if a.i1 < 0 else min(a.i1, a.n)
    for i in range(max(a.i0, 0), i1):
        f = os.path.join(OUT, f"{a.tag}_{i:03d}.tsv")
        if os.path.exists(f):
            continue
        p = make_params(a.pseed + i, a.G, float(stop[i]), float(A[i]), float(t0[i]), float(ypu[i]), a.tau, a.cap, a.vmode, a.t_max, a.PB, a.pattern, a.where, a.mix, a.a_scale)
        t = time.time(); rows, _, _ = run(p, do_exercise=False)
        if rows[-1]["healthy"] > stop[i] + 0.01:
            # the run ended before the drawn stage (t_max reached, or a stop criterion of sim.run);
            # 7 of 60 cases of the adv2 cohort of the manuscript end this way (t_max = 3.0)
            print(f"  warning: case {i} stopped at healthy = {rows[-1]['healthy']:.3f}, drawn stage {stop[i]:.3f}", flush=True)
        df = pd.DataFrame(rows); df["vasc_A"] = A[i]; df["t_vasc0"] = t0[i]; df["ypu"] = ypu[i]; df["stop"] = stop[i]; df["G"] = a.G; df["PB"] = a.PB; df["pattern"] = a.pattern; df["where"] = a.where; df["mix"] = a.mix; df["a_scale"] = a.a_scale; df["vasc_A_scaled"] = A[i] * a.a_scale
        # how the case was generated, so that no reader has to assume it
        for k in ("tau", "cap", "p_vasc", "a_med", "a_sig", "seed", "pseed", "lo", "hi", "vmode", "n", "t_max"):
            df[k] = getattr(a, k)
        df["lupnet_version"] = lupnet.__version__
        tmp = f + f".tmp{os.getpid()}"
        df.to_csv(tmp, sep="\t", index=False); os.replace(tmp, f)
        print(a.tag, i, round(time.time() - t), "s", flush=True)


def cmd_analyze(a):
    if a.table:
        T = pd.read_csv(a.table, sep="\t")
        R = T[T.cohort == a.tag][["mPAP", "FVC", "DLCO"]].reset_index(drop=True) if a.tag else T[["mPAP", "FVC", "DLCO"]]
    else:
        from cohort_io import final_visits
        R = final_visits(OUT, a.tag)
        if len(R) and R.died.any():
            print(f"  {int(R.died.sum())} of {len(R)} runs ended in a non-viable state (oxygen demand not met or output "
                  f"below 2.5 L/min); their last viable visit is used")
    if len(R) == 0:
        raise SystemExit("no patients found: run the cohort first, or pass --table data/virtual_ph_cohort.tsv")
    Rm = R[(R.FVC >= 50) & (R.FVC <= 90)].reset_index(drop=True)

    from ph_stats import bootstrap, point, TARGETS as T_
    print(f"{a.tag}: {len(R)} simulated cases, {len(Rm)} with FVC 50-90%")
    if len(Rm) >= 10:
        ci, pt = bootstrap(Rm), point(Rm)
        for k, v in T_.items():
            lo, hi = ci[k]
            print(f"  {k:7s} {pt[k]:7.2f}  (95% CI {lo:.2f} to {hi:.2f})  target {v}  {'inside' if lo <= v <= hi else 'OUTSIDE'}")
    else:
        print("  fewer than 10 cases in the FVC 50-90% band (advanced disease): the seven mild-to-moderate targets do not apply")
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
    r.add_argument("--p_vasc", type=float, default=0.7); r.add_argument("--a_med", type=float, default=1.5); r.add_argument("--a_sig", type=float, default=1.1)
    r.add_argument("--tau", type=float, default=2.0); r.add_argument("--cap", type=float, default=0.05)
    r.add_argument("--vmode", default="independent", choices=["independent", "fibrosis"])
    r.add_argument("--t_max", type=float, default=3.0, help="model time limit of a run (years = t_max x years_per_unit x ypu); 3.0 in the manuscript cohorts")
    r.add_argument("--PB", type=float, default=760.0, help="barometric pressure, mmHg (760 sea level, the manuscript cohorts; 585 Mexico City, 2,240 m)")
    r.add_argument("--pattern", default="calibrated", choices=list(PATTERNS), help="histological pattern: calibrated (manuscript), uip, nsip")
    r.add_argument("--where", default="diffuse", choices=["diffuse", "scar", "field", "mixed"],
                   help="where the lesion sits: diffuse (manuscript), scar only, field (scar and adjacent lung), mixed")
    r.add_argument("--mix", type=float, default=0.5, help="mixed: share of the lesion that is diffuse")
    r.add_argument("--a_scale", type=float, default=1.0, help="multiplies the drawn lesion severity")
    r.add_argument("--i0", type=int, default=0, help="first case of this shard")
    r.add_argument("--i1", type=int, default=-1, help="one past the last case of this shard; -1 means --n")
    z = sub.add_parser("analyze"); z.add_argument("--tag", default="")
    z.add_argument("--table", default="", help="summary table (data/virtual_ph_cohort.tsv); --tag selects its cohort column")
    a = ap.parse_args()
    cmd_run(a) if a.cmd == "run" else cmd_analyze(a)
