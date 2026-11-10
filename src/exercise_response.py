#!/usr/bin/env python3
"""Pressure-flow response of the calibrated cohort (Figure 6, table S16).

The first N cases of the calibrated cohort (big2: same seeds, lesions and stages) are simulated again up to their
visit, and the circulation is then re-equilibrated at higher cardiac outputs with the vessels held at the calibre
they had reached (no remodelling) and the breathing pattern, oxygen uptake and hypoxic vasoconstriction free to
respond. The flow is imposed: the right ventricle is not allowed to limit it. For reference, each case also records
what the resting right ventricle delivers when 10 L/min are demanded (`Q_rv`, `mPAP_rv`); versions of the manuscript
up to v21 reported that state as "10 L/min" (see docs/AUDIT.md, finding 21).

    python3 src/exercise_response.py run --n 70 [--i0 0 --i1 70]     # one file per case, resumable, parallel-safe
    python3 src/exercise_response.py analyze                         # pressure_flow_cases.tsv, table S16, Figure 6
On several cores: WHAT=flow bash scripts/cohorts_local.sh
"""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
WORK = os.path.join(ROOT, "results", "pressure_flow"); FIG = os.path.join(ROOT, "results", "figures")
CASES = os.path.join(FIG, "pressure_flow_cases.tsv")
FLOWS = (7.5, 10.0)                                   # L/min
BIG2 = dict(n=200, seed=5000, pseed=20000, lo=0.25, hi=0.85, p_vasc=0.7, a_med=1.5, a_sig=1.1, tau=2.0, cap=0.05)


def load_flow(L, Q, imposed=True):
    """mPAP and delivered output when a cardiac output Q is imposed (or only demanded of the resting right
    ventricle), without remodelling; the state of the lung is restored afterwards."""
    p = L.p
    keep = (L.h.copy(), L.ru_a.copy(), L.ru_v.copy(), L.ru_ta.copy(), L.x.copy(), L.rd.copy(), L._rem_budget, L.Ees, L.EDV, p.rv)
    L._rem_budget = None
    try:
        if imposed:
            p.rv = False
        st = L.equilibrate(Qtot=Q, VO2=p.VO2 * Q / p.Q, VA=p.VA * Q / p.Q, remodel=False) if imposed else L.equilibrate(Qtot=Q, remodel=False)
    finally:
        L.h, L.ru_a, L.ru_v, L.ru_ta, L.x, L.rd, L._rem_budget, L.Ees, L.EDV, p.rv = keep
    return float(st["vs"]["PPA"]), float(st["Q"]), float(st["gas"]["SaO2"])


def cmd_run(a):
    from ph_cohort import make_params, draws
    from lupnet.sim import run
    os.makedirs(WORK, exist_ok=True)
    stop, A, t0, ypu = draws(BIG2["n"], BIG2["seed"], BIG2["lo"], BIG2["hi"], BIG2["p_vasc"], BIG2["a_med"], BIG2["a_sig"])
    for i in range(max(a.i0, 0), min(a.n if a.i1 < 0 else a.i1, a.n)):
        f = os.path.join(WORK, f"case_{i:03d}.tsv")
        if os.path.exists(f):
            continue
        p = make_params(BIG2["pseed"] + i, 8, float(stop[i]), float(A[i]), float(t0[i]), float(ypu[i]), BIG2["tau"], BIG2["cap"])
        rows, _, L = run(p, do_exercise=False); x = rows[-1]
        r = dict(case=i, vasc_A=A[i], healthy=x["healthy"], FVC=100 * x["FVC"], DLCO=100 * x["DLCO"], Q_rest=x["Q"], mPAP_rest=x["mPAP"], PVR_rest=x["PVR"])
        for Q in FLOWS:
            m, q, s = load_flow(L, Q); r[f"mPAP_{Q:g}"] = m; r[f"Q_{Q:g}"] = q; r[f"SaO2_{Q:g}"] = s
        r["mPAP_rv"], r["Q_rv"], _ = load_flow(L, 10.0, imposed=False)
        tmp = f + f".tmp{os.getpid()}"; pd.DataFrame([r]).to_csv(tmp, sep="\t", index=False); os.replace(tmp, f)
        print("pressure-flow case", i, flush=True)


def cmd_analyze(a):
    import glob, logging, matplotlib
    matplotlib.use("Agg"); logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    import matplotlib.pyplot as plt
    files = sorted(glob.glob(os.path.join(WORK, "case_*.tsv"))); os.makedirs(FIG, exist_ok=True)
    if files:
        pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True).sort_values("case").round(4).to_csv(CASES, sep="\t", index=False)
    if not os.path.exists(CASES):
        raise SystemExit("no cases: run `exercise_response.py run` first")
    D = pd.read_csv(CASES, sep="\t"); D["vasc"] = D.vasc_A > 0
    D["slope"] = (D["mPAP_10"] - D.mPAP_rest) / (D["Q_10"] - D.Q_rest)
    D["slope_v21"] = (D.mPAP_rv - D.mPAP_rest) / 5.0                 # as reported up to manuscript v21
    rng = np.random.default_rng(1)
    def prop(x):
        if len(x) == 0:
            return "–"
        b = [100 * rng.choice(x, len(x)).mean() for _ in range(2000)]
        return f"{100 * x.mean():.0f}% ({np.percentile(b, 2.5):.0f}-{np.percentile(b, 97.5):.0f}); {int(x.sum())} of {len(x)}"
    iqr = lambda s: f"{s.median():.1f} ({s.quantile(.25):.1f}-{s.quantile(.75):.1f})"
    rows = []
    for v, lab in ((False, "No vasculopathy"), (True, "Vasculopathy")):
        g = D[D.vasc == v]; n = g[g.mPAP_rest < 25]
        rows.append({"Group": lab, "Simulated cases": len(g), "Resting mPAP (mmHg)": iqr(g.mPAP_rest), "mPAP at 7.5 L/min": iqr(g["mPAP_7.5"]),
                     "mPAP at 10 L/min": iqr(g["mPAP_10"]), "Slope (mmHg per L/min)": f"{g.slope.median():.2f} ({g.slope.quantile(.25):.2f}-{g.slope.quantile(.75):.2f})",
                     "Cases with resting mPAP < 25": len(n), "Of these, >= 25 mmHg at 7.5 L/min": prop((n["mPAP_7.5"] >= 25).values),
                     "Of these, >= 25 mmHg at 10 L/min": prop((n["mPAP_10"] >= 25).values), "Slope > 3": prop((g.slope > 3).values),
                     "Output delivered by the resting RV when 10 L/min are demanded": iqr(g.Q_rv),
                     "As in manuscript v21: mPAP, slope": f"{g.mPAP_rv.median():.1f}, {g.slope_v21.median():.2f}"})
    T = pd.DataFrame(rows); T.to_csv(os.path.join(FIG, "TableS16_pressure_flow.tsv"), sep="\t", index=False)
    print(T.T.to_string(header=False))
    # ---- Figure 6
    MM = 1 / 25.4
    plt.rcParams.update({"font.family": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False, "axes.linewidth": 0.6, "savefig.dpi": 600, "pdf.fonttype": 42})
    RED, GREY, DARK = "#B2182B", "#9E9E9E", "#252525"
    N = D[D.mPAP_rest < 25]
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 62 * MM), gridspec_kw=dict(wspace=0.5))
    tag = lambda ax, s: ax.text(-0.27, 1.06, s, transform=ax.transAxes, fontsize=9, fontweight="bold")
    ax = axs[0]
    for v, c in ((False, GREY), (True, RED)):
        for _, r in N[N.vasc == v].iterrows():
            ax.plot([r.Q_rest, r["Q_7.5"], r["Q_10"]], [r.mPAP_rest, r["mPAP_7.5"], r["mPAP_10"]], "-o", ms=1.6, lw=0.5, color=c, alpha=0.6, mew=0)
        ax.plot([], [], "-o", ms=2.5, lw=0.8, color=c, label=("Vasculopathy" if v else "No vasculopathy") + f" (n = {(N.vasc == v).sum()})")
    ax.axhline(25, color=DARK, lw=0.5, ls=":"); ax.set(xlabel="Cardiac output (L/min)", ylabel="mPAP (mmHg)", xticks=[5, 7.5, 10])
    ax.legend(loc="upper left", fontsize=5.6, handlelength=1.4); tag(ax, "A")
    ax = axs[1]; x = 0
    for Q in FLOWS:
        for v, c in ((False, GREY), (True, RED)):
            y = (N[N.vasc == v][f"mPAP_{Q:g}"] >= 25).values
            b = [100 * rng.choice(y, len(y)).mean() for _ in range(2000)] if len(y) else [0]
            ax.bar(x, 100 * y.mean() if len(y) else 0, 0.8, color=c)
            ax.errorbar(x, 100 * y.mean() if len(y) else 0, yerr=[[100 * y.mean() - np.percentile(b, 2.5)], [np.percentile(b, 97.5) - 100 * y.mean()]] if len(y) else None,
                        fmt="none", ecolor=DARK, elinewidth=0.7, capsize=1.5)
            x += 1
        x += 0.6
    ax.set_xticks([0.5, 3.1]); ax.set_xticklabels(["7.5 L/min", "10 L/min"]); ax.set(ylabel="Cases reaching 25 mmHg (%)", ylim=(0, 105)); tag(ax, "B")
    ax = axs[2]
    for v, c in ((False, GREY), (True, RED)):
        g = D[D.vasc == v]; ax.plot(g.mPAP_rest, g.slope, "o", ms=2.6, color=c, mew=0, alpha=0.8)
    co = np.polyfit(D.mPAP_rest, D.slope, 1); xx = np.linspace(D.mPAP_rest.min(), D.mPAP_rest.max(), 20); ax.plot(xx, np.polyval(co, xx), color=DARK, lw=0.9)
    ax.axhline(3, color=DARK, lw=0.5, ls=":"); ax.set(xlabel="Resting mPAP (mmHg)", ylabel="mPAP / cardiac output slope\n(mmHg per L/min)"); tag(ax, "C")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"V21_Figure6.{ext}"), bbox_inches="tight", pad_inches=0.04)
    print(f"{len(D)} cases; {CASES}, TableS16_pressure_flow.tsv and V21_Figure6 written")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--n", type=int, default=70); r.add_argument("--i0", type=int, default=0); r.add_argument("--i1", type=int, default=-1)
    sub.add_parser("analyze"); a = ap.parse_args()
    cmd_run(a) if a.cmd == "run" else cmd_analyze(a)
