"""Sensitivity of the clinical benchmarks to the parameters of the lesion.

The lesion has five parameters, chosen jointly against seven aggregate targets, so
other combinations may fit those targets equally well. This script maps that
degeneracy for the two parameters the conclusion rests on: the time constant of the
lesion (T_les, --tau) and the proportion of cases that carry one (--p_vasc).

It does not simulate. It reads cohorts already generated at each parameter value and
reports, for each one, the quantities the manuscript calibrates against. Generate the
cohorts first (about 25 s per simulated case on one core):

    for tau in 1.0 3.0; do
        python3 src/ph_cohort.py run --tag tau${tau} --n 200 --seed 5000 --pseed 20000 \
            --lo 0.25 --hi 0.85 --tau $tau
    done
    for p in 0.50 0.85; do
        python3 src/ph_cohort.py run --tag pv${p} --n 200 --seed 5000 --pseed 20000 \
            --lo 0.25 --hi 0.85 --p_vasc $p
    done
    python3 src/sensitivity.py --ref big2 --tags tau1.0 tau3.0 pv0.50 pv0.85 \
        --labels "T_les 1 y" "T_les 3 y" "70% to 50%" "70% to 85%"

Keeping --seed and --pseed fixed means every cohort has the same fibrosis courses and
the same severities, so a difference between rows comes from the parameter and not
from a different draw of cases.

Writes results/figures/TableS16_sensitivity.tsv and FigureS2.{png,pdf}.
"""
import argparse
import os

import numpy as np
import pandas as pd

from cohort_io import cases
import ph_stats

MM = 1 / 25.4
BLUE, RED, GREY, DARK = "#2166AC", "#B2182B", "#8C8C8C", "#252525"


def mild(tag, cohorts):
    """One row per simulated case, last visit, restricted to FVC 50-90% predicted."""
    rows = []
    for _, d in cases(cohorts, tag):
        r0, x = d.iloc[0], d.iloc[-1]
        rows.append(dict(mPAP=x.mPAP, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO,
                         A=float(x.vasc_A)))
    R = pd.DataFrame(rows)
    if R.empty:
        raise SystemExit(f"cohort {tag!r} not found in {cohorts}")
    R["vasc"] = R.A > 0
    R["PH"] = R.mPAP >= 25
    return R[(R.FVC >= 50) & (R.FVC <= 90)].reset_index(drop=True)


def row(tag, label, cohorts):
    M = mild(tag, cohorts)
    if len(M) < 10:
        raise SystemExit(f"cohort {tag!r} has only {len(M)} cases with FVC 50-90% predicted: "
                         "too few for the benchmarks. Generate more cases (--n) and re-run.")
    p = ph_stats.point(M)
    _, lo, hi = ph_stats.prevalence_ci(M, lambda D: D.mPAP >= 25)
    return dict(label=label, tag=tag, n=len(M), vasc_pct=100 * M.vasc.mean(),
                ph_pct=p["ge25"], ph_lo=lo, ph_hi=hi, severe_pct=p["ge35"],
                rFVC=p["rFVC"], rDLCO=p["rDLCO"], dlco_ratio=p["ratio"],
                mPAP_no_lesion=M.mPAP[~M.vasc].mean(), mPAP_lesion=M.mPAP[M.vasc].mean())


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--ref", default="big2", help="calibrated cohort, the reference row")
    a.add_argument("--tags", nargs="+", required=True, help="cohorts generated at other parameter values")
    a.add_argument("--labels", nargs="+", help="one label per tag; defaults to the tags")
    a.add_argument("--cohorts", default="results/ph_cohorts")
    a.add_argument("--out", default="results/figures")
    a = a.parse_args()
    labels = a.labels or a.tags
    if len(labels) != len(a.tags):
        raise SystemExit("--labels must give one label per tag")

    R = pd.DataFrame([row(a.ref, "Calibrated", a.cohorts)] +
                     [row(t, l, a.cohorts) for t, l in zip(a.tags, labels)])
    os.makedirs(a.out, exist_ok=True)
    cols = ["label", "n", "vasc_pct", "ph_pct", "severe_pct", "rFVC", "rDLCO", "dlco_ratio"]
    R[cols].round(2).to_csv(os.path.join(a.out, "TableS16_sensitivity.tsv"), sep="\t", index=False)
    print(R[cols].round(2).to_string(index=False))

    # PH prevalence against the clinical target, anchored at zero
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 600, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 62 * MM), constrained_layout=True)
    y = np.arange(len(R))[::-1]
    ax = axes[0]
    ax.axvline(19, color=DARK, lw=0.8, ls=":", label="ARTEMIS-IPF")
    ax.barh(y, R.ph_pct, 0.6, color=[DARK] + [BLUE] * (len(R) - 1))
    ax.set(yticks=y, xlabel="PH (mPAP ≥25 mmHg), % of simulated cases", xlim=(0, max(40, R.ph_pct.max() * 1.15)))
    ax.set_yticklabels(R.label)
    ax.legend(loc="lower right", handlelength=1.2)
    ax.text(-0.18, 1.02, "A", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    ax = axes[1]
    ax.axvline(0, color=DARK, lw=0.6)
    ax.axvline(-0.30, color=DARK, lw=0.8, ls=":")
    ax.barh(y + 0.17, R.rDLCO, 0.33, color=RED, label="r(mPAP, DLCO)")
    ax.barh(y - 0.17, R.rFVC, 0.33, color=GREY, label="r(mPAP, FVC)")
    ax.set(yticks=y, xlabel="Correlation coefficient", xlim=(-0.8, 0.3))
    ax.set_yticklabels(R.label)
    ax.legend(loc="lower right", handlelength=1.2)
    ax.text(-0.18, 1.02, "B", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, f"FigureS2.{ext}"), bbox_inches="tight", pad_inches=0.04)
    print(f"wrote {a.out}/TableS16_sensitivity.tsv and FigureS2")


if __name__ == "__main__":
    main()
