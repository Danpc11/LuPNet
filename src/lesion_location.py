#!/usr/bin/env python3
"""Where the vascular lesion sits, and how much of it can be a local consequence of fibrosis (tables S25-S26).

The cohorts share the fibrosis course and the lesion draws of big2 and adv2 and differ only in the location of the
lesion (`ph_cohort.py run --where`):
  diffuse  big2, adv2              every unit (the manuscript cohorts)
  mix50    big2_mix50, adv2_mix50  half of the lesion diffuse, half following the field (--where mixed --mix 0.25 --a_scale 2.0)
  field    big2_field, adv2_field  scar and the healthy units next to it (--where field --a_scale 5.8)
  scar     big2_scar               fibrotic and honeycomb units only
`a_scale` gives each mechanism the PH prevalence of the diffuse cohort in mild-to-moderate disease, so that the
comparison is at equal prevalence. Runs that end in a non-viable state contribute their last viable visit and are
counted as deaths.

S25  the seven mild-to-moderate benchmarks, advanced disease (FVC < 50%, pooled) and deaths, by mechanism.
S26  PH (mPAP >= 25 mmHg) by FVC band, by mechanism, and in catheterised patients when --clinical is given.

    python3 src/lesion_location.py [--clinical Base_HP_con_Cate.xlsx]
"""
import argparse, os, sys, warnings
import pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import final_visits
from ph_stats import bootstrap, point, prevalence_ci, TARGETS

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
F = os.path.join(ROOT, "results", "ph_cohorts"); OUT = os.path.join(ROOT, "results", "figures")
MECH = {"diffuse": ("big2", "adv2"), "mix50": ("big2_mix50", "adv2_mix50"), "field": ("big2_field", "adv2_field"), "scar": ("big2_scar",)}
BANDS = [(20, 50), (50, 70), (70, 100)]

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--clinical", default=""); a = ap.parse_args(); os.makedirs(OUT, exist_ok=True)
    rows, bands = [], []
    for name, tags in MECH.items():
        if not os.path.exists(os.path.join(F, tags[0] + ".tsv.gz")):
            print(f"{name}: cohort {tags[0]} not found, skipped"); continue
        D = final_visits(F, tags[0]); M = D[(D.FVC >= 50) & (D.FVC <= 90)].reset_index(drop=True)
        pt, ci = point(M), bootstrap(M)
        r = dict(mechanism=name, n_mild=len(M), **{k: round(pt[k], 2) for k in TARGETS},
                 **{k + "_ci": f"{ci[k][0]:.2f} to {ci[k][1]:.2f}" for k in TARGETS},
                 targets_inside=int(sum(ci[k][0] <= v <= ci[k][1] for k, v in TARGETS.items())), deaths=int(D.died.sum()))
        P = pd.concat([D] + [final_visits(F, t) for t in tags[1:]], ignore_index=True); A5 = P[P.FVC < 50]
        if len(tags) > 1:
            ph = prevalence_ci(A5, lambda d: d.mPAP >= 25); sv = prevalence_ci(A5, lambda d: d.mPAP > 40)
            r.update(n_adv=len(A5), PH_adv=f"{ph[0]:.0f} ({ph[1]:.0f}-{ph[2]:.0f})", severe_adv=f"{sv[0]:.0f} ({sv[1]:.0f}-{sv[2]:.0f})",
                     deaths=int(P.died.sum()))
        rows.append(r)
        bands.append(dict(mechanism=name, **{f"FVC {lo}-{hi}": f"{100 * (P[(P.FVC >= lo) & (P.FVC < hi)].mPAP >= 25).mean():.0f}% (n = {((P.FVC >= lo) & (P.FVC < hi)).sum()})"
                                             for lo, hi in BANDS}))
    if a.clinical:
        from noninvasive_tool import load_clinical
        V = load_clinical(a.clinical).dropna(subset=["FVC", "mPAP"])
        bands.append(dict(mechanism="patients", **{f"FVC {lo}-{hi}": f"{int((V[(V.FVC >= lo) & (V.FVC < hi)].mPAP >= 25).sum())} of {((V.FVC >= lo) & (V.FVC < hi)).sum()}"
                                                   for lo, hi in BANDS}))
    S25, S26 = pd.DataFrame(rows), pd.DataFrame(bands)
    S25.to_csv(os.path.join(OUT, "TableS25_lesion_location.tsv"), sep="\t", index=False)
    S26.to_csv(os.path.join(OUT, "TableS26_PH_by_FVC_and_location.tsv"), sep="\t", index=False)
    show = ["mechanism", "n_mild", "ge25", "ge35", "rDLCO", "rFVC", "ratio", "targets_inside", "PH_adv", "severe_adv", "deaths"]
    print(S25[[c for c in show if c in S25]].to_string(index=False)); print(S26.to_string(index=False))
    print("targets:", TARGETS, "| advanced disease: PH 46%, mPAP > 40 mmHg about 9%")
