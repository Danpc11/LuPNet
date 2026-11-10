#!/usr/bin/env python3
"""Synthetic clinical file with the columns that noninvasive_tool.py expects (data/synthetic_patients.tsv).

The rows are NOT patients. They are drawn from the simulated visits of the calibrated and advanced cohorts, with a
referral bias towards PH, baseline spread and measurement noise, and optional measurements generated from fixed
linear models of mPAP. The file exists so that the scripts that need a clinical file can be run and tested without
patient data; no number obtained from it is a result.

    python3 src/synthetic_clinical.py            # writes data/synthetic_patients.tsv (60 rows, seed 0)
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import noninvasive_tool as nt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "synthetic_patients.tsv")
# measurement = intercept + slope x mPAP + N(0, sd); round numbers, not fitted to anyone
MODELS = {"sPAP_echo": (30.0, 0.6, 12.0), "PA_diam": (29.0, 0.25, 4.0), "logBNP": (3.5, 0.03, 0.9), "RV_FAC": (45.0, -0.25, 7.0), "RA_area": (14.0, 0.12, 4.0)}

if __name__ == "__main__":
    rng = np.random.default_rng(0); n = 60
    R = nt.reference(nt.trajectory_points(tags=("big2", "adv2")))
    w = np.where(R.mPAP >= 25, 3.0, 1.0); V = R.sample(n, weights=w, random_state=0).reset_index(drop=True)     # referral bias
    z = rng.multivariate_normal([0, 0], [[1, 0.5], [0.5, 1]], n); m = V.mPAP.values + rng.normal(0, 2.5, n)
    D = pd.DataFrame({nt.COLS["ID"]: [f"SYN-{i:03d}" for i in range(n)], nt.COLS["dx"]: "synthetic",
                      nt.COLS["pattern"]: rng.choice(["NIU", "NINE", "otro"], n, p=[0.25, 0.5, 0.25]),
                      nt.COLS["FVC"]: np.clip(V.FVC * (1 + 0.13 * z[:, 0]) + rng.normal(0, 5, n), 15, 120).round(0),
                      nt.COLS["DLCO"]: np.clip(V.DLCO * (1 + 0.15 * z[:, 1]) + rng.normal(0, 5, n), 8, 110).round(0),
                      nt.COLS["mPAP"]: m.round(0), nt.COLS["PVR"]: np.clip(V.PVR + rng.normal(0, 0.5, n), 0.5, None).round(1),
                      nt.COLS["PAWP"]: rng.integers(5, 14, n)})
    for k, (a, b, sd) in MODELS.items():
        x = a + b * m + rng.normal(0, sd, n)
        D[nt.COLS["BNP" if k == "logBNP" else k]] = np.exp(x).round(0) if k == "logBNP" else x.round(1)
    D.loc[rng.choice(n, 8, replace=False), nt.COLS["PA_diam"]] = np.nan          # some measurements are missing, as in practice
    D.loc[rng.choice(n, 6, replace=False), nt.COLS["BNP"]] = np.nan
    D.to_csv(OUT, sep="\t", index=False); print(f"{OUT}: {n} synthetic rows, {int((D[nt.COLS['mPAP']] >= 25).sum())} with mPAP >= 25")
