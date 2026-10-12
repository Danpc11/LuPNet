#!/usr/bin/env python3
"""Validation of the research calculator (index.html) on simulated cases held out from its reference set.

For 130 held-out cases (one random visit each, with measurement noise), the calculator's likelihood weighting is
applied with FVC and DLCO and with each optional input (resting SpO2, alveolar-arterial gradient, echocardiographic
systolic pressure). Prints the areas under the ROC curve for PH and relevant vasculopathy and the mean error of mPAP;
these are the figures quoted in the calculator.

    python3 src/validate_calculator.py
"""
import json, os, sys
import numpy as np
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import cases

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
C = json.load(open(os.path.join(ROOT, "results", "ph_calculator.json")))
R = np.array(C["reference"], float)
ids = []
for tag in ("big2", "adv2"):
    for i, d in cases(os.path.join(ROOT, "results", "ph_cohorts"), tag):
        s = d[(d.converged == 1) & (d.healthy >= 0.12) & (d.healthy <= 0.85)]
        if len(s):
            ids += [f"{tag}{i}"] * min(25, len(s))
ids = np.array(ids)
assert len(ids) == len(R)
sp = (R[:, 2] - 2) / 0.61
SDs = dict(FVC=5, DLCO=5, sPAP=8, AaDO2=6, SpO2=2)
test = np.random.default_rng(0).choice(np.unique(ids), 130, replace=False)
for name, use in (("FVC + DLCO", ()), ("+ SpO2", ("SpO2",)), ("+ A-a gradient", ("AaDO2",)), ("+ echo sPAP", ("sPAP",)),
                  ("+ echo sPAP + A-a gradient", ("sPAP", "AaDO2")), ("+ all three", ("sPAP", "AaDO2", "SpO2"))):
    g = np.random.default_rng(1)
    P, V, M, T, Tv, Tm = [], [], [], [], [], []
    for c in test:
        k = g.choice(np.where(ids == c)[0]); ref = ids != c
        obs = dict(FVC=R[k, 0] + g.normal(0, 5), DLCO=R[k, 1] + g.normal(0, 5), sPAP=sp[k] + g.normal(0, 8),
                   AaDO2=R[k, 3] + g.normal(0, 6), SpO2=R[k, 5] + g.normal(0, 2))
        ll = -0.5 * (((R[ref, 0] - obs["FVC"]) / 5) ** 2 + ((R[ref, 1] - obs["DLCO"]) / 5) ** 2)
        col = dict(sPAP=sp[ref], AaDO2=R[ref, 3], SpO2=R[ref, 5])
        for u in use:
            ll += -0.5 * ((col[u] - obs[u]) / SDs[u]) ** 2
        w = np.exp(ll - ll.max()); w /= w.sum()
        P.append((w * (R[ref, 2] >= 25)).sum()); V.append((w * R[ref, 4]).sum()); M.append((w * R[ref, 2]).sum())
        T.append(R[k, 2] >= 25); Tv.append(R[k, 4]); Tm.append(R[k, 2])
    print(f"{name:28s} AUC PH {roc_auc_score(T, P):.2f} | AUC relevant vasculopathy {roc_auc_score(Tv, V):.2f} | "
          f"mean error mPAP {np.mean(np.abs(np.array(M) - np.array(Tm))):.1f} mmHg")
