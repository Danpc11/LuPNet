#!/usr/bin/env python3
"""Build the inputs of the PH-IPF calculator from the in silico cohorts.

1. data/virtual_ph_cohort.tsv: one row per simulated case (final visit) from the calibrated cohort (200 patients,
   tag big2), the advanced-disease cohort (60, tag adv2) and the cohort without vasculopathy (60, tag hom).
2. results/ph_calculator.json:
   - tissue_curve: expected mPAP from tissue loss alone as a function of FVC (simulated cases without vasculopathy);
   - dlco_curve: expected log DLCO for a given FVC without vasculopathy, and its residual SD;
   - mixture: two-component Gaussian mixture of the DLCO residual (all patients with vasculopathy model);
   - logistic: P(mPAP >= 25 | FVC, DLCO), fitted with measurement noise (mPAP 2.5 mmHg, FVC and DLCO 5 points).

    python3 src/fit_calculator.py --cohorts results/ph_cohorts
Cohort files are produced by src/ph_cohort.py.
"""
import argparse, glob, json, os
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ap = argparse.ArgumentParser(); ap.add_argument("--cohorts", default=os.path.join(ROOT, "results", "ph_cohorts")); a = ap.parse_args()


def load(tag, group):
    from cohort_io import cases
    rows = []
    for _, d in cases(a.cohorts, tag):
        r0, x = d.iloc[0], d.iloc[-1]
        A = float(x["vasc_A"]) if "vasc_A" in d else 0.0
        rows.append(dict(cohort=group, mPAP=round(float(x.mPAP), 2), FVC=round(float(100 * x.FVC / r0.FVC), 1),
                         DLCO=round(float(100 * x.DLCO / r0.DLCO), 1), vasculopathy=int(A > 0), vasc_A=round(A, 3)))
    return pd.DataFrame(rows)


C = pd.concat([load("big2", "calibrated"), load("adv2", "advanced")], ignore_index=True)
C.to_csv(os.path.join(ROOT, "data", "virtual_ph_cohort.tsv"), sep="\t", index=False)
V = C.reset_index(drop=True)
H = V[V.vasculopathy == 0]
tissue = np.polyfit(H.FVC, H.mPAP, 2)
novasc = V[V.vasculopathy == 0]
dl = np.polyfit(novasc.FVC, np.log(novasc.DLCO), 2)
res = np.log(V.DLCO) - np.polyval(dl, V.FVC)
from sklearn.mixture import GaussianMixture
from sklearn.linear_model import LogisticRegression
gm = GaussianMixture(2, random_state=0).fit(res.values.reshape(-1, 1))
low = int(np.argmin(gm.means_.ravel()))
mix = dict(weights=gm.weights_.round(4).tolist(), means=gm.means_.ravel().round(4).tolist(),
           sds=np.sqrt(gm.covariances_.ravel()).round(4).tolist(), vascular_component=low)
rng = np.random.default_rng(0); Xs, ys = [], []
for _ in range(10):
    fv = V.FVC + rng.normal(0, 5, len(V)); d = (V.DLCO + rng.normal(0, 5, len(V))).clip(lower=5)
    Xs.append(np.c_[fv, np.log(fv / d)]); ys.append((V.mPAP + rng.normal(0, 2.5, len(V))) >= 25)
X = np.vstack(Xs); y = np.concatenate(ys)
lr = LogisticRegression().fit(X, y)
from sklearn.metrics import roc_auc_score
auc = roc_auc_score(y, lr.predict_proba(X)[:, 1])
out = dict(tissue_curve=tissue.round(6).tolist(), dlco_curve=dl.round(6).tolist(), dlco_residual_sd=round(float(res[V.vasculopathy == 0].std()), 4),
           mixture=mix, logistic=dict(intercept=round(float(lr.intercept_[0]), 4), coef_FVC=round(float(lr.coef_[0][0]), 5),
                                      coef_log_ratio=round(float(lr.coef_[0][1]), 4), auc_with_noise=round(float(auc), 3)),
           n_patients=int(len(V)), fvc_range=[float(V.FVC.min()), float(V.FVC.max())], dlco_range=[float(V.DLCO.min()), float(V.DLCO.max())])
json.dump(out, open(os.path.join(ROOT, "results", "ph_calculator.json"), "w"), indent=1)
print(json.dumps(out, indent=1)); print(f"cohort table: {len(C)} rows")
