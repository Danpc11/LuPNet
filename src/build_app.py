#!/usr/bin/env python3
"""Build index.html (the PH-IPF calculator) from results/ph_calculator.json and data/virtual_ph_cohort.tsv.

    python3 src/fit_calculator.py && python3 src/build_app.py
"""
import json, os
import pandas as pd
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
C = json.load(open(os.path.join(ROOT, "results", "ph_calculator.json")))
V = pd.read_csv(os.path.join(ROOT, "data", "virtual_ph_cohort.tsv"), sep="\t")
pts = [[round(r.FVC, 1), round(r.DLCO, 1), round(r.mPAP, 1), int(r.vasculopathy)] for r in V.itertuples()]
t = open(os.path.join(ROOT, "src", "app_template.html")).read()
rep = {"__COEF__": json.dumps(C), "__COHORT__": json.dumps(pts), "__NPAT__": str(len(pts)), "__AUC__": f"{C['logistic']['auc_with_noise']:.2f}",
       "__FMIN__": f"{C['fvc_range'][0]:.0f}", "__FMAX__": f"{C['fvc_range'][1]:.0f}", "__DMIN__": f"{C['dlco_range'][0]:.0f}", "__DMAX__": f"{C['dlco_range'][1]:.0f}"}
for k, v in rep.items():
    t = t.replace(k, v)
assert "__" not in t, "unfilled placeholder"
open(os.path.join(ROOT, "index.html"), "w").write(t); print("index.html written,", len(pts), "simulated cases")
