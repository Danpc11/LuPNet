#!/usr/bin/env python3
"""Checks added by the audit of v0.3.7 (tables S17-S20 in results/figures/). None of them changes a published number;
they state what the published numbers rest on.

S17  Where the healthy pressure drop comes from. Of mPAP - LAP in the healthy 256-unit lung, the shared law sets the
     arterial conduits, the venous conduits and the terminal sub-tree; the capillary sheet and the venules keep the
     resistances of the reference calibration (dP0 x split). The table gives each part in mmHg.
S18  The aged lung without fibrosis. The calibrated parameters age the vessels (age_M = 0.3, age_sigma = 0.3) after the
     healthy baseline, so a case with no fibrosis and no lesion already sits above the healthy lung. This is the
     reference that separates ageing from fibrosis in Figure 2g.
S19  How far the lesion has run at the visit, case by case (1 - exp(-(years + t0) / tau)), for big2 and adv2. With
     tau = 2 years and t0 ~ U(0, 10), the lesion is near its final value in almost every case, so the cohort results
     identify A and the prevalence, not tau or t0.
S20  The low-DLCO phenotype under the inter-individual spread of the healthy baseline used in table S6 (SD 13% FVC,
     15% DLCO, correlation 0.5) plus measurement noise (5 points): delta BIC, size of the low component and its
     sensitivity and specificity for PH, one row per replicate.

    python3 src/audit_extra.py
Importable: aged_baseline() is used by make_figures_erj.py for the reference line of Figure 2g.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lupnet import Params
from lupnet.lung import Lung
from lupnet.septal import make_lung

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results", "figures")
BASE = os.path.join(ROOT, "data", "calibrated_ipf_params.json")
VASC_TAU = 2.0


def calibrated_params(G=8):
    d = json.load(open(BASE)); p = Params()
    for k, v in d.items():
        if hasattr(p, k):
            setattr(p, k, tuple(v) if isinstance(getattr(p, k), tuple) else v)
    p.law, p.b, p.G = "shared", 0.925, G
    return p


def pressure_decomposition(G=8):
    """mmHg of the healthy pressure drop in each part of the circuit (shared law, no ageing)."""
    L = Lung(Params(law="shared", b=0.925, G=G)); vs = L.base_state["vs"]; p = L.p; N = L.N
    Ra = L._vres(vs["r_a"], L.Lv); Rv = L._vres(vs["r_v"], L.Lv)
    w = (L.nleaf / N) ** 2
    parts = dict(arterial_conduits=float((Ra * w).sum()) * p.Q, venous_conduits=float((Rv * w).sum()) * p.Q,
                 terminal_subtree=float(L._vres(vs["r_ta"], p.L_term).mean()) / N * p.Q,
                 capillary_sheet=L.R_cap0 / N * p.Q, venules=L.R_vn / N * p.Q)
    parts["total"] = sum(parts.values()); parts["mPAP"] = vs["PPA"]
    parts["set_by_law_pct"] = 100 * (parts["arterial_conduits"] + parts["venous_conduits"] + parts["terminal_subtree"]) / parts["total"]
    return parts


def aged_baseline(G=8):
    """mPAP of the calibrated lung with the vessels aged (age_M, age_sigma) and no fibrosis, no surfactant loss, no lesion."""
    q = calibrated_params(G)
    q.lam0_direct = 0.0; q.lam_ci = 0.0; q.lam_ci_cyc = 0.0; q.lam_ae = 0.0; q.lam_s = 0.0
    L = make_lung(q); b = L.base_state
    L._rem_budget = None
    st = L.equilibrate()                     # aging_on is True after the baseline: vessels move to z = 1 + age_M
    st = L.equilibrate()
    return dict(mPAP_healthy=b["vs"]["PPA"], mPAP=st["vs"]["PPA"], PVR=st["vs"]["dP"] / st["Q"], zerr=st["zerr"],
                age_M=q.age_M, age_sigma=q.age_sigma)


def lesion_saturation(cohorts=("big2", "adv2")):
    from cohort_io import cases
    rows = []
    for tag in cohorts:
        for i, d in cases(os.path.join(ROOT, "results", "ph_cohorts"), tag):
            x = d.iloc[-1]
            if float(x.vasc_A) <= 0:
                continue
            rows.append(dict(cohort=tag, case=i, vasc_A=float(x.vasc_A), t_vasc0=float(x.t_vasc0), years=float(x.years),
                             saturation=1 - np.exp(-(float(x.years) + float(x.t_vasc0)) / VASC_TAU), mPAP=float(x.mPAP)))
    return pd.DataFrame(rows)


def phenotype(D):
    """Two-component mixture of the DLCO residual (expected DLCO fitted to all cases), as in table S9."""
    from sklearn.mixture import GaussianMixture
    co = np.polyfit(D.FVC, np.log(D.DLCO), 2)
    res = (np.log(D.DLCO) - np.polyval(co, D.FVC)).values.reshape(-1, 1)
    g1 = GaussianMixture(1, random_state=0).fit(res); g2 = GaussianMixture(2, random_state=0).fit(res)
    low = g2.predict(res) == int(np.argmin(g2.means_.ravel())); ph = (D.mPAP >= 25).values
    return dict(dBIC=g1.bic(res) - g2.bic(res), n_low=int(low.sum()), PH_in_low_pct=100 * ph[low].mean() if low.any() else np.nan,
                sens_pct=100 * (low & ph).sum() / max(ph.sum(), 1), spec_pct=100 * (~low & ~ph).sum() / max((~ph).sum(), 1),
                residual_sd=float(res.std()))


def phenotype_with_baseline_spread(n_rep=20, seed=0):
    C = pd.read_csv(os.path.join(ROOT, "data", "virtual_ph_cohort.tsv"), sep="\t")
    C = C[C.cohort == "calibrated"].reset_index(drop=True)
    rows = [dict(replicate="noiseless", **phenotype(C))]
    rng = np.random.default_rng(seed)
    for rep in range(n_rep):
        z = rng.multivariate_normal([0, 0], [[1, 0.5], [0.5, 1]], len(C)); D = C.copy()
        D["FVC"] = D.FVC * (1 + 0.13 * z[:, 0]) + rng.normal(0, 5, len(D))
        D["DLCO"] = (D.DLCO * (1 + 0.15 * z[:, 1]) + rng.normal(0, 5, len(D))).clip(lower=3)
        rows.append(dict(replicate=rep, **phenotype(D)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    P = pressure_decomposition()
    pd.DataFrame([dict(part=k, mmHg=round(v, 3)) for k, v in P.items()]).to_csv(os.path.join(OUT, "TableS17_pressure_decomposition.tsv"), sep="\t", index=False)
    print("S17 healthy pressure drop (mmHg):", {k: round(v, 2) for k, v in P.items()})
    A = aged_baseline()
    pd.DataFrame([A]).round(3).to_csv(os.path.join(OUT, "TableS18_aged_baseline.tsv"), sep="\t", index=False)
    print(f"S18 aged lung, no fibrosis: mPAP {A['mPAP']:.2f} (healthy {A['mPAP_healthy']:.2f}), PVR {A['PVR']:.2f} WU")
    S = lesion_saturation()
    S.round(4).to_csv(os.path.join(OUT, "TableS19_lesion_saturation.tsv"), sep="\t", index=False)
    for tag, g in S.groupby("cohort"):
        print(f"S19 {tag}: lesion at visit, saturation min {g.saturation.min():.2f}, 5th pct {g.saturation.quantile(0.05):.2f}, "
              f"median {g.saturation.median():.2f}; years to visit median {g.years.median():.1f} (IQR {g.years.quantile(.25):.1f}-{g.years.quantile(.75):.1f})")
    R = phenotype_with_baseline_spread()
    R.round(3).to_csv(os.path.join(OUT, "TableS20_phenotype_baseline_spread.tsv"), sep="\t", index=False)
    noiseless, reps = R.iloc[0], R.iloc[1:]
    print(f"S20 phenotype: noiseless dBIC {noiseless.dBIC:.1f}, sens {noiseless.sens_pct:.0f}%, spec {noiseless.spec_pct:.0f}%; "
          f"with baseline spread + noise ({len(reps)} replicates): dBIC median {reps.dBIC.median():.1f} "
          f"[{reps.dBIC.min():.1f}, {reps.dBIC.max():.1f}], dBIC > 10 in {100 * (reps.dBIC > 10).mean():.0f}% of replicates, "
          f"sens {reps.sens_pct.mean():.0f}%, spec {reps.spec_pct.mean():.0f}%")
    print("tables S17-S20 written to results/figures/")
