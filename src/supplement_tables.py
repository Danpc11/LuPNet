#!/usr/bin/env python3
"""Supplementary tables of the manuscript: S2 (morphometry by order), S4 (alternative vascular mechanism), S5 (2022
haemodynamic definitions) and the fibrosis-only cohorts of section S5. Requires the cohorts big2, adv2, altF, hom and het
in results/ph_cohorts.

    python3 src/supplement_tables.py
"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import pearsonr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lupnet.lung import womersley_G
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results", "figures"); os.makedirs(OUT, exist_ok=True)


def cohort(tag):
    from cohort_io import cases
    rows = []
    for _, d in cases(os.path.join(ROOT, "results", "ph_cohorts"), tag):
        r0, x = d.iloc[0], d.iloc[-1]
        rows.append(dict(mPAP=x.mPAP, PVR=x.PVR, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO, A=float(x.vasc_A)))
    return pd.DataFrame(rows)


def _unused_stats(D, rng):
    mp = D.mPAP + rng.normal(0, 2.5, len(D)); ph = D.mPAP >= 25
    return dict(lt20=100 * (D.mPAP < 20).mean(), b20=100 * ((D.mPAP >= 20) & (D.mPAP < 25)).mean(), ge25=100 * ph.mean(),
                ge35=100 * (D.mPAP >= 35).mean(), rD=pearsonr(mp, D.DLCO + rng.normal(0, 5, len(D)))[0],
                rF=pearsonr(mp, D.FVC + rng.normal(0, 5, len(D)))[0], ratio=D.DLCO[ph].mean() / D.DLCO[~ph].mean())


def noise_avg_r(D, col, seed=0):
    rng = np.random.default_rng(seed)
    return float(np.mean([pearsonr(D.mPAP + rng.normal(0, 2.5, len(D)), D[col] + rng.normal(0, 5, len(D)))[0] for _ in range(1000)]))


from ph_stats import bootstrap, point, prevalence_ci
from scipy.stats import spearmanr
# S4: calibrated (independent) vs vasculopathy linked to fibrosis
rows = []
for tag, name in (("big2", "independent"), ("altF", "linked to fibrosis")):
    D = cohort(tag); M = D[(D.FVC >= 50) & (D.FVC <= 90)].reset_index(drop=True)
    ci, pt = bootstrap(M), point(M)
    for k in pt:
        rows.append(dict(mechanism=name, n=len(M), n_vasculopathy=int((M.A > 0).sum()), quantity=k, value=round(pt[k], 6),
                         ci_low=round(ci[k][0], 6), ci_high=round(ci[k][1], 6)))
    print(f"{name}: realised prevalence of vasculopathy {100 * (D.A > 0).mean():.0f}% of {len(D)} cases ({int((M.A > 0).sum())}/{len(M)} with FVC 50-90%)")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "TableS4_alternative_mechanism.tsv"), sep="\t", index=False)
# S5: 2022 definitions (mPAP >= 25 in FVC 50-90% uses the same bootstrap as Table 1)
C, A = cohort("big2"), cohort("adv2"); M = C[(C.FVC >= 50) & (C.FVC <= 90)].reset_index(drop=True)
Adv = pd.concat([C, A])[lambda D: D.FVC < 50].reset_index(drop=True)
rows = []
for lab, D in (("FVC 50-90%", M), ("FVC < 50% (pooled big2 + adv2)", Adv)):
    for name, f in (("mPAP >= 25", lambda D: D.mPAP >= 25), ("mPAP > 20", lambda D: D.mPAP > 20),
                    ("mPAP > 20 and PVR > 2 WU", lambda D: (D.mPAP > 20) & (D.PVR > 2)), ("PVR > 5 WU", lambda D: D.PVR > 5)):
        v, lo, hi = prevalence_ci(D, f)
        if lab.startswith("FVC 50") and name == "mPAP >= 25":
            v, (lo, hi) = point(D)["ge25"], bootstrap(D)["ge25"]
        rows.append(dict(group=lab, n=len(D), definition=name, pct=round(v, 4), ci_low=round(lo, 4), ci_high=round(hi, 4)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "TableS5_2022_definitions.tsv"), sep="\t", index=False)
# Fibrosis alone (hom) and variable vascular responses (het): supplement S5 text
rows = []
for tag in ("hom", "het"):
    D = []
    from cohort_io import cases
    for _, d in cases(os.path.join(ROOT, "results", "ph_cohorts"), tag):
        r0, x = d.iloc[0], d.iloc[-1]
        D.append(dict(mPAP=x.mPAP, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO, eta=x.eta_p, hpv=x.hpv_max))
    D = pd.DataFrame(D); M = D[(D.FVC >= 50) & (D.FVC <= 90)]
    rows.append(dict(cohort=tag, n=len(D), n_fvc_50_90=len(M), mPAP_mean=round(M.mPAP.mean(), 2), mPAP_sd=round(M.mPAP.std(), 2),
                     mPAP_ge25=int((M.mPAP >= 25).sum()), r_FVC=round(pearsonr(M.mPAP, M.FVC)[0], 3), r_DLCO=round(pearsonr(M.mPAP, M.DLCO)[0], 3),
                     spearman_eta=round(spearmanr(M.mPAP, M.eta)[0], 3) if tag == "het" else None,
                     spearman_hpv=round(spearmanr(M.mPAP, M.hpv)[0], 3) if tag == "het" else None))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "TableS_fibrosis_only.tsv"), sep="\t", index=False)
# S2: morphometry by order (b = 0.90, phi1 = 0.8 arteries / 0.4 veins)
Q, RHO, OM, W1 = 5 / 60 * 1e-3, 1060.0, 2 * np.pi * 70 / 60, 2.865
def tree(f):
    H = pd.read_csv(os.path.join(ROOT, "data", f), sep="\t"); return H, Q / H.N.values, H.L_cm.values / 100, H.mu_cP.values * 1e-3, H.D_cm.values / 100
def radius(q, L, mu, b, t0, phi):
    g = lambda x: np.log(4 * mu * q / (np.pi * np.exp(3 * x)) * (1 + W1 * phi * womersley_G(np.exp(x) * np.sqrt(OM * RHO / mu))[0])) - np.log(t0 * np.exp(x) ** (b - 1) * L ** ((b - 1) / 2))
    return float(np.exp(brentq(g, np.log(1e-9), np.log(5.0))))
pred = lambda T, b, t0, phi: np.array([2 * radius(q, l, m, b, t0, phi) for q, l, m in zip(T[1], T[2], T[3])])
err = lambda T, D: float(np.sqrt(np.mean((np.log10(D) - np.log10(T[4])) ** 2)))
Ta, Tv = tree("huang1996_pulmonary_arteries.tsv"), tree("huang1996_pulmonary_veins.tsv")
t0 = 10 ** minimize_scalar(lambda lt: err(Ta, pred(Ta, 0.9, 10 ** lt, 0.8)) ** 2, bounds=(-3, 2), method="bounded", options=dict(xatol=1e-4)).x
rows = []
for name, T, phi in (("artery", Ta, 0.8), ("vein", Tv, 0.4)):
    D = pred(T, 0.9, t0, phi)
    for o, n, dm, dp in zip(T[0].order, T[0].N, T[4], D):
        rows.append(dict(tree=name, order=o, N=n, D_meas_mm=round(dm * 1e3, 3), D_pred_mm=round(dp * 1e3, 3), dev_pct=round(100 * (dp / dm - 1), 1)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "TableS2_morphometry.tsv"), sep="\t", index=False)
print("tables S2, S4, S5 and fibrosis-only written to results/figures/")
