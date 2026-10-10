#!/usr/bin/env python3
"""Supplementary figure S1 and tables S6-S11: range of lung function and sensitivity analyses of the audit.

a) FVC and DLCO (% of each case's healthy baseline) along the full trajectories of the calibrated cohort, against the
   fraction of healthy tissue; the shaded band is the window from which the single cross-sectional visit was drawn.
b) DLCO against FVC for all trajectory points (grey) and for the cross-sectional visits (coloured by mPAP).
Table S6: calibration targets after adding inter-individual variation of the healthy baseline in % predicted
(SD 13% for FVC and 15% for DLCO, correlation 0.5; mean of 20 replicates).

    python3 src/make_supplement_figures.py
"""
import logging, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg"); logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import cases
from ph_stats import point

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results", "figures"); os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4
plt.rcParams.update({"font.family": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False, "axes.linewidth": 0.6, "savefig.dpi": 600, "pdf.fonttype": 42})
BLUE, GREEN, GREY, DARK = "#2166AC", "#1B7837", "#BDBDBD", "#252525"
traj, visits = [], []
for i, d in cases(os.path.join(ROOT, "results", "ph_cohorts"), "big2"):
    r0, x = d.iloc[0], d.iloc[-1]
    traj.append(pd.DataFrame(dict(case=i, healthy=d.healthy, FVC=100 * d.FVC / r0.FVC, DLCO=100 * d.DLCO / r0.DLCO)))
    visits.append(dict(FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO, mPAP=x.mPAP, A=float(x.vasc_A)))
T = pd.concat(traj, ignore_index=True); V = pd.DataFrame(visits)
edges = np.linspace(0.25, 1.0, 16); mid = 0.5 * (edges[1:] + edges[:-1])
q = lambda col, p: [T[col][(T.healthy >= a) & (T.healthy < b)].quantile(p) for a, b in zip(edges[:-1], edges[1:])]
fig, axs = plt.subplots(1, 2, figsize=(180 * MM, 68 * MM), gridspec_kw=dict(wspace=0.35))
ax = axs[0]
ax.axvspan(0.25, 0.85, color="#F0F0F0", lw=0, zorder=0)
ax.text(0.55, 104, "window of the cross-sectional visit", ha="center", fontsize=6, color=DARK)
for col, c in (("FVC", BLUE), ("DLCO", GREEN)):
    ax.fill_between(mid, q(col, 0.25), q(col, 0.75), color=c, alpha=0.2, lw=0)
    ax.plot(mid, q(col, 0.5), color=c, label=f"{col}, median and IQR")
ax.set(xlim=(1.0, 0.25), ylim=(0, 110), xlabel="Healthy tissue remaining (fraction)", ylabel="% of healthy baseline")
ax.legend(loc="lower left"); ax.text(-0.15, 1.04, "a", transform=ax.transAxes, fontsize=9, fontweight="bold")
ax = axs[1]
ax.plot(T.FVC, T.DLCO, ".", ms=1.2, color=GREY, alpha=0.4, zorder=0, label="Trajectory points")
sc = ax.scatter(V.FVC, V.DLCO, c=V.mPAP.clip(15, 35), cmap="RdBu_r", vmin=15, vmax=35, s=7, lw=0, label="Cross-sectional visits")
cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03, ticks=[15, 20, 25, 30, 35]); cb.set_label("mPAP (mmHg)", fontsize=6.5)
cb.ax.tick_params(labelsize=6, width=0.5, length=2); cb.outline.set_linewidth(0.5)
ax.set(xlim=(35, 102), ylim=(0, 105), xlabel="FVC (% of healthy baseline)", ylabel="DLCO (% of healthy baseline)")
ax.legend(loc="upper left", markerscale=3); ax.text(-0.15, 1.04, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"FigureS1.{ext}"), bbox_inches="tight", pad_inches=0.04)
# table S6: baseline variation in % predicted
C = V.copy(); rows = []
for rep in range(20):
    g = np.random.default_rng(100 + rep); z = g.multivariate_normal([0, 0], [[1, 0.5], [0.5, 1]], len(C))
    D = C.copy(); D["FVC"] = D.FVC * (1 + 0.13 * z[:, 0]); D["DLCO"] = D.DLCO * (1 + 0.15 * z[:, 1])
    M = D[(D.FVC >= 50) & (D.FVC <= 90)].reset_index(drop=True); p = point(M); p["n"] = len(M); rows.append(p)
R = pd.DataFrame(rows)
S6 = pd.DataFrame(dict(quantity=R.columns, mean=R.mean().round(3).values, min=R.min().round(3).values, max=R.max().round(3).values))
S6.to_csv(os.path.join(OUT, "TableS6_baseline_variation.tsv"), sep="\t", index=False)
bins = [(0.95, 1.01), (0.85, 0.95), (0.75, 0.85), (0.6, 0.75), (0.45, 0.6), (0.25, 0.45)]
S7 = pd.DataFrame([dict(healthy=f"{a:.2f}-{min(b, 1):.2f}", FVC_median=round(T.FVC[(T.healthy >= a) & (T.healthy < b)].median(), 1),
                        DLCO_median=round(T.DLCO[(T.healthy >= a) & (T.healthy < b)].median(), 1)) for a, b in bins])
S7.to_csv(os.path.join(OUT, "TableS7_trajectory_ranges.tsv"), sep="\t", index=False)
print(S6.to_string(index=False)); print(S7.to_string(index=False))


# ---------------------------------------------------------------- audit sensitivity analyses (tables S8-S11)
from sklearn.mixture import GaussianMixture
from cohort_io import cases as _cases


def _visits(tag):
    R = []
    for _, d in _cases(os.path.join(ROOT, "results", "ph_cohorts"), tag):
        r0, x = d.iloc[0], d.iloc[-1]
        R.append(dict(mPAP=x.mPAP, PVR=x.PVR, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO, A=float(x.vasc_A), hc=x.honeycomb))
    return pd.DataFrame(R)


CAL, ADV = _visits("big2"), _visits("adv2")
MILD = CAL[(CAL.FVC >= 50) & (CAL.FVC <= 90)].reset_index(drop=True)
# S8: pulmonary artery wedge pressure (left atrial pressure fixed at 8 mmHg in the model)
rows = [dict(scenario="Fixed at 8 mmHg (model)", **point(MILD))]
for lab, mu, sd, post in (("Mean 8, SD 2.5 mmHg", 8, 2.5, 0.0), ("Mean 7, SD 2.5 mmHg", 7, 2.5, 0.0), ("Fixed at 9 mmHg", 9, 0.0, 0.0),
                          ("Mean 9, SD 2.5 mmHg, 5% post-capillary (16-22 mmHg)", 9, 2.5, 0.05)):
    out = []
    for r in range(20):
        g = np.random.default_rng(300 + r); M = MILD.copy()
        pawp = np.clip(g.normal(mu, sd, len(M)), 4, 15)
        post_ = g.uniform(0, 1, len(M)) < post; pawp[post_] = g.uniform(16, 22, post_.sum())
        M["mPAP"] = M.mPAP + (pawp - 8); out.append(point(M))
    rows.append(dict(scenario=lab, **pd.DataFrame(out).mean().to_dict()))
pd.DataFrame(rows).round(3).to_csv(os.path.join(OUT, "TableS8_wedge_pressure.tsv"), sep="\t", index=False)
# S9: low-DLCO phenotype without using the vasculopathy labels (expected DLCO from a trimmed fit to all cases)
x, y = CAL.FVC.values, np.log(CAL.DLCO.values); ph = (CAL.mPAP >= 25).values; rows = []
masks = {"Expected DLCO from cases without vasculopathy (manuscript)": CAL.A.values == 0, "Expected DLCO from all cases": np.ones(len(CAL), bool)}
m = np.ones(len(CAL), bool)
for _ in range(5):
    c = np.polyfit(x[m], y[m], 2); r = y - np.polyval(c, x); mad = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = r > -2 * mad
masks["Expected DLCO from a trimmed fit to all cases"] = m
for lab, mk in masks.items():
    c = np.polyfit(x[mk], y[mk], 2); res = (y - np.polyval(c, x)).reshape(-1, 1)
    g2 = GaussianMixture(2, random_state=0).fit(res); g1 = GaussianMixture(1, random_state=0).fit(res)
    L = g2.predict(res) == int(np.argmin(g2.means_.ravel()))
    rows.append(dict(method=lab, dBIC=round(g1.bic(res) - g2.bic(res), 1), n_low=int(L.sum()), PH_in_low_pct=round(100 * ph[L].mean(), 1),
                     sensitivity_pct=round(100 * (L & ph).sum() / ph.sum(), 1), specificity_pct=round(100 * ((~L) & (~ph)).sum() / (~ph).sum(), 1)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "TableS9_phenotype_without_labels.tsv"), sep="\t", index=False)
# S10: compliance of honeycomb units (hc_C = 0.60 in the model) - first-order post hoc correction of FVC to hc_C = 0.15
P = pd.concat([CAL, ADV], ignore_index=True); rows = []
for lab, shift in (("Model (honeycomb compliance 0.60)", 0.0), ("Honeycomb compliance 0.15 (post hoc approximation)", 0.45)):
    D = P.copy(); D["FVC"] = D.FVC - 100 * D.hc * shift; a = D[D.FVC < 50]
    Cm = CAL.copy(); Cm["FVC"] = Cm.FVC - 100 * Cm.hc * shift; Mm = Cm[(Cm.FVC >= 50) & (Cm.FVC <= 90)].reset_index(drop=True); pm = point(Mm)
    rows.append(dict(scenario=lab, FVC_min=round(D.FVC.min(), 1), n_FVC_lt50=len(a), PH_FVC_lt50_pct=round(100 * (a.mPAP >= 25).mean(), 1),
                     severe_FVC_lt50_pct=round(100 * (a.mPAP > 40).mean(), 1), n_mild=len(Mm), **{f"mild_{k}": round(v, 3) for k, v in pm.items()}))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "TableS10_honeycomb_compliance.tsv"), sep="\t", index=False)
# S11: mPAP categories with measurement noise added to mPAP (SD 2.5 mmHg), mean of 1,000 realisations
g = np.random.default_rng(0); Pn = []
for _ in range(1000):
    mp = MILD.mPAP + g.normal(0, 2.5, len(MILD)); Pn.append([100 * (mp < 20).mean(), 100 * ((mp >= 20) & (mp < 25)).mean(), 100 * (mp >= 25).mean(), 100 * (mp >= 35).mean()])
pd.DataFrame([dict(version="Without noise (manuscript)", lt20=100 * (MILD.mPAP < 20).mean(), b20_25=100 * ((MILD.mPAP >= 20) & (MILD.mPAP < 25)).mean(),
                   ge25=100 * (MILD.mPAP >= 25).mean(), ge35=100 * (MILD.mPAP >= 35).mean()),
              dict(version="With measurement noise", **dict(zip(["lt20", "b20_25", "ge25", "ge35"], np.mean(Pn, 0))))]).round(1).to_csv(
    os.path.join(OUT, "TableS11_categories_with_noise.tsv"), sep="\t", index=False)
print("tables S8-S11 written")
