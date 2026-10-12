#!/usr/bin/env python3
"""Non-invasive identification of vasculopathy and PH: echocardiographic systolic pressure and arterial blood gases.

1. Clinical observables from every point of the simulated trajectories: FVC and DLCO (% of healthy baseline), resting PaO2,
   alveolar-arterial O2 gradient (A-a, sea level: PAO2 = 150 - PaCO2/0.8), and echocardiographic systolic pulmonary artery
   pressure (sPAP, from mPAP = 0.61 sPAP + 2 mmHg).
2. Emulator: gradient-boosted trees mapping the latent disease state to the observables (validated by case).
3. Detection: logistic models for relevant vasculopathy (effective lesion >= baseline arteriolar resistance) and PH,
   with measurement noise; areas under the ROC curve by case-grouped cross-validation (table S13).
4. Individual inference: for held-out simulated cases, a posterior over the latent state given the observables, from
   the emulator, gives each case's mPAP with a 90% interval and its probability of relevant vasculopathy (table S14).
Writes Figure 5 and tables S13-S14 to results/figures/.
"""
import os, sys, logging
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg"); logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, roc_curve
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import cases

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
F = os.path.join(ROOT, "results", "ph_cohorts"); OUT = os.path.join(ROOT, "results", "figures"); os.makedirs(OUT, exist_ok=True)
NOISE = dict(FVC=5.0, DLCO=5.0, PaO2=5.0, AaDO2=6.0, sPAP=8.0, SpO2=2.0)       # measurement SDs (percentage points; mmHg)
LAT = ["healthy", "fibrotic", "honeycomb", "years", "Aeff", "eta", "hpv", "lost_base", "lost_apex"]


def trajectory_points():
    rows = []
    for tag in ("big2", "adv2", "hom", "het"):
        for i, d in cases(F, tag):
            r0 = d.iloc[0]; A = float(d.vasc_A.iloc[0]); t0 = float(d.t_vasc0.iloc[0]) if "t_vasc0" in d else 0.0
            for _, x in d.iterrows():
                if not x.converged:
                    continue
                rows.append(dict(case=f"{tag}{i}", healthy=x.healthy, fibrotic=x.fibrotic, honeycomb=x.honeycomb, years=x.years,
                                 Aeff=A * (1 - np.exp(-(x.years + t0) / 2.0)) if A > 0 else 0.0,
                                 eta=float(x.get("eta_p", 0.0)) if tag == "het" else 0.0, hpv=float(x.get("hpv_max", 3.0)) if tag == "het" else 3.0,
                                 lost_base=x.lost_base, lost_apex=x.lost_apex, mPAP=x.mPAP, FVC=100 * x.FVC / r0.FVC,
                                 DLCO=100 * x.DLCO / r0.DLCO, PaO2=x.PaO2, PaCO2=x.PaCO2, SpO2=100 * x.SaO2))
    D = pd.DataFrame(rows).replace([np.inf, -np.inf], np.nan).dropna(subset=["mPAP", "FVC", "DLCO", "PaO2", "PaCO2"]).reset_index(drop=True)
    D["AaDO2"] = 150 - D.PaCO2 / 0.8 - D.PaO2
    D["sPAP"] = (D.mPAP - 2) / 0.61
    return D


def gbr():
    return HistGradientBoostingRegressor(max_iter=400, learning_rate=0.06, max_leaf_nodes=31, min_samples_leaf=20, random_state=0)


D = trajectory_points()
# ---------------------------------------------------------------- emulator validation (held out by case)
emu_rows = []
for y in ["mPAP", "FVC", "DLCO", "PaO2", "AaDO2", "SpO2"]:
    pred = np.zeros(len(D))
    for tr, te in GroupKFold(5).split(D, groups=D.case):
        pred[te] = gbr().fit(D[LAT].iloc[tr], D[y].iloc[tr]).predict(D[LAT].iloc[te])
    e = pred - D[y].values
    emu_rows.append(dict(output=y, R2=round(1 - np.var(e) / np.var(D[y]), 3), MAE=round(float(np.mean(np.abs(e))), 2)))
pd.DataFrame(emu_rows).to_csv(os.path.join(OUT, "TableS13a_emulator.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- detection by case-grouped cross-validation
V = D[D.case.str.match(r"^(big2|adv2)") & (D.healthy >= 0.12) & (D.healthy <= 0.85)].reset_index(drop=True)
V["vasc"] = (V.Aeff >= 1.0).astype(int); V["PH"] = (V.mPAP >= 25).astype(int)
SETS = {"FVC + DLCO": ["FVC", "DLCO"], "+ SpO$_2$": ["FVC", "DLCO", "SpO2"], "+ PaO$_2$": ["FVC", "DLCO", "PaO2"], "+ A–a gradient": ["FVC", "DLCO", "AaDO2"],
        "+ echo sPAP": ["FVC", "DLCO", "sPAP"], "+ echo sPAP + A–a gradient": ["FVC", "DLCO", "sPAP", "AaDO2"]}


def noisy(df, cols, g):
    X = pd.DataFrame({c: df[c] + g.normal(0, NOISE[c], len(df)) for c in cols}, index=df.index)
    X["FVC_DLCO"] = X.FVC / X.DLCO.clip(lower=3)
    return X


det, rocs = [], {}
for target in ("vasc", "PH"):
    for name, cols in SETS.items():
        aucs = []
        for rep in range(5):
            g = np.random.default_rng(rep); X = noisy(V, cols, g); y = V[target].values; pr = np.zeros(len(V))
            for tr, te in GroupKFold(5).split(X, y, groups=V.case):
                pr[te] = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(X.iloc[tr], y[tr]).predict_proba(X.iloc[te])[:, 1]
            aucs.append(roc_auc_score(y, pr))
            if rep == 0:
                rocs[(target, name)] = roc_curve(y, pr)
        det.append(dict(target=target, observations=name, AUC=round(float(np.mean(aucs)), 3), AUC_min=round(min(aucs), 3), AUC_max=round(max(aucs), 3)))
DET = pd.DataFrame(det); DET.to_csv(os.path.join(OUT, "TableS13_detection.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- individual inference on held-out simulated cases
rng = np.random.default_rng(0)
pool = [c for c in D.case.unique() if c.startswith(("big2", "adv2"))]
test = set(rng.choice(pool, 100, replace=False)); TR = D[~D.case.isin(test)]
EM = {y: gbr().fit(TR[LAT], TR[y]) for y in ["FVC", "DLCO", "PaO2", "AaDO2", "SpO2", "mPAP"]}
emu_sd = {}
for y in EM:                                                    # emulator error, added to the measurement error
    pred = np.zeros(len(TR))
    for tr, te in GroupKFold(5).split(TR, groups=TR.case):
        pred[te] = gbr().fit(TR[LAT].iloc[tr], TR[y].iloc[tr]).predict(TR[LAT].iloc[te])
    emu_sd[y] = float(np.std(pred - TR[y].values))
emu_sd["sPAP"] = emu_sd["mPAP"] / 0.61
# prior: disease states of the training trajectories, with lesion severity and onset re-drawn from their priors
P = TR[TR.case.str.match(r"^(big2|adv2)") & (TR.healthy >= 0.1)][LAT].sample(200000, replace=True, random_state=1).reset_index(drop=True)
g = np.random.default_rng(2); A = np.exp(np.log(1.5) + 1.1 * g.normal(size=len(P))) * (g.uniform(size=len(P)) < 0.7)
P["Aeff"] = A * (1 - np.exp(-(P.years.values + g.uniform(0, 10, len(P))) / 2))
PR = {y: EM[y].predict(P[LAT]) for y in EM}; PR["sPAP"] = (PR["mPAP"] - 2) / 0.61
T = D[D.case.isin(test) & (D.healthy >= 0.12) & (D.healthy <= 0.85)].groupby("case").sample(3, random_state=3, replace=True).reset_index(drop=True)
INF, keep = [], {}
for name, cols in SETS.items():
    gg = np.random.default_rng(5); est, lo, hi, pv = [], [], [], []
    for _, row in T.iterrows():
        ll = np.zeros(len(P))
        for o in cols:
            sd = np.hypot(NOISE[o], emu_sd[o]); ll += -0.5 * ((PR[o] - (row[o] + gg.normal(0, NOISE[o]))) / sd) ** 2
        w = np.exp(ll - ll.max()); w /= w.sum(); order = np.argsort(PR["mPAP"]); cw = np.cumsum(w[order])
        est.append(float((w * PR["mPAP"]).sum())); lo.append(PR["mPAP"][order][np.searchsorted(cw, 0.05)]); hi.append(PR["mPAP"][order][np.searchsorted(cw, 0.95)])
        pv.append(float((w * (P.Aeff.values >= 1)).sum()))
    est, lo, hi, pv = map(np.array, (est, lo, hi, pv)); tm = T.mPAP.values
    keep[name] = (est, lo, hi)
    INF.append(dict(observations=name, n=len(T), MAE_mPAP=round(float(np.mean(np.abs(est - tm))), 2), interval90_width=round(float(np.mean(hi - lo)), 1),
                    coverage90=round(float(np.mean((tm >= lo) & (tm <= hi))), 2), AUC_vasculopathy=round(roc_auc_score(T.Aeff >= 1, pv), 3),
                    AUC_PH=round(roc_auc_score(T.mPAP >= 25, est), 3)))
INFD = pd.DataFrame(INF); INFD.to_csv(os.path.join(OUT, "TableS14_individual_inference.tsv"), sep="\t", index=False)
print(pd.DataFrame(emu_rows).to_string(index=False)); print(DET.to_string(index=False)); print(INFD.to_string(index=False))

# ---------------------------------------------------------------- Figure 5
MM = 1 / 25.4
plt.rcParams.update({"font.family": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False, "axes.linewidth": 0.6, "savefig.dpi": 600, "pdf.fonttype": 42})
BLUE, RED, GREY, GREEN, DARK = "#2166AC", "#B2182B", "#9E9E9E", "#1B7837", "#252525"
fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 66 * MM), gridspec_kw=dict(wspace=0.45))
ax = axs[0]; names = list(SETS); yv = np.arange(len(names))
for k, (tg, col, lab) in enumerate((("vasc", RED, "Vasculopathy"), ("PH", BLUE, "PH"))):
    s = DET[DET.target == tg].set_index("observations").loc[names]
    ax.barh(yv + (0.2 if k == 0 else -0.2), s.AUC, height=0.38, color=col, label=lab)
ax.set_yticks(yv); ax.set_yticklabels(names); ax.invert_yaxis(); ax.set_xlim(0.5, 1.0); ax.set_xlabel("Area under the ROC curve")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6); ax.text(-0.75, 1.03, "a", transform=ax.transAxes, fontsize=9, fontweight="bold")
ax = axs[1]
for name, col in (("FVC + DLCO", GREY), ("+ echo sPAP + A–a gradient", BLUE)):
    fpr, tpr, _ = rocs[("vasc", name)]; ax.plot(fpr, tpr, color=col, label=f"{name}, AUC {DET[(DET.target == 'vasc') & (DET.observations == name)].AUC.iloc[0]:.2f}")
ax.plot([0, 1], [0, 1], color="#DDDDDD", lw=0.6); ax.set(xlabel="1 − specificity", ylabel="Sensitivity", xlim=(0, 1), ylim=(0, 1))
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=5.6, handlelength=1.2); ax.set_title("Relevant vasculopathy", fontsize=7)
ax.text(-0.25, 1.03, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")
ax = axs[2]; first = ~T.case.duplicated().values; tm = T.mPAP.values[first]
est, lo, hi = (a[first] for a in keep["FVC + DLCO"])
ax.plot(tm, est, "o", ms=2.6, color=GREY, alpha=0.8, label=f"FVC + DLCO (error {INFD.set_index('observations').loc['FVC + DLCO', 'MAE_mPAP']:.1f} mmHg)")
est, lo, hi = (a[first] for a in keep["+ echo sPAP + A–a gradient"])
ax.errorbar(tm, est, yerr=[est - lo, hi - est], fmt="o", ms=2.6, color=BLUE, ecolor=BLUE, elinewidth=0.5, alpha=0.85, capsize=0,
            label=f"+ echo sPAP + A–a gradient (error {INFD.set_index('observations').loc['+ echo sPAP + A–a gradient', 'MAE_mPAP']:.1f} mmHg; 90% interval)")
ax.plot([12, 45], [12, 45], color=DARK, lw=0.6, ls=":"); ax.set(xlabel="Simulated mPAP (mmHg)", ylabel="Estimated mPAP (mmHg)", xlim=(12, 45), ylim=(10, 50))
ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.22), fontsize=5.6, handlelength=1.2); ax.text(-0.25, 1.03, "c", transform=ax.transAxes, fontsize=9, fontweight="bold")
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"Figure5.{ext}"), bbox_inches="tight", pad_inches=0.04)
fig.savefig(os.path.join(OUT, "Figure5.tif"), bbox_inches="tight", pad_inches=0.04, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
print("Figure 5 and tables S13-S14 written")
