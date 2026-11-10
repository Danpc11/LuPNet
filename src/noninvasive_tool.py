#!/usr/bin/env python3
"""Individual, non-invasive estimate of pulmonary hypertension in fibrosing interstitial lung disease.

Inputs: FVC and DLCO (% predicted) are required; echocardiographic systolic pressure, pulmonary artery diameter on CT,
BNP, right ventricular fractional area change and right atrial area are optional, in any combination.
Right heart catheterisation is never an input: it is used only to fit the measurement models and to score the tool.

Structure (Bayes):
    p(state | data)  ∝  prior_LuPNet(state)  x  p(FVC, DLCO | state)  x  prod_k p(x_k | mPAP(state))
  - prior: the simulated visits of the manuscript cohorts (big2 and adv2, established disease), each with its lung
    lost, effective lesion, FVC, DLCO, mPAP and PVR. The share of visits with PH is re-weighted to the prevalence of
    the population the tool is used in (--prevalence, or the training folds in cross-validation).
  - FVC and DLCO: measurement noise (5 points) plus the inter-individual spread of the healthy baseline (13% and 15%,
    supplement table S8), because the model gives % of the case's own baseline and patients give % predicted.
  - optional inputs: x_k = a_k + b_k mPAP + noise, fitted on catheterised patients (training folds only).
What the model adds to the estimate: mPAP with an interval, the probability of a relevant vasculopathy, and the split
of the pressure into what the lost lung explains and what the lesion adds, from the partial dependence surface
mPAP = F(lung lost, effective lesion) of the simulations (R2 0.998).

    python3 src/noninvasive_tool.py --clinical patients.xlsx        # cross-validation and tables
Columns of the clinical file are mapped in COLS; the file is not part of the repository.
"""
import argparse, os, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import cases, usable, effective_lesion

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results", "figures")
SYN_OUT = os.path.join(ROOT, "results", "synthetic")      # outputs of a synthetic clinical file never share a folder with real ones
from ph_stats import SD_PFT, SPREAD
OPTIONAL = ["sPAP_echo", "PA_diam", "logBNP", "RV_FAC", "RA_area"]
COLS = {"FVC": "FVC (%)", "DLCO": "DLCO aj (%)", "mPAP": "PmAP", "PVR": "RVP (UW)", "PAWP": "Cuña", "sPAP_echo": "PSAP",
        "PA_diam": "Diámetro AP (mm)", "BNP": "BNP (pg/mL)", "RV_FAC": "FACVD (%)", "RA_area": "ÁREA_AD(CM2)",
        "pattern": "Patrón tomográfico", "dx": "Diagnóstico", "ID": "ID"}


# ------------------------------------------------------------------ reference (prior) and partial dependence surface
# reference cohorts of each tomographic pattern; "other" and any pattern without its own cohort use the calibrated set
REFERENCE_TAGS = {"UIP": ("big2_uip", "adv2_uip"), "NSIP": ("big2_nsip", "adv2_nsip"), "other": ("big2", "adv2")}


def trajectory_points(tags=("big2", "adv2", "hom", "big2_uip", "big2_nsip", "adv2_uip", "adv2_nsip")):
    rows = []
    for tag in tags:
        for i, d in cases(os.path.join(ROOT, "results", "ph_cohorts"), tag):
            r0 = d.iloc[0]; d = d.assign(_Aeff=effective_lesion(d))                         # severity as applied
            for _, x in d.iloc[1:][usable(d).iloc[1:].values].iterrows():
                rows.append(dict(case=f"{tag}{i}", tag=tag, lost=1 - x.healthy, Aeff=float(x._Aeff),
                                 mPAP=x.mPAP, PVR=x.PVR, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO, SpO2=100 * x.SaO2))
    return pd.DataFrame(rows)


def surface(T):
    """mPAP = F(lung lost, effective lesion): gradient-boosted fit to every simulated visit."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    em = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.06, min_samples_leaf=20, random_state=0).fit(T[["lost", "Aeff"]], T.mPAP)
    return em, float(em.score(T[["lost", "Aeff"]], T.mPAP))


def reference(T, per_case=25, tags=("big2", "adv2")):
    """Visits in established disease of the cohorts with vasculopathy prior (as the published calculator)."""
    R = T[T.tag.isin(list(tags)) & (T.lost >= 0.15) & (T.lost <= 0.88)]
    keep = []
    for _, g in R.groupby("case", sort=False):
        keep.extend(g.index[np.linspace(0, len(g) - 1, min(per_case, len(g))).round().astype(int)])
    return R.loc[keep].reset_index(drop=True)


# ------------------------------------------------------------------ measurement models and posterior
def fit_obs(train):
    """x = a + b mPAP + N(0, sd), per optional input, from catheterised patients."""
    M = {}
    for k in OPTIONAL:
        s = train[k].notna() & train.mPAP.notna()
        if s.sum() < 8:
            continue
        b, a = np.polyfit(train.mPAP[s], train[k][s], 1)
        M[k] = (a, b, float(np.std(train[k][s] - (a + b * train.mPAP[s]), ddof=2)))
    return M


EDGES = np.array([-np.inf, 20, 25, 35, 45, np.inf])       # mPAP classes of the severity prior, mmHg


def severity_weights(R, mpap_population, pseudo=0.5):
    """Prior weights that give the reference visits the distribution of mPAP of the population the tool is used in
    (five classes; `pseudo` counts per class keep an empty class from being excluded)."""
    cp = np.histogram(mpap_population, EDGES)[0] + pseudo; cp = cp / cp.sum()
    cr = np.histogram(R.mPAP.values, EDGES)[0] / len(R)
    k = np.digitize(R.mPAP.values, EDGES[1:-1])
    return cp[k] / np.maximum(cr[k], 1e-9)


def posterior(R, row, use, M, prevalence=None, prior_w=None):
    """Weights of the reference visits given a patient's data. `use`: optional inputs to include."""
    ll = np.zeros(len(R))
    for k in ("FVC", "DLCO"):
        sd = np.hypot(SD_PFT, SPREAD[k] * R[k].values)
        ll += -0.5 * ((R[k].values - row[k]) / sd) ** 2 - np.log(sd)
    for k in use:
        if k in M and not pd.isna(row[k]):
            a, b, sd = M[k]; ll += -0.5 * ((a + b * R.mPAP.values - row[k]) / sd) ** 2
    w = np.exp(ll - ll.max())
    if prior_w is not None:
        w = w * prior_w
    elif prevalence is not None:
        ph = R.mPAP.values >= 25; p0 = ph.mean()
        w = w * np.where(ph, prevalence / p0, (1 - prevalence) / (1 - p0))
    return w / w.sum()


def summarise(R, w, F=None):
    o = np.argsort(R.mPAP.values); cw = np.cumsum(w[o]); m = R.mPAP.values
    out = dict(mPAP=float((w * m).sum()), lo=float(m[o][np.searchsorted(cw, 0.05)]), hi=float(m[o][min(np.searchsorted(cw, 0.95), len(m) - 1)]),
               pPH=float((w * (m >= 25)).sum()), pPH20=float((w * (m > 20)).sum()), pVasc=float((w * (R.Aeff.values >= 1)).sum()),
               PVR=float((w * R.PVR.values).sum()), lost=float((w * R.lost.values).sum()))
    if F is not None:
        f0 = F.predict(np.c_[R.lost.values, np.zeros(len(R))]); base = float(F.predict(np.array([[0.05, 0.0]]))[0])
        out.update(baseline=base, from_lung_lost=float((w * (f0 - base)).sum()), from_lesion=float((w * (m - f0)).sum()))
    return out


# ------------------------------------------------------------------ clinical file
def load_clinical(path):
    ext = os.path.splitext(path)[1].lower()
    d = pd.read_csv(path, sep="\t" if ext in (".tsv", ".txt") else ",") if ext in (".tsv", ".txt", ".csv") else pd.read_excel(path)
    d.columns = [str(c).strip() for c in d.columns]
    P = pd.DataFrame({k: d[v] for k, v in COLS.items() if v in d.columns})
    for k in P.columns.difference(["pattern", "dx", "ID"]):
        P[k] = pd.to_numeric(P[k], errors="coerce")
    P["PA_diam"] = P.PA_diam.where(P.PA_diam > 0); P["logBNP"] = np.log(P.BNP)
    pat = P.pattern.astype(str).str.upper()
    P["group"] = np.where(pat.str.contains("NIU"), "UIP", np.where(pat.str.contains("NINE"), "NSIP", "other"))
    return P


def is_synthetic(P):
    """True for the synthetic clinical file (every ID starts with SYN-)."""
    return bool(len(P)) and bool(P.ID.astype(str).str.startswith("SYN-").all())


def output_dir(P):
    d = SYN_OUT if is_synthetic(P) else OUT; os.makedirs(d, exist_ok=True)
    return d


def references_by_pattern(T):
    """One reference per tomographic pattern; falls back to the calibrated set when a pattern cohort is absent."""
    have = set(T.tag.unique()); out = {}
    for g, tags in REFERENCE_TAGS.items():
        present = tuple(t for t in tags if t in have)          # the advanced cohort of a pattern is used when it exists
        out[g] = reference(T, tags=present if present else REFERENCE_TAGS["other"])
    return out


def cross_validate(P, R, sets, n_rep=20, k=5, seed=0, prior="prevalence"):
    """Repeated k-fold: measurement models and prevalence from the training folds; scores on the held-out patients.
    Returns one row per input set and the out-of-fold estimates of every patient (mean over repeats)."""
    from sklearn.metrics import roc_auc_score
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    V = P.dropna(subset=["FVC", "DLCO", "mPAP"]).reset_index(drop=True); y = (V.mPAP >= 25).values
    rng = np.random.default_rng(seed); res, per = [], {}
    for name, use in sets.items():
        est = np.zeros((n_rep, len(V))); pph = np.zeros_like(est); lo = np.zeros_like(est); hi = np.zeros_like(est); plr = np.full_like(est, np.nan)
        cc = V[["FVC", "DLCO"] + use].notna().all(axis=1).values                 # complete cases, for the regression
        for r in range(n_rep):
            fold = rng.permutation(len(V)) % k
            for f in range(k):
                tr, te = np.where(fold != f)[0], np.where(fold == f)[0]
                M = fit_obs(V.iloc[tr]); prev = float(y[tr].mean())
                PW = {g: severity_weights(Rg, V.mPAP.values[tr]) for g, Rg in R.items()} if prior == "severity" else None
                for j in te:
                    g = V.group.iloc[j]; Rj = R[g]                                  # reference of the patient's pattern
                    s = summarise(Rj, posterior(Rj, V.iloc[j], use, M, prev, PW[g] if PW else None))
                    est[r, j], pph[r, j], lo[r, j], hi[r, j] = s["mPAP"], s["pPH"], s["lo"], s["hi"]
                trc = tr[cc[tr]]; tec = te[cc[te]]
                if len(tec) and len(np.unique(y[trc])) == 2:
                    lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000)).fit(V.iloc[trc][["FVC", "DLCO"] + use], y[trc])
                    plr[r, tec] = lr.predict_proba(V.iloc[tec][["FVC", "DLCO"] + use])[:, 1]
        m = V.mPAP.values
        auc = [roc_auc_score(y, pph[r]) for r in range(n_rep)]
        auc_cc = [roc_auc_score(y[cc], pph[r][cc]) for r in range(n_rep)]
        auc_lr = [roc_auc_score(y[cc], plr[r][cc]) for r in range(n_rep)]
        res.append(dict(inputs=name, n=len(V), AUC_LuPNet=np.mean(auc), AUC_min=np.min(auc), AUC_max=np.max(auc),
                        MAE_mPAP=np.mean(np.abs(est - m)), bias=np.mean(est - m), coverage90=np.mean((m >= lo) & (m <= hi)),
                        n_complete=int(cc.sum()), AUC_LuPNet_complete=np.mean(auc_cc), AUC_logistic_complete=np.mean(auc_lr)))
        per[name] = pd.DataFrame(dict(ID=V.ID, group=V.group, mPAP=m, est=est.mean(0), lo=lo.mean(0), hi=hi.mean(0), pPH=pph.mean(0)))
    return pd.DataFrame(res), per, V


SETS = {"FVC + DLCO": [], "+ echo sPAP": ["sPAP_echo"], "+ PA diameter": ["PA_diam"], "+ BNP": ["logBNP"], "+ RV FAC": ["RV_FAC"],
        "+ PA diameter + BNP": ["PA_diam", "logBNP"], "+ all optional": OPTIONAL}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--clinical", required=True); ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--reference", default="pattern", choices=["pattern", "calibrated"],
                    help="pattern: UIP and NSIP patients use the cohort of their pattern; calibrated: every patient uses the manuscript cohorts")
    ap.add_argument("--prior", default="severity", choices=["severity", "prevalence"],
                    help="severity: the reference takes the mPAP distribution of the catheterised population (training folds); "
                         "prevalence: only the share of visits with PH is matched")
    a = ap.parse_args()
    suf = ("" if a.reference == "pattern" else "_calibrated") + ("" if a.prior == "severity" else "_prevalence")
    T = trajectory_points()
    RP = references_by_pattern(T) if a.reference == "pattern" else {g: reference(T) for g in REFERENCE_TAGS}
    have = set(T.tag.unique())
    FP = {g: surface(T[T.tag.isin([t for t in tags if t in have] + (["hom"] if g == "other" else []))])
          for g, tags in REFERENCE_TAGS.items() if any(t in have for t in tags)}
    for g in RP:
        print(f"reference {g}: {len(RP[g])} simulated visits of {RP[g].case.nunique()} cases" + (f"; surface R2 {FP[g][1]:.3f}" if g in FP else ""))
    P = load_clinical(a.clinical); OUT = output_dir(P)
    if is_synthetic(P):
        print("synthetic clinical file: outputs go to results/synthetic/ and are not results")
    S, per, V = cross_validate(P, RP, SETS, n_rep=a.reps, prior=a.prior)
    S.round(3).to_csv(os.path.join(OUT, f"TableS22_noninvasive_tool_cv{suf}.tsv"), sep="\t", index=False); print(S.round(2).to_string(index=False))
    M = fit_obs(V)
    pd.DataFrame([dict(input=k, intercept=v[0], slope_per_mmHg=v[1], sd=v[2]) for k, v in M.items()]).round(3).to_csv(os.path.join(OUT, "TableS23_measurement_models.tsv"), sep="\t", index=False)
    rep = []
    for j, row in V.iterrows():
        g = row.group; F = FP.get(g, FP["other"])[0]; rest = V.drop(index=j)          # models and prior without this patient
        Mj = fit_obs(rest); prevj = float((rest.mPAP >= 25).mean())
        pw = severity_weights(RP[g], rest.mPAP.values) if a.prior == "severity" else None
        s = summarise(RP[g], posterior(RP[g], row, OPTIONAL, Mj, prevj, pw), F)
        rep.append(dict(ID=row.ID, group=g, FVC=row.FVC, DLCO=row.DLCO, mPAP_RHC=row.mPAP, PVR_RHC=row.PVR, **s))
    pd.DataFrame(rep).round(2).to_csv(os.path.join(OUT, f"TableS24_individual_report{suf}.tsv"), sep="\t", index=False)
    for name, df in per.items():
        df.round(2).to_csv(os.path.join(OUT, f"cv{suf}_" + name.replace("+", "plus").replace(" ", "_") + ".tsv"), sep="\t", index=False)
    print(f"tables S22-S24 written to {os.path.relpath(OUT, ROOT)}/")
