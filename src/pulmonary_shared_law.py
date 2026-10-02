#!/usr/bin/env python3
"""Step 1-2 of the lung programme: does the shared target law of InFlow predict the human pulmonary arterial tree?

Data: Huang et al. 1996 (J Appl Physiol 81:2123), 16 diameter-defined Strahler orders of the human pulmonary arterial
tree (number of elements N_n, diameter D_n, length L_n, apparent viscosity mu_n), as tabulated by Shi (2010).
Each order carries Q_n = Q / N_n (Q = 5 L/min). Model: target tau* = tau0 r^(b-1) l^((b-1)/2) (Eq. 1 of InFlow), sensed
shear S = tau_mean [1 + w1 phi1 G(alpha)] with the harmonic law (w1 = 2.865 from Feaver 2013; phi1 = first-harmonic
amplitude of pulmonary flow relative to its mean). For each order the radius is the root of S(r) = tau*(r, L_n).
tau0 is the single free scale (one constant for all 16 orders), fitted by least squares in log D; b is profiled.
Variants: phi1 = 1.0 / 1.6 / 1.8, pulsatility undamped or damped along the tree, mean-only sensing, Murray (b = 1).
"""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import brentq, minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lupnet.lung import womersley_G as G_tube

H = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'huang1996_pulmonary_arteries.tsv'), sep='\t')
Q = 5.0 / 60 * 1e-3; RHO = 1060.0; HR = 70.0; W1 = 2.865
Qn = Q / H.N.values; Ln = H.L_cm.values / 100; mun = H.mu_cP.values * 1e-3; Dm = H.D_cm.values / 100


def radius(Qv, L, mu, b, tau0, phi):
    w = 2 * np.pi * HR / 60
    f = lambda x: np.log(4 * mu * Qv / (np.pi * np.exp(3 * x)) * (1 + W1 * phi * G_tube(np.exp(x) * np.sqrt(w * RHO / mu))[0])) - \
        np.log(tau0 * np.exp(x) ** (b - 1) * L ** ((b - 1) / 2))
    return float(np.exp(brentq(f, np.log(1e-9), np.log(5.0))))


def predict(b, tau0, phi_main, damp):
    phis = phi_main * (np.exp(-(16 - H.order.values) / damp) if damp else np.ones(len(H)))
    return np.array([2 * radius(q, l, m, b, tau0, p) for q, l, m, p in zip(Qn, Ln, mun, phis)])


def fit(b, phi, damp):
    obj = lambda lt: np.mean((np.log10(predict(b, 10 ** lt, phi, damp)) - np.log10(Dm)) ** 2)
    r = minimize_scalar(obj, bounds=(-3, 2), method='bounded', options=dict(xatol=1e-4))
    return float(np.sqrt(r.fun)), 10 ** r.x


rows = []
bgrid = np.round(np.arange(0.40, 1.101, 0.025), 3)
for lab, phi, damp in (("pulsatile phi1=1.6, undamped", 1.6, 0), ("pulsatile phi1=1.6, damped (e-fold 5 orders)", 1.6, 5),
                       ("pulsatile phi1=1.0", 1.0, 0), ("pulsatile phi1=1.8", 1.8, 0), ("mean-only", 0.0, 0)):
    res = [fit(b, phi, damp) for b in bgrid]
    rms = np.array([r[0] for r in res]); ib = int(rms.argmin())
    i1 = int(np.where(np.isclose(bgrid, 1.0))[0][0]); i675 = int(np.where(np.isclose(bgrid, 0.675))[0][0])
    rows.append(dict(model=lab, best_b=bgrid[ib], rms_log10D_best=round(rms[ib], 4), factor_best=round(10 ** rms[ib], 3),
                     rms_at_b0675=round(rms[i675], 4), rms_at_b1=round(rms[i1], 4)))
    print(rows[-1], flush=True)
R = pd.DataFrame(rows)
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'), exist_ok=True); R.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'pulmonary_shared_law.tsv'), sep='\t', index=False)
# per-order comparison at the main model, b = 0.675 and best b
for b in (0.675, float(R.iloc[0].best_b), 1.0):
    _, t0 = fit(b, 1.6 if b != 1.0 else 0.0, 0)
    D = predict(b, t0, 1.6 if b != 1.0 else 0.0, 0)
    print(f"b={b}: order  D_meas(mm)  D_pred(mm)")
    for o, dm, dp in zip(H.order, Dm, D): print(f"   {o:2d}  {dm*1e3:8.3f}  {dp*1e3:8.3f}")
