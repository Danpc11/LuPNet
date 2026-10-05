#!/usr/bin/env python3
"""Venous test of the shared target law (Huang et al. 1996 pulmonary veins, 15 orders). Two tests:
(1) free scale tau0 and profiled b, as for the arteries; (2) out-of-sample: the arterial tau0 and b (fitted on the
arterial tree) applied unchanged to the veins, for several venous pulsatilities phi1_v."""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import brentq, minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lupnet.lung import womersley_G as G_tube
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
Q = 5.0 / 60 * 1e-3; RHO = 1060.0; HR = 70.0; W1 = 2.865; w = 2 * np.pi * HR / 60


def tree(fname):
    H = pd.read_csv(os.path.join(ROOT, 'data', fname), sep='\t')
    return H, Q / H.N.values, H.L_cm.values / 100, H.mu_cP.values * 1e-3, H.D_cm.values / 100


def radius(Qv, L, mu, b, tau0, phi):
    f = lambda x: np.log(4 * mu * Qv / (np.pi * np.exp(3 * x)) * (1 + W1 * phi * G_tube(np.exp(x) * np.sqrt(w * RHO / mu))[0])) - \
        np.log(tau0 * np.exp(x) ** (b - 1) * L ** ((b - 1) / 2))
    return float(np.exp(brentq(f, np.log(1e-9), np.log(5.0))))


def predict(T, b, tau0, phi):
    _, Qn, Ln, mun, _ = T
    return np.array([2 * radius(q, l, m, b, tau0, phi) for q, l, m in zip(Qn, Ln, mun)])


def err(T, D):
    return float(np.sqrt(np.mean((np.log10(D) - np.log10(T[4])) ** 2)))


def fit(T, b, phi):
    r = minimize_scalar(lambda lt: err(T, predict(T, b, 10 ** lt, phi)) ** 2, bounds=(-3, 2), method='bounded', options=dict(xatol=1e-4))
    return float(np.sqrt(r.fun)), 10 ** r.x


A = tree('huang1996_pulmonary_arteries.tsv'); V = tree('huang1996_pulmonary_veins.tsv')
bgrid = np.round(np.arange(0.50, 1.201, 0.025), 3)
rows = []
_, tau0_art = fit(A, 0.90, 0.8)
print(f"arterial fit: b = 0.90, phi1 = 0.8, tau0 = {tau0_art:.3f} Pa")
for phi in (0.0, 0.4, 0.8):
    res = [fit(V, b, phi) for b in bgrid]; rms = np.array([r[0] for r in res]); ib = int(rms.argmin())
    e_oos = err(V, predict(V, 0.90, tau0_art, phi))
    rows.append(dict(phi1_v=phi, best_b_veins=bgrid[ib], factor_best=round(10 ** rms[ib], 3), tau0_veins=round(res[ib][1], 3),
                     factor_out_of_sample_arterial_law=round(10 ** e_oos, 3)))
    print(rows[-1], flush=True)
R = pd.DataFrame(rows); R.to_csv(os.path.join(ROOT, 'results', 'pulmonary_veins_shared_law.tsv'), sep='\t', index=False)
D = predict(V, 0.90, tau0_art, 0.4)
print("order  D_meas(mm)  D_pred_oos(mm)  [arterial tau0 and b = 0.90, phi1_v = 0.4]")
for o, dm, dp in zip(V[0].order, V[4], D): print(f"  {o:2d}  {dm*1e3:8.3f}  {dp*1e3:8.3f}")
