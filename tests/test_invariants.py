"""Mathematical invariants of the lung model (added by the audit)."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from lupnet.lung import Lung, LPM, MMHG, pois
from lupnet.params import Params

L = Lung(Params(law="shared", b=0.925, G=8)); B = L.base_state; vs = B["vs"]; p = L.p


def test_flow_conservation():
    # every unit's flow adds up to the cardiac output, and each bifurcation conserves flow
    assert abs(vs["Ql"].sum() - B["Q"]) < 1e-9 * B["Q"]
    for k in L.inner_levels:
        assert np.allclose(vs["Q"][k], vs["Q"][L.c1[k]] + vs["Q"][L.c2[k]], rtol=1e-9)


def test_pressure_identity():
    # mPAP = LAP + Q x equivalent resistance, and PVR (Wood units) = (mPAP - LAP) / Q
    assert abs(vs["PPA"] - (p.LAP + B["Q"] * vs["Req"])) < 1e-9
    assert abs(vs["dP"] / B["Q"] - vs["Req"]) < 1e-9


def test_poiseuille_units():
    # 1 mmHg min/L = 1 Wood unit = MMHG / LPM Pa s m^-3
    R = pois(3.5e-3, 0.01, 1e-3) * LPM / MMHG
    assert abs(R - 8 * 3.5e-3 * 0.01 / (np.pi * 1e-12) * LPM / MMHG) < 1e-12


def test_shared_law_unique_root():
    # the mismatch between sensed shear and target decreases monotonically with radius, so the root is unique
    f, Lseg = p.Q * LPM / 2 ** 5, p.L_pa * 2 ** (-5 / 3)
    tau0 = L.tau0_shared
    r = np.exp(np.linspace(np.log(1e-6), np.log(0.05), 400))
    S = 4 * p.mu_b * f / (np.pi * r ** 3) * L.pulse_factor(r, True)
    g = np.log(S) - np.log(tau0 * r ** (p.b - 1) * Lseg ** ((p.b - 1) / 2))
    assert np.all(np.diff(g) < 0)


def test_healthy_rest_point_and_normalisation():
    assert B["zerr"] < 1e-4
    assert abs(B["DL"] / L.DL0 - 1) < 1e-12
