import json, os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from lupnet.lung import Lung, womersley_G
from lupnet.params import Params


def test_womersley_limits():
    g = womersley_G(np.array([1e-4, 0.5, 10.0]))
    assert abs(g[0] - 1) < 1e-9 and abs(g[1] - 1) < 1e-2 and abs(g[2] - 2.78) < 0.02


def test_healthy_lung_shared_law():
    L = Lung(Params(law="shared", b=0.925, G=8))
    vs = L.base_state["vs"]
    assert 12.0 < vs["PPA"] < 17.5          # normal mean pulmonary artery pressure
    assert L.base_state["zerr"] < 1e-4      # the healthy tree is a rest point of the remodelling rule


def test_vasculopathy_raises_resistance():
    L = Lung(Params(law="shared", b=0.925, G=8))
    p0 = L.base_state["vs"]["PPA"]
    L.R_lesion = np.full(L.N, 2.0 * L._Rta0)
    p1 = L.equilibrate()["vs"]["PPA"]
    assert p1 > p0 + 1.0


def test_calculator_coefficients():
    C = json.load(open(os.path.join(ROOT, "results", "ph_calculator.json")))
    assert C["logistic"]["coef_log_ratio"] > 0     # a disproportionately low DLCO raises the probability of PH
    assert C["mixture"]["means"][C["mixture"]["vascular_component"]] < min(C["mixture"]["means"]) + 1e-12
