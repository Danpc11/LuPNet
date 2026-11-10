"""Pressure-flow response: the flow is imposed, and the lung is left as it was."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import lupnet
from lupnet import Params
from lupnet.lung import Lung
from exercise_response import load_flow


def test_version():
    import json
    assert lupnet.__version__ == json.load(open(os.path.join(ROOT, ".zenodo.json")))["version"]


def test_imposed_flow_and_state_restored():
    L = Lung(Params(law="shared", b=0.925, G=6, rv=True))
    rest = L.equilibrate(); p0 = rest["vs"]["PPA"]
    before = (L.h.copy(), L.ru_a.copy(), L.ru_ta.copy(), L.Ees, L.EDV, L.p.rv)
    m10, q10, _ = load_flow(L, 10.0)
    m75, q75, _ = load_flow(L, 7.5)
    mrv, qrv, _ = load_flow(L, 10.0, imposed=False)
    assert abs(q10 - 10.0) < 1e-6 and abs(q75 - 7.5) < 1e-6          # the flow asked for is the flow obtained
    assert p0 < m75 < m10
    assert qrv < 10.0 - 0.5 and mrv < m10                              # the resting right ventricle does not deliver 10 L/min
    after = (L.h, L.ru_a, L.ru_ta, L.Ees, L.EDV, L.p.rv)
    assert all(np.array_equal(a, b) for a, b in zip(before, after))
    assert abs(L.equilibrate()["vs"]["PPA"] - p0) < 1e-6
