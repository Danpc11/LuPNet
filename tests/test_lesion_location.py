"""Where the vascular lesion sits (v0.3.13): the default must be the diffuse lesion of the manuscript."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from lupnet import Params
from lupnet.lung import Lung
from lupnet.sim import lesion_weights, apply_lesion


def lung():
    L = Lung(Params(law="shared", b=0.925, G=6)); L.state[:20] = 1; L.state[20:24] = 2
    return L


def test_weights():
    L = lung(); scar = L.state > 0
    assert np.all(lesion_weights(L, "diffuse") == 1)
    w = lesion_weights(L, "scar"); assert np.all(w[scar] == 1) and np.all(w[~scar] == 0)
    f = lesion_weights(L, "field"); assert np.all(f[scar] == 1) and f[~scar].max() <= 1 and (f[~scar] > 0).any() and (f[~scar] == 0).any()
    m = lesion_weights(L, "mixed", 0.25); assert np.allclose(m, 0.25 + 0.75 * f)


def test_default_is_diffuse_and_scar_lesion_adds_little():
    L = lung(); assert L.p.vasc_where == "diffuse"
    p0 = L.equilibrate()["vs"]["PPA"]
    out = {}
    for where in ("scar", "field", "diffuse"):
        L.p.vasc_where = where; apply_lesion(L, 4.0); out[where] = L.equilibrate()["vs"]["PPA"]
    assert out["scar"] - p0 < 0.5 * (out["diffuse"] - p0)      # blood bypasses the scar: a lesion there adds little
    assert out["scar"] <= out["field"] <= out["diffuse"]
