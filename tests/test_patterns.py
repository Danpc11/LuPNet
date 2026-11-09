"""Histological patterns (v0.3.9): the default must leave the published cohorts untouched."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from lupnet import Params
from lupnet.lung import Lung
from lupnet.sim import subpleural_factor
from ph_cohort import make_params, PATTERNS


def test_default_pattern_is_neutral():
    L = Lung(Params(law="shared", b=0.925, G=6))
    assert np.all(subpleural_factor(L) == 1.0)
    p = make_params(1, 8, 0.5, 0.0, 0.0, 1.0, 2.0, 0.05)
    assert p.pattern == "calibrated" and p.sub_gain == 0.0 and p.PB == 760.0


def test_uip_seed_is_subpleural_and_basal():
    L = Lung(Params(law="shared", b=0.925, G=6, sub_gain=3.0)); f = subpleural_factor(L)
    assert np.all(f[L.pleural_dist > 1] == 1.0) and f[L.pleural_dist <= 1].min() > 1.0
    ring = L.pleural_dist <= 1
    assert f[ring & (L.y < 0.25)].mean() > f[ring & (L.y > 0.75)].mean()


def test_pattern_parameter_sets():
    u, n = make_params(1, 8, 0.5, 0, 0, 1, 2, 0.05, pattern="uip"), make_params(1, 8, 0.5, 0, 0, 1, 2, 0.05, pattern="nsip")
    assert u.sub_gain > 0 and u.kn > n.kn and u.lam_hc > n.lam_hc and n.at2_sigma < u.at2_sigma
    assert set(PATTERNS) == {"calibrated", "uip", "nsip"}
