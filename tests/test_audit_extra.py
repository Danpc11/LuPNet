"""Guards for the findings of the v0.3.7 audit (see docs/AUDIT.md, findings 9-13)."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from cohort_io import cases
import audit_extra


def test_years_column_is_t_times_years_per_unit():
    # Figure 2f of v0.3.7 rebuilt the lesion clock from t x ypu and lost the 6.71 years per model unit;
    # the stored `years` column is the clock to use
    for i, d in cases(os.path.join(ROOT, "results", "ph_cohorts"), "big2"):
        if i > 5:
            break
        ratio = (d.years / (d.t * d.ypu)).iloc[1:]
        assert np.allclose(ratio, 6.71, rtol=1e-6)


def test_aged_baseline_without_fibrosis():
    A = audit_extra.aged_baseline()
    assert 14.0 < A["mPAP_healthy"] < 15.0
    assert 16.0 < A["mPAP"] < 17.0           # ageing alone adds about 2.2 mmHg
    assert A["zerr"] < 1e-4


def test_pressure_decomposition_adds_up():
    P = audit_extra.pressure_decomposition()
    assert abs(P["total"] - (P["mPAP"] - 8.0)) < 0.05
    assert 40 < P["set_by_law_pct"] < 55     # about half of the healthy drop comes from the law


def test_lesion_saturated_at_visit():
    S = audit_extra.lesion_saturation(("big2",))
    assert S.saturation.min() > 0.5 and S.saturation.median() > 0.95
