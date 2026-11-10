"""Provenance and the shared readers: cohorts match their manifest; one definition of the final visit and the lesion."""
import os, sys
import numpy as np, pandas as pd
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import cohort_manifest
from cohort_io import last_visit, usable, applied_severity, effective_lesion, MANUSCRIPT


def test_cohorts_match_the_manifest():
    M = pd.read_csv(cohort_manifest.OUT, sep="\t")
    assert set(M.cohort) == set(cohort_manifest.SPEC) == {f[:-7] for f in os.listdir(cohort_manifest.F) if f.endswith(".tsv.gz")}
    for r in M.itertuples():
        assert cohort_manifest.md5(os.path.join(cohort_manifest.F, r.cohort + ".tsv.gz")) == r.md5, r.cohort


def test_last_visit_is_the_last_viable_one():
    d = pd.DataFrame(dict(mPAP=[14.0, 20.0, 30.0, 80.0], gas_ok=[True, True, True, False], Q=[5.0, 5.0, 4.0, 2.0], converged=[1, 1, 0, 1]))
    r0, x = last_visit(d)
    assert r0.mPAP == 14.0 and x.mPAP == 30.0
    assert usable(d).tolist() == [True, True, False, False]


def test_one_definition_of_the_lesion():
    d = pd.DataFrame(dict(vasc_A=[2.0, 2.0], t_vasc0=[1.0, 1.0], years=[0.0, 3.0], a_scale=[3.0, 3.0], tau=[4.0, 4.0]))
    assert applied_severity(d) == 6.0
    assert np.allclose(effective_lesion(d), 6.0 * (1 - np.exp(-(np.array([0.0, 3.0]) + 1.0) / 4.0)))
    assert MANUSCRIPT["tau"] == 2.0 and MANUSCRIPT["p_vasc"] == 0.7
