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


def test_partial_verification_does_not_replace_the_record(tmp_path, monkeypatch):
    import verify_cohorts
    record = os.path.join(ROOT, "data", "cohorts_verification.tsv"); before = open(record, "rb").read()
    monkeypatch.setattr(verify_cohorts, "WORK", str(tmp_path))
    pd.DataFrame([dict(cohort="hom", case=0, visits_stored=5, visits_new=5, columns=70, max_rel_diff=0.0, column="", mPAP_stored=17.0, mPAP_new=17.0)]
                 ).to_csv(tmp_path / "hom__000.tsv", sep="\t", index=False)
    verify_cohorts.cmd_summary(None)
    assert open(record, "rb").read() == before and (tmp_path / "summary.tsv").exists()
    R = pd.read_csv(record, sep="\t"); M = pd.read_csv(cohort_manifest.OUT, sep="\t")
    assert len(R) == len(M) and (R.identical == R.cases).all() and R.cases.sum() == M.cases.sum()


def test_table_record(tmp_path):
    import table_checksums as tc
    R = pd.read_csv(tc.OUT, sep="\t")
    assert len(R) == 24 and set(R.table[R["how"] == "numeric"]) == set(os.listdir(tc.REF))
    ref = os.path.join(tc.REF, "TableS13_detection.tsv"); D = pd.read_csv(ref, sep="\t")
    num = [c for c in D.columns if pd.api.types.is_numeric_dtype(D[c])]
    a = D.copy(); a[num[-1]] = a[num[-1]] + 0.005; a.to_csv(tmp_path / "a.tsv", sep="\t", index=False)
    b = D.copy(); b[num[-1]] = b[num[-1]] + 0.2; b.to_csv(tmp_path / "b.tsv", sep="\t", index=False)
    assert tc.close(str(tmp_path / "a.tsv"), ref) and not tc.close(str(tmp_path / "b.tsv"), ref)
