"""A cohort tag must not read the per-case files of another cohort whose tag starts with it."""
import os, sys
import pandas as pd
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from cohort_io import cases


def test_tag_prefix_is_not_confused(tmp_path):
    for name in ("big2_000.tsv", "big2_001.tsv", "big2_uip_000.tsv", "big2_cdmx_000.tsv"):
        pd.DataFrame(dict(x=[1])).to_csv(tmp_path / name, sep="\t", index=False)
    assert [i for i, _ in cases(str(tmp_path), "big2")] == [0, 1]
    assert [i for i, _ in cases(str(tmp_path), "big2_uip")] == [0]


def test_packed_cohort_is_read_when_only_other_tags_have_case_files(tmp_path):
    pd.DataFrame(dict(case=[0, 0, 1], x=[1, 2, 3])).to_csv(tmp_path / "big2.tsv.gz", sep="\t", index=False)
    pd.DataFrame(dict(x=[9])).to_csv(tmp_path / "big2_uip_000.tsv", sep="\t", index=False)
    assert [len(d) for _, d in cases(str(tmp_path), "big2")] == [2, 1]


def test_packed_cohort_wins_over_a_partial_rerun(tmp_path):
    import pytest
    pd.DataFrame(dict(case=[0, 1, 2], x=[1, 2, 3])).to_csv(tmp_path / "hom.tsv.gz", sep="\t", index=False)
    pd.DataFrame(dict(x=[9])).to_csv(tmp_path / "hom_00.tsv", sep="\t", index=False)        # trace of an interrupted re-run
    with pytest.warns(UserWarning, match="using the packed cohort"):
        got = [i for i, _ in cases(str(tmp_path), "hom")]
    assert got == [0, 1, 2]
    for k in (1, 2):                                                                          # a complete per-case cohort is used
        pd.DataFrame(dict(x=[9])).to_csv(tmp_path / f"hom_0{k}.tsv", sep="\t", index=False)
    assert [int(d.x.iloc[0]) for _, d in cases(str(tmp_path), "hom")] == [9, 9, 9]


def test_generating_parameters():
    from cohort_io import meta, lesion_tau, effective_lesion
    import numpy as np
    old = pd.DataFrame(dict(vasc_A=[2.0, 2.0], t_vasc0=[0.0, 0.0], years=[0.0, 2.0]))          # written before 0.4.0
    assert lesion_tau(old) == 2.0 and meta(old, "a_scale", 1.0) == 1.0
    assert np.allclose(effective_lesion(old), [0.0, 2.0 * (1 - np.exp(-1.0))])
    new = old.assign(tau=4.0, a_scale=3.0)
    assert lesion_tau(new) == 4.0 and np.allclose(effective_lesion(new)[1], 6.0 * (1 - np.exp(-0.5)))
    assert np.all(effective_lesion(old.assign(vasc_A=0.0)) == 0)
