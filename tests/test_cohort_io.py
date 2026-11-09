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
