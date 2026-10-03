"""Read in silico cohorts stored either as one file per simulated case (tag_000.tsv, ... written by ph_cohort.py)
or as one compressed file per cohort (tag.tsv.gz, with a 'case' column), which is how the cohorts of the manuscript
are distributed."""
import glob, os
import pandas as pd


def cases(folder, tag):
    """Yield (case_index, DataFrame of the run) for every simulated case of a cohort."""
    files = sorted(glob.glob(os.path.join(folder, f"{tag}_*.tsv")))
    if files:
        for f in files:
            yield int(os.path.basename(f).rsplit("_", 1)[1].split(".")[0]), pd.read_csv(f, sep="\t")
        return
    packed = os.path.join(folder, f"{tag}.tsv.gz")
    if os.path.exists(packed):
        D = pd.read_csv(packed, sep="\t")
        for i, d in D.groupby("case", sort=True):
            yield int(i), d.drop(columns="case").reset_index(drop=True)


def pack(folder, tag):
    """Combine tag_*.tsv into tag.tsv.gz."""
    parts = [d.assign(case=i) for i, d in cases(folder, tag)]
    pd.concat(parts, ignore_index=True).to_csv(os.path.join(folder, f"{tag}.tsv.gz"), sep="\t", index=False)
