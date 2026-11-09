"""Read in silico cohorts stored either as one file per simulated case (tag_000.tsv, ... written by ph_cohort.py)
or as one compressed file per cohort (tag.tsv.gz, with a 'case' column), which is how the cohorts of the manuscript
are distributed."""
import glob, os, re
import pandas as pd


def cases(folder, tag):
    """Yield (case_index, DataFrame of the run) for every simulated case of a cohort."""
    # exactly tag_<digits>.tsv: "big2_*.tsv" would also pick up the cases of big2_uip, big2_nsip or big2_cdmx
    pat = re.compile(re.escape(tag) + r"_\d+\.tsv$")
    files = sorted(f for f in glob.glob(os.path.join(folder, f"{tag}_*.tsv")) if pat.fullmatch(os.path.basename(f)))
    if files:
        for f in files:
            yield int(os.path.basename(f).rsplit("_", 1)[1].split(".")[0]), pd.read_csv(f, sep="\t")
        return
    packed = os.path.join(folder, f"{tag}.tsv.gz")
    if os.path.exists(packed):
        D = pd.read_csv(packed, sep="\t")
        for i, d in D.groupby("case", sort=True):
            yield int(i), d.drop(columns="case").reset_index(drop=True)


def viable(d):
    """Boolean mask of the visits of a run that a living patient could have: oxygen demand met (gas_ok) and resting
    output above the stop threshold. Runs that end in a state that fails this are reported as deaths, not as visits."""
    ok = d.gas_ok.astype(bool) if "gas_ok" in d else pd.Series(True, index=d.index)
    if "Q" in d:
        ok = ok & (d.Q > 2.5)
    return ok


def final_visits(folder, tag):
    """One row per simulated case: the last viable visit, as % of the case's own baseline, and whether the run
    ended in a non-viable state (died)."""
    rows = []
    for i, d in cases(folder, tag):
        r0 = d.iloc[0]; v = viable(d); x = d[v].iloc[-1] if v.any() else d.iloc[-1]
        sc = float(d["a_scale"].iloc[0]) if "a_scale" in d else 1.0
        rows.append(dict(case=i, mPAP=x.mPAP, PVR=x.PVR, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO,
                         vasc_A=float(x.vasc_A), vasc_A_scaled=float(x.vasc_A) * sc, died=bool(not v.iloc[-1])))
    return pd.DataFrame(rows)


def pack(folder, tag):
    """Combine tag_*.tsv into tag.tsv.gz."""
    parts = [d.assign(case=i) for i, d in cases(folder, tag)]
    pd.concat(parts, ignore_index=True).to_csv(os.path.join(folder, f"{tag}.tsv.gz"), sep="\t", index=False)
