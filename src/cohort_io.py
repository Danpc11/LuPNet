"""Read in silico cohorts stored either as one file per simulated case (tag_000.tsv, ... written by ph_cohort.py)
or as one compressed file per cohort (tag.tsv.gz, with a 'case' column), which is how the cohorts of the manuscript
are distributed."""
import glob, os, re, warnings
import numpy as np, pandas as pd


def cases(folder, tag):
    """Yield (case_index, DataFrame of the run) for every simulated case of a cohort.

    Per-case files (tag_000.tsv, ...) are used when they are the only form of the cohort or when they are at least
    as many as the cases of the packed file; otherwise the packed file (tag.tsv.gz) is used and a warning says so,
    because a few per-case files next to a packed cohort are the trace of an interrupted re-run, not the cohort."""
    # exactly tag_<digits>.tsv: "big2_*.tsv" would also pick up the cases of big2_uip, big2_nsip or big2_cdmx
    pat = re.compile(re.escape(tag) + r"_\d+\.tsv$")
    files = sorted(f for f in glob.glob(os.path.join(folder, f"{tag}_*.tsv")) if pat.fullmatch(os.path.basename(f)))
    packed = os.path.join(folder, f"{tag}.tsv.gz")
    D = pd.read_csv(packed, sep="\t") if os.path.exists(packed) else None
    if files and D is not None and len(files) < D["case"].nunique():
        warnings.warn(f"cohort {tag}: {len(files)} per-case files but {D['case'].nunique()} cases in {tag}.tsv.gz; "
                      f"using the packed cohort (remove the per-case files or finish the run)")
        files = []
    if files:
        for f in files:
            yield int(os.path.basename(f).rsplit("_", 1)[1].split(".")[0]), pd.read_csv(f, sep="\t")
        return
    if D is not None:
        for i, d in D.groupby("case", sort=True):
            yield int(i), d.drop(columns="case").reset_index(drop=True)


def meta(d, key, default):
    """Generating parameter of a run (tau, cap, p_vasc, a_scale, ...). Cohorts written before 0.4.0 do not store
    them; `default` is then the value of the manuscript cohorts."""
    return type(default)(d[key].iloc[0]) if key in d and pd.notna(d[key].iloc[0]) else default


def lesion_tau(d):
    """Time constant of the lesion of a run, years (2.0 in every cohort written before 0.4.0)."""
    return meta(d, "tau", 2.0)


def effective_lesion(d, years=None):
    """Severity of the lesion in place along a run: A x a_scale x (1 - exp(-(years + t0) / tau)); 0 without lesion."""
    A = applied_severity(d); t0 = float(d.t_vasc0.iloc[0]) if "t_vasc0" in d else 0.0
    y = d.years.to_numpy() if years is None else years
    return A * (1.0 - np.exp(-(y + t0) / lesion_tau(d))) if A > 0 else np.zeros(np.shape(y))


def viable(d):
    """Boolean mask of the visits of a run that a living patient could have: oxygen demand met (gas_ok) and resting
    output above the stop threshold. Runs that end in a state that fails this are reported as deaths, not as visits."""
    ok = d.gas_ok.astype(bool) if "gas_ok" in d else pd.Series(True, index=d.index)
    if "Q" in d:
        ok = ok & (d.Q > 2.5)
    return ok


MANUSCRIPT = dict(p_vasc=0.7, a_med=1.5, a_sig=1.1, tau=2.0, t0_max=10.0, cap=0.05)   # lesion priors of the manuscript cohorts


def usable(d):
    """Visits that are viable and whose equilibration converged."""
    ok = viable(d)
    return ok & d.converged.astype(bool) if "converged" in d else ok


def last_visit(d):
    """(baseline row, last viable visit) of a run: the visit every script reports."""
    v = viable(d)
    return d.iloc[0], (d[v].iloc[-1] if v.any() else d.iloc[-1])


def applied_severity(d):
    """Severity of the lesion as applied: drawn severity x a_scale."""
    return float(d.vasc_A.iloc[0]) * meta(d, "a_scale", 1.0)


def final_visits(folder, tag):
    """One row per simulated case: the last viable visit, as % of the case's own baseline, and whether the run
    ended in a non-viable state (died)."""
    rows = []
    for i, d in cases(folder, tag):
        (r0, x), v = last_visit(d), viable(d)
        sc = meta(d, "a_scale", 1.0)
        rows.append(dict(case=i, mPAP=x.mPAP, PVR=x.PVR, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO,
                         vasc_A=float(x.vasc_A), vasc_A_scaled=float(x.vasc_A) * sc, died=bool(not v.iloc[-1])))
    return pd.DataFrame(rows)


def pack(folder, tag):
    """Combine tag_*.tsv into tag.tsv.gz."""
    parts = [d.assign(case=i) for i, d in cases(folder, tag)]
    pd.concat(parts, ignore_index=True).to_csv(os.path.join(folder, f"{tag}.tsv.gz"), sep="\t", index=False)
