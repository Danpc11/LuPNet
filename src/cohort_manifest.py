#!/usr/bin/env python3
"""Manifest of the simulated cohorts (data/cohorts_manifest.tsv): tag, cases, command, seeds, code version, md5.

    python3 src/cohort_manifest.py            # write the manifest
    python3 src/cohort_manifest.py --check    # verify that every shipped cohort matches its md5
Cohorts written before 0.4.0 do not store their generating parameters or the code version; the manifest records the
command, and `written_with` is "not recorded" for them. `reproduced_with` is the version with which every case was
simulated again and found identical (src/verify_cohorts.py, data/cohorts_verification.tsv).
"""
import argparse, hashlib, os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lupnet

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); F = os.path.join(ROOT, "results", "ph_cohorts")
OUT = os.path.join(ROOT, "data", "cohorts_manifest.tsv"); VER = os.path.join(ROOT, "data", "cohorts_verification.tsv")
from cohort_io import DESIGN
_args = lambda d: f"--n {d['n']} --seed {d['seed']} --pseed {d['pseed']} --lo {d['lo']} --hi {d['hi']}"
BIG, ADV, ALT = _args(DESIGN["big2"]), _args(DESIGN["adv2"]), f"--n {DESIGN['altF']['n']} --seed {DESIGN['altF']['seed']} --pseed {DESIGN['altF']['pseed']} --vmode fibrosis"
RUN = "ph_cohort.py run --tag {tag} "
SPEC = {  # tag: arguments
    "big2": BIG, "adv2": ADV, "altF": ALT,
    "hom": "fibrosis_only_cohorts.py --arm hom --n 60 (seeds 100 + case)", "het": "fibrosis_only_cohorts.py --arm het --n 38 (seeds 100 + case)",
    "big2_uip": BIG + " --pattern uip", "big2_nsip": BIG + " --pattern nsip",
    "adv2_uip": ADV + " --pattern uip", "adv2_nsip": ADV + " --pattern nsip",
    "big2_cdmx": BIG + " --PB 585", "adv2_cdmx": ADV + " --PB 585", "altF_cdmx": ALT + " --PB 585",
    "big2_scar": BIG + " --where scar", "big2_field": BIG + " --where field --a_scale 5.8", "adv2_field": ADV + " --where field --a_scale 5.8",
    "big2_mix50": BIG + " --where mixed --mix 0.25 --a_scale 2.0", "adv2_mix50": ADV + " --where mixed --mix 0.25 --a_scale 2.0",
    "tau1.0": BIG + " --tau 1.0", "tau3.0": BIG + " --tau 3.0", "pv0.50": BIG + " --p_vasc 0.50", "pv0.85": BIG + " --p_vasc 0.85",
}


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def build():
    rows = []
    V = pd.read_csv(VER, sep="\t").set_index("cohort") if os.path.exists(VER) else None
    for tag, args in SPEC.items():
        f = os.path.join(F, tag + ".tsv.gz"); d = pd.read_csv(f, sep="\t", low_memory=False)
        ver = str(d["lupnet_version"].iloc[0]) if "lupnet_version" in d else "not recorded"
        rep = lupnet.__version__ if V is not None and tag in V.index and V.loc[tag, "identical"] == V.loc[tag, "cases"] == d.case.nunique() else "not verified"
        rows.append(dict(cohort=tag, cases=int(d.case.nunique()), command=args if args.startswith("fibrosis_only") else RUN.format(tag=tag) + args,
                         lesion_priors="none" if tag in ("hom", "het") else "p_vasc 0.7, A lognormal(1.5, 1.1), t0 U(0, 10) y, tau 2 y, cap 0.05 unless the command says otherwise",
                         tree_units=256, written_with=ver, reproduced_with=rep, md5=md5(f)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
    if a.check:
        M = pd.read_csv(OUT, sep="\t"); bad = [r.cohort for r in M.itertuples() if md5(os.path.join(F, r.cohort + ".tsv.gz")) != r.md5]
        extra = sorted(set(f[:-7] for f in os.listdir(F) if f.endswith(".tsv.gz")) - set(M.cohort))
        if bad or extra:
            raise SystemExit(f"cohorts that differ from the manifest: {bad}; cohorts not in the manifest: {extra}")
        print(f"{len(M)} cohorts match data/cohorts_manifest.tsv")
    else:
        build().to_csv(OUT, sep="\t", index=False); print(OUT, "written")
