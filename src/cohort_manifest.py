#!/usr/bin/env python3
"""Manifest of the simulated cohorts (data/cohorts_manifest.tsv): tag, cases, command, seeds, code version, md5.

    python3 src/cohort_manifest.py            # write the manifest
    python3 src/cohort_manifest.py --check    # verify that every shipped cohort matches its md5
Cohorts written before 0.4.0 do not store their generating parameters; the manifest records them. Sampled cases of
the manuscript cohorts were re-simulated with 0.4.0 and are identical to the stored ones.
"""
import argparse, hashlib, os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); F = os.path.join(ROOT, "results", "ph_cohorts")
OUT = os.path.join(ROOT, "data", "cohorts_manifest.tsv")
BIG, ADV, ALT = "--n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85", "--n 60 --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45", "--n 100 --seed 7000 --pseed 40000 --vmode fibrosis"
RUN = "ph_cohort.py run --tag {tag} "
SPEC = {  # tag: (arguments, written with)
    "big2": (BIG, "0.3.4"), "adv2": (ADV, "0.3.4"), "altF": (ALT, "0.3.4"),
    "hom": ("fibrosis_only_cohorts.py --arm hom --n 60 (seeds 100 + case)", "0.3.4"), "het": ("fibrosis_only_cohorts.py --arm het --n 38 (seeds 100 + case)", "0.3.4"),
    "big2_uip": (BIG + " --pattern uip", "0.3.10"), "big2_nsip": (BIG + " --pattern nsip", "0.3.10"),
    "adv2_uip": (ADV + " --pattern uip", "0.3.11"), "adv2_nsip": (ADV + " --pattern nsip", "0.3.11"),
    "big2_cdmx": (BIG + " --PB 585", "0.3.10"), "adv2_cdmx": (ADV + " --PB 585", "0.3.10"), "altF_cdmx": (ALT + " --PB 585", "0.3.10"),
    "big2_scar": (BIG + " --where scar", "0.3.13"), "big2_field": (BIG + " --where field --a_scale 5.8", "0.3.13"), "adv2_field": (ADV + " --where field --a_scale 5.8", "0.3.13"),
    "big2_mix50": (BIG + " --where mixed --mix 0.25 --a_scale 2.0", "0.3.13"), "adv2_mix50": (ADV + " --where mixed --mix 0.25 --a_scale 2.0", "0.3.13"),
    "tau1.0": (BIG + " --tau 1.0", "0.4.0"), "tau3.0": (BIG + " --tau 3.0", "0.4.0"), "pv0.50": (BIG + " --p_vasc 0.50", "0.4.0"), "pv0.85": (BIG + " --p_vasc 0.85", "0.4.0"),
}


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def build():
    rows = []
    for tag, (args, ver) in SPEC.items():
        f = os.path.join(F, tag + ".tsv.gz"); d = pd.read_csv(f, sep="\t", usecols=["case"])
        rows.append(dict(cohort=tag, cases=int(d.case.nunique()), command=args if args.startswith("fibrosis_only") else RUN.format(tag=tag) + args,
                         lesion_priors="none" if tag in ("hom", "het") else "p_vasc 0.7, A lognormal(1.5, 1.1), t0 U(0, 10) y, tau 2 y, cap 0.05 unless the command says otherwise",
                         tree_units=256, written_with=ver, md5=md5(f)))
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
