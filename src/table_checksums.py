#!/usr/bin/env python3
"""Checksums of the tables that run_all.sh writes to results/figures/ (data/table_checksums.tsv).

    python3 src/table_checksums.py            # record the checksums of the tables now in results/figures/
    python3 src/table_checksums.py --check    # compare them with the recorded ones
The recorded checksums are those of the locked environment (requirements-lock.txt); continuous integration runs
run_all.sh there and checks them. Tables S22-S24 and the patients row of S26 need a clinical file and are not
covered. Tables are text with rounded numbers, so a different checksum means a different number.
"""
import argparse, glob, hashlib, os
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); FIG = os.path.join(ROOT, "results", "figures")
OUT = os.path.join(ROOT, "data", "table_checksums.tsv")
SKIP = ("TableS22", "TableS23", "TableS24")


def current():
    files = sorted(glob.glob(os.path.join(FIG, "Table*.tsv")) + glob.glob(os.path.join(FIG, "manuscript_numbers.tsv")))
    return pd.DataFrame([dict(table=os.path.basename(f), md5=hashlib.md5(open(f, "rb").read()).hexdigest())
                         for f in files if not os.path.basename(f).startswith(SKIP)])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
    C = current()
    if not a.check:
        C.to_csv(OUT, sep="\t", index=False); print(f"{len(C)} tables recorded in data/table_checksums.tsv")
    else:
        R = pd.read_csv(OUT, sep="\t").merge(C, on="table", how="left", suffixes=("_recorded", "_now"))
        bad = R[R.md5_recorded != R.md5_now]
        if len(bad):
            raise SystemExit("tables that differ from data/table_checksums.tsv or are missing:\n" + bad.to_string(index=False))
        print(f"{len(R)} tables match data/table_checksums.tsv")
