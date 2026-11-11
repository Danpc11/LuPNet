#!/usr/bin/env python3
"""Reference of the tables that run_all.sh writes to results/figures/ (data/table_checksums.tsv).

    python3 src/table_checksums.py            # record the tables now in results/figures/
    python3 src/table_checksums.py --check    # compare them with the record
The record is that of the locked environment (requirements-lock.txt); continuous integration runs run_all.sh there
and checks it. Tables are text with rounded numbers, so a different checksum means a different number. Tables
S13-S15 come from a gradient-boosted emulator whose sums depend on the processor and the number of threads; they
are compared number by number with a tolerance against the copies in data/reference_tables/. Tables S22-S24 and the
patients row of S26 need a clinical file and are not covered.
"""
import argparse, glob, hashlib, os, shutil
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); FIG = os.path.join(ROOT, "results", "figures")
OUT = os.path.join(ROOT, "data", "table_checksums.tsv"); REF = os.path.join(ROOT, "data", "reference_tables")
SKIP = ("TableS22", "TableS23", "TableS24")
NUMERIC = ("TableS13", "TableS14", "TableS15")          # emulator tables: compared with a tolerance
ABS_TOL, REL_TOL = 0.02, 0.03                           # 0.02 in an AUC; 3% in an error in mmHg


def md5(f):
    return hashlib.md5(open(f, "rb").read()).hexdigest()


def tables():
    files = sorted(glob.glob(os.path.join(FIG, "Table*.tsv")) + glob.glob(os.path.join(FIG, "manuscript_numbers.tsv")))
    return [f for f in files if not os.path.basename(f).startswith(SKIP)]


def close(now, ref):
    """Same shape, same text, numbers within the tolerance."""
    A, B = pd.read_csv(now, sep="\t"), pd.read_csv(ref, sep="\t")
    if A.shape != B.shape or list(A.columns) != list(B.columns):
        return False
    for c in A.columns:
        if pd.api.types.is_numeric_dtype(B[c]):
            x, y = A[c].to_numpy(float), B[c].to_numpy(float)
            if not np.all((np.abs(x - y) <= np.maximum(ABS_TOL, REL_TOL * np.abs(y))) | (np.isnan(x) & np.isnan(y))):
                return False
        elif not A[c].astype(str).equals(B[c].astype(str)):
            return False
    return True


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
    if not a.check:
        os.makedirs(REF, exist_ok=True); rows = []
        for f in tables():
            n = os.path.basename(f); num = n.startswith(NUMERIC)
            if num:
                shutil.copyfile(f, os.path.join(REF, n))
            rows.append(dict(table=n, how="numeric" if num else "exact", md5=md5(f)))
        pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False); print(f"{len(rows)} tables recorded in data/table_checksums.tsv")
    else:
        R = pd.read_csv(OUT, sep="\t"); bad = []
        for r in R.itertuples():
            f = os.path.join(FIG, r.table)
            if not os.path.exists(f):
                bad.append(f"{r.table}: missing")
            elif r.how == "exact" and md5(f) != r.md5:
                bad.append(f"{r.table}: checksum differs")
            elif r.how == "numeric" and not close(f, os.path.join(REF, r.table)):
                bad.append(f"{r.table}: numbers differ beyond the tolerance")
        if bad:
            raise SystemExit("tables that differ from the record:\n  " + "\n  ".join(bad))
        n_exact = int((R["how"] == "exact").sum())
        print(f"{len(R)} tables match the record ({n_exact} exactly, {len(R) - n_exact} within tolerance)")
