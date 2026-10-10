#!/usr/bin/env bash
# Reproduces every number, figure and table of the manuscript and rebuilds the calculator.
# Figures and tables use the cohorts shipped in results/ph_cohorts/*.tsv.gz; the commented lines regenerate them
# from scratch (about 4 hours on one core; runs are resumable).
set -e
python3 src/pulmonary_shared_law.py           # shared law vs human arterial morphometry (Huang 1996)
python3 src/pulmonary_veins_shared_law.py     # the same law vs the venous tree
# python3 src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85   # calibrated cohort
# python3 src/ph_cohort.py run --tag adv2 --n 60  --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45   # advanced disease
# python3 src/ph_cohort.py run --tag altF --n 100 --seed 7000 --pseed 40000 --vmode fibrosis       # alternative mechanism
# python3 src/fibrosis_only_cohorts.py --arm hom --n 60                                            # fibrosis alone
# python3 src/fibrosis_only_cohorts.py --arm het --n 38                                            # variable vascular responses
python3 src/make_figures_erj.py               # Figures 1-4 and Table 1
python3 src/supplement_tables.py              # Tables S2, S4, S5 and the fibrosis-only cohorts
python3 src/make_supplement_figures.py        # Figure S1 and tables S6-S11 (lung-function range, baseline variation, audit sensitivity analyses)
python3 src/fit_calculator.py                 # data/virtual_ph_cohort.tsv and results/ph_calculator.json
python3 src/build_app.py                      # index.html
