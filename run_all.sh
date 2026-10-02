#!/usr/bin/env bash
# Reproduces the results of LuPNet's pulmonary-vascular analyses and rebuilds the calculator.
# The two virtual cohorts take about 2 hours on one core; runs are resumable.
set -e
python3 src/pulmonary_shared_law.py                         # shared target law vs human arterial morphometry (Huang 1996)
python3 src/pulmonary_veins_shared_law.py                   # the same law vs the venous tree
python3 src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85   # calibrated cohort
python3 src/ph_cohort.py run --tag adv2 --n 60  --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45   # advanced disease
python3 src/ph_cohort.py analyze --tag big2
python3 src/ph_cohort.py analyze --tag adv2
python3 src/fit_calculator.py                               # data/virtual_ph_cohort.tsv and results/ph_calculator.json
python3 src/build_app.py                                    # index.html
