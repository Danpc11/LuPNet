#!/usr/bin/env bash
# Reproduces every number, figure and table of the manuscript and rebuilds the calculator.
# Figures and tables use the cohorts shipped in results/ph_cohorts/*.tsv.gz; the commented lines regenerate them
# from scratch (about 5 hours on one core; runs are resumable). Without simulations the script takes about 6 minutes.
set -e
python3 src/pulmonary_shared_law.py           # shared law vs human arterial morphometry (Huang 1996)
python3 src/pulmonary_veins_shared_law.py     # the same law vs the venous tree
# python3 src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85   # calibrated cohort
# python3 src/ph_cohort.py run --tag adv2 --n 60  --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45   # advanced disease
# python3 src/ph_cohort.py run --tag altF --n 100 --seed 7000 --pseed 40000 --vmode fibrosis       # alternative mechanism
# Sensitivity of the benchmarks to the two parameters the conclusion rests on (same seeds, so only the
# parameter changes). About 25 s per simulated case on one core; the runs are resumable.
# for t in 1.0 3.0; do python3 src/ph_cohort.py run --tag tau$t --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85 --tau $t; done       # sensitivity
# for p in 0.50 0.85; do python3 src/ph_cohort.py run --tag pv$p --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85 --p_vasc $p; done
# Pattern, altitude, lesion-location and sensitivity cohorts on several cores: scripts/cohorts_local.sh
# python3 src/fibrosis_only_cohorts.py --arm hom --n 60                                            # fibrosis alone
# python3 src/fibrosis_only_cohorts.py --arm het --n 38                                            # variable vascular responses
python3 src/cohort_manifest.py --check        # the shipped cohorts match data/cohorts_manifest.tsv
python3 src/make_figures_erj.py               # V21_Figure1-4 (manuscript v21) and Table 1
python3 src/supplement_tables.py              # Tables S2, S4, S5 and the fibrosis-only cohorts
python3 src/exercise_response.py analyze || echo "skipped: no pressure-flow cases (WHAT=flow bash scripts/cohorts_local.sh)"   # table S16, V21_Figure6
python3 src/tree_size.py analyze || echo "skipped: no 1,024-unit cases (WHAT=tree bash scripts/cohorts_local.sh)"            # table S6
python3 src/make_supplement_figures.py        # Figure S1 and tables S8-S12 (lung-function range, baseline variation, audit sensitivity analyses)
python3 src/noninvasive.py                    # V21_Figure5 and tables S13-S15 (echocardiography, blood gases, emulator)
python3 src/sensitivity.py --ref big2 --tags tau1.0 tau3.0 pv0.50 pv0.85 --labels 'T_les 1 y' 'T_les 3 y' '70% to 50%' '70% to 85%'   # table S17, figure S2
python3 src/make_figures_story.py             # Figure1-6 of the article (add --clinical FILE for Figure 5 and Figure 6D-F)
python3 src/manuscript_numbers.py            # numbers quoted in the text and table S11
python3 src/lesion_location.py                # tables S25-S26: where the lesion sits (add --clinical FILE for the patients row)
# python3 src/noninvasive_tool.py --clinical FILE   # tables S22-S24: individual estimate; needs the clinical file, not in the repository
python3 src/audit_extra.py                    # tables S18-S21: pressure decomposition, aged baseline, lesion saturation, phenotype robustness
python3 src/fit_calculator.py                 # data/virtual_ph_cohort.tsv and results/ph_calculator.json
python3 src/build_app.py                      # index.html
python3 src/validate_calculator.py            # calculator validation on held-out simulated cases
