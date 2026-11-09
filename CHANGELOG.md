# Changelog

## 0.3.11 (2026-11)

- Cohorts of 200 cases with the UIP and the fibrotic NSIP pattern (`big2_uip`, `big2_nsip`) and the three manuscript cohorts at the barometric pressure of Mexico City (`big2_cdmx`, `adv2_cdmx`, `altF_cdmx`), same seeds as the manuscript. PH (mPAP >= 25, FVC 50-90%) is 17% in both patterns (19% calibrated); UIP meets the seven targets, NSIP six (DLCO ratio 0.85, a target measured in IPF).
- `noninvasive_tool.py --reference pattern` (default): UIP and NSIP patients are scored against the cohort of their pattern; the partial dependence surface is fitted per pattern.
- `cohort_io.cases` reads only `tag_<digits>.tsv`; `ph_cohort.py analyze` no longer fails on cohorts without cases in the FVC 50-90% band; `scripts/cohorts_local.sh` gains `WHAT=advanced` (adv2_uip, adv2_nsip).

## 0.3.10 (2026-11)

- `src/noninvasive_tool.py`: individual estimate of PH without catheterisation. FVC and DLCO required; echocardiographic sPAP, pulmonary artery diameter on CT, BNP, RV fractional area change and right atrial area optional in any combination. Prior from the simulated visits of the manuscript cohorts (5,036 visits), measurement models and prevalence fitted on catheterised patients, repeated 5-fold cross-validation, comparison with logistic regression on the same inputs (tables S21-S23). The clinical file is passed with `--clinical` and is not part of the repository.
- Partial dependence: in the simulations mPAP is a surface of two variables, lung lost and effective lesion (R2 0.998). The tool reports for each patient the pressure of the aged lung, what the lost lung adds and what the lesion adds.

## 0.3.9 (2026-11)

- Histological patterns. `ph_cohort.py run --pattern {calibrated,uip,nsip}` (default `calibrated`, the manuscript cohorts, unchanged). UIP: subpleural, basal-weighted seed (`sub_gain = 3`), spread by contiguity (`kn = 1.5`, `kn_s = 3`), heterogeneous vulnerability, honeycombing. Fibrotic NSIP: no seed and no contiguity, near-uniform vulnerability (`at2_sigma = 0.3`), almost no honeycombing (`lam_hc = 0.05`), softer scar (`fib_C = 0.30`). The calibrated set is patchy without subpleural predominance (lost units 0.46 subpleural vs 0.44 core).
- New parameter `sub_gain` and `sim.subpleural_factor`: multiplier on surfactant decline and on the direct fibrosis hazard in the subpleural ring, weighted 1.5 at the base to 0.5 at the apex; 1 everywhere at the default `sub_gain = 0`.
- `ph_cohort.py run --PB` (barometric pressure; 585 mmHg for Mexico City) and `scripts/cohorts_local.sh`, which generates the pattern and altitude cohorts on the cores of one machine without a scheduler.
- In 12 paired cases (same seeds and lesions), mPAP differs by +0.3 mmHg between NSIP and UIP (Wilcoxon P = 1.00): the pattern sets where tissue is lost and the DLCO-FVC relation, not the pressure.

## 0.3.8 (2026-11)

Audit of 0.3.7 before the comparison with clinical data. No published number changes; the regression tests pass unchanged.

- Fix: Figure 2f computed the lesion clock from `t x ypu`, without the 6.71 years per model unit, so the lesion ran 6.7 times too slowly along the trajectories. The stored `years` column is now used.
- Figure 2g draws the mPAP of the aged lung without fibrosis (16.5 mmHg, `age_M = 0.3`) as a reference line, and the script prints the contribution of fibrosis above it (1.8, 1.0 and 0.4 mmHg by band) alongside the published values above the healthy lung (3.5, 2.7, 2.2).
- `src/audit_extra.py` and tables S17-S20: decomposition of the healthy pressure drop (47% set by the law), the aged baseline, the saturation of the lesion at the visit (median 1.00), and the low-DLCO phenotype under the baseline spread of table S6 (delta BIC median 7, bimodal in 40% of replicates). `tests/test_audit_extra.py` guards them.
- `ph_cohort.py run` takes `--t_max` and reports cases that stop before the drawn stage (7 of 60 in adv2).
- `Lung._baseline` and `sim.run` restore the fields of `Params` they alter in `finally` blocks.
- THEORY.md, docs/AUDIT.md and README corrected: the cohorts use the product form of the unit diffusing capacity (`rf = False`), not Roughton-Forster; DLNO is not computed; the morphometric fits use a free tau0 while the lung fixes it on the main pulmonary artery; healthy mPAP is an output of the law for its arterial half; ventilation is fixed and a quarter of the calibrated cases are hypercapnic; `beta = 0` (no strain feedback) in the calibrated set.

## 0.3.7 (2026-10)

- `ph_cohort.py run` takes `--i0` and `--i1`, so one cohort can be split across cores or cluster tasks. The draws depend only on `--n` and `--seed`, so a shard computes exactly the cases the whole run would have computed; a shard was checked to be identical, row by row, to the serial run of the same case.
- Each case is written atomically (`.tmp` then rename), so a task killed mid-write cannot leave a truncated file that the resume logic would accept as finished.
- `scripts/sensitivity_parallel.sh` runs the four sensitivity cohorts on every available core and then builds Table S16 and Figure S2; `scripts/sensitivity.slurm` does the same as a SLURM array. Both pin numpy to one thread per worker and are resumable.
- `sensitivity.py` reports a cohort with too few cases in the FVC 50-90% band instead of failing inside the statistics.
- The README citation gave v0.3.3 while the repository was at v0.3.4; it now gives the current version and the concept DOI.

## 0.3.6 (2026-10)

- `src/sensitivity.py`: the seven benchmarks recomputed across cohorts generated at other values of the lesion time constant (T_les) and of the proportion of cases that carry a lesion, with the generating commands in `run_all.sh`. Writes `TableS16_sensitivity.tsv` and `FigureS2`. Keeping the seeds fixed means a difference between rows comes from the parameter and not from a different draw of cases. The script reads cohorts and does not simulate; the four runs take about 5 h of one core.

## 0.3.5 (2026-10)

- Figure 2g: inside each FVC band the reference is now the cases of the same cohort that carry no lesion, instead of a fit to the fibrosis-only cohort. Those cases are the same fibrosis without vasculopathy, so the difference is what the lesion adds. The published values follow: fibrosis 3.5 to 2.1 mmHg, vasculopathy 6.1, 5.3 and 6.7 mmHg.
- `tests/test_manuscript_numbers.py`: the contributions of Figure 2g are now part of the regression test.
- The shear law is cited to its preprint (Pérez-Calixto et al., bioRxiv 2026, doi:10.64898/2026.09.29.755319; code at github.com/Danpc11/InFlow) in README, THEORY and `.zenodo.json`, instead of an unpublished manuscript.

## 0.3.4 (2026-10)

- `src/make_figures_erj.py`: Figure 2 gains panels f, g and h, the two clocks of fibrosis and vasculopathy. f) how far the lesion has run at each extent of fibrosis, along the trajectories of 14 simulated cases with an independent lesion and 14 with a lesion tied to fibrosis; g) what each process contributes to mPAP by FVC band, with the contribution of every case shown; h) mPAP against FVC in the alternative mechanism, with the pressure that fibrosis alone produces. The panel prints its own numbers so that the text can be checked against them.
- `mPAP_fibrosis_only()`: the pressure of the fibrosis-only cohort as a function of FVC, now the single baseline behind panels g and h.
- `tests/test_manuscript_numbers.py`: regression test of the published quantities (the seven calibration targets, the fibrosis-only and variable-response cohorts with their subset sizes, advanced disease out of sample, and the alternative mechanism). A change to the model that moves a number in the paper now fails the tests instead of passing unnoticed.
- `PANEL_UPPER` switches every panel letter to upper case for journals that require it.
- Tests run on Python 3.10 and 3.12; `requirements-dev.txt` and `requirements-lock.txt` added, the latter pinning the versions that produced the published numbers.
- `params.py`: comments on `vasc_tau` and `cap_frac` corrected; they described a time constant and a functional form the cohorts do not use.

## 0.3.3 (2026-10)

- Figure 5 and tables S13-S15 include resting SpO2; legends of Figure 5 placed below the panels; Figure 1a states the 70% prevalence of vasculopathy.
- `src/manuscript_numbers.py`: every number quoted in the text of the manuscript and supplementary table S10 (PH by FVC range), recomputed from the stored cohorts.
- `src/validate_calculator.py`: validation of the calculator on held-out simulated cases (the figures quoted in the calculator).
- README calculator section, THEORY and the audit aligned with the manuscript.

## 0.3.2 (2026-10)

- Calculator: resting SpO2 added as an optional input (with echocardiographic systolic pressure and the A-a gradient); responsive layout for phones, tablets and desktops; the pressure estimate is shown as a scale with its 90% interval and the pressure expected from fibrosis alone, and the probabilities as bars; the FVC-DLCO map is drawn as SVG at the size of the screen.

## 0.3.1 (2026-10)

- `src/noninvasive.py`: emulator of the simulator (gradient-boosted trees, ~4 µs per evaluation), identification of relevant vasculopathy and PH with echocardiographic systolic pressure and the alveolar-arterial O2 gradient, and individual Bayesian estimation of mPAP with 90% intervals (Figure 5, tables S13-S15).
- Calculator: FVC and DLCO required; echocardiographic systolic pressure and A-a gradient optional. Estimates weight 5,036 simulated visits by their likelihood given the patient's data and report the estimated mPAP with its 90% interval and the probabilities of PH and of relevant vasculopathy.

## 0.3.0 (2026-10)

- Honeycomb units keep 15% of normal compliance (`hc_C = 0.15`, previously 0.60), so FVC keeps falling as honeycombing extends; this removes an artificial FVC plateau near 42% and the overestimate of PH prevalence at low FVC.
- Vasculopathy prevalence recalibrated to 70% (`--p_vasc 0.7`); the other four vasculopathy parameters are unchanged. All seven mild-to-moderate targets are inside their 95% intervals; PH with FVC < 50% is 48% (transplant candidates 46%).
- All five cohorts re-simulated with the final parameters; figures, tables and the calculator rebuilt. Figure 3 covers FVC down to 20%.
- Supplementary table S10 is the prevalence of PH by FVC range.

## 0.2.3

- Mathematical invariants of the model as tests (`tests/test_invariants.py`): flow conservation, pressure identity, Wood-unit conversion, uniqueness of the shared-law radius, rest point and normalisation.
- Audit sensitivity analyses in `src/make_supplement_figures.py`: wedge pressure (table S8), low-DLCO phenotype without vasculopathy labels (S9), honeycomb compliance and the FVC plateau (S10), pressure categories with measurement noise (S11).
- `docs/AUDIT.md`: mathematical and methodological audit.

## 0.2.2 (2026-10)

- Fibrosis-only cohorts of the manuscript (`hom`, 60 cases; `het`, 38 cases) added to `results/ph_cohorts/` with `src/fibrosis_only_cohorts.py`, which regenerates them exactly.
- `src/ph_stats.py`: one statistics routine for Table 1, the supplementary tables and `ph_cohort.py analyze`, so identical quantities have identical intervals.
- Figure 4c reports the mean AUC over 1,000 realisations of measurement noise; Figure 2a y-axis renamed "Share (%)"; Figure 1a labels the sub-tree limit as a 10 µm radius.
- `pulmonary_shared_law.py` and `pulmonary_veins_shared_law.py` use the effective pulsatilities of the manuscript (φ₁ = 0.8 arteries, 0.4 veins).
- `run_all.sh` reproduces every figure and table from the shipped cohorts; `data/virtual_ph_cohort.tsv` keeps three decimals.
- THEORY.md numbers aligned with the manuscript; font fallback for systems without Liberation Sans.
- `src/make_supplement_figures.py`: figure S1 (FVC and DLCO along the trajectories and the window of the cross-sectional visit) and tables S6-S7 (inter-individual variation of the healthy baseline; lung function by stage).
- Figure 1: less empty space below the network schematic.
- `docs/AUDIT.md`: reproducibility audit (independent recomputation of every number in the manuscript).
- **License changed** from MIT to the PolyForm Noncommercial License 1.0.0: noncommercial use only; commercial use requires a separate license from the authors. Versions up to 0.2.1 remain under MIT.

## 0.2.1 (2026-10)

- Archival release for Zenodo; code, data and simulations are identical to 0.2.0.
- Manuscript status in the README updated to "in preparation".

## 0.2.0 (2026-10)

- Alternative vascular mechanism: `vasc_mode="fibrosis"` (lesion scaled by the extent of fibrosis), `--vmode` in `src/ph_cohort.py`.
- Figures and tables of the manuscript (`src/make_figures_erj.py`, `src/supplement_tables.py`), including the network schematic and the 2022 haemodynamic definitions.
- The simulations of the three cohorts used in the manuscript are included as `results/ph_cohorts/{big2,adv2,altF}.tsv.gz`; `src/cohort_io.py` reads cohorts stored per case or packed per cohort.
- Terminology: model output is described as "simulated cases" and "in silico cohorts"; "patients" is reserved for clinical data.

## 0.1.0 (2026-10)

- Lung network model (`src/lupnet`): parenchymal units with surfactant and collapse, gas exchange with Roughton–Forster
  membrane and capillary components, fibrosis progression with ageing, right-ventricular adaptation, exacerbations.
- Shared target law for the pulmonary vessels (`law="shared"`): one shear constant for every vessel, a target that depends on
  radius and segment length, exact Womersley sensing of the first harmonic. Validated against human arterial and venous
  morphometry (Huang et al. 1996).
- Structural vasculopathy independent of fibrosis (`vasc_A`, `vasc_tau`, `t_vasc0`, `cap_frac`) and virtual PH cohorts
  (`src/ph_cohort.py`).
- PH-IPF research calculator (`index.html`).
- Fix: the Womersley factor of the shared law was half its exact value; the effective pulsatilities (`phi1_a = 0.8`,
  `phi1_v = 0.4`) reproduce all published numbers exactly.
