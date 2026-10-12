# Changelog

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
