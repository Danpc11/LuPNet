# Changelog

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
