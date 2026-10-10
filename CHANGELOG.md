# Changelog

## 0.2.2 (2026-10)

- Fibrosis-only cohorts of the manuscript (`hom`, 60 cases; `het`, 38 cases) added to `results/ph_cohorts/` with `src/fibrosis_only_cohorts.py`, which regenerates them exactly.
- `src/ph_stats.py`: one statistics routine for Table 1, the supplementary tables and `ph_cohort.py analyze`, so identical quantities have identical intervals.
- Figure 4c reports the mean AUC over 1,000 realisations of measurement noise; Figure 2a y-axis renamed "Share (%)"; Figure 1a labels the sub-tree limit as a 10 µm radius.
- `pulmonary_shared_law.py` and `pulmonary_veins_shared_law.py` use the effective pulsatilities of the manuscript (φ₁ = 0.8 arteries, 0.4 veins).
- `run_all.sh` reproduces every figure and table from the shipped cohorts; `data/virtual_ph_cohort.tsv` keeps three decimals.
- THEORY.md numbers aligned with the manuscript; font fallback for systems without Liberation Sans.
- `src/make_supplement_figures.py`: figure S1 (FVC and DLCO along the trajectories and the window of the cross-sectional visit) and tables S6-S7 (inter-individual variation of the healthy baseline; lung function by stage).
- Figure 1: less empty space below the network schematic.
- `manuscript/`: current draft, online supplement, figures and tables; `docs/AUDIT.md`: reproducibility audit.

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
