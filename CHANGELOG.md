# Changelog

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
