# LuPNet
## Lung Perfusion Network model

### A physics-based model of the fibrotic lung and its pulmonary circulation

![Python](https://img.shields.io/badge/Python-3.9%2B%20%28locked%20versions%3A%203.11%2B%29-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/LuPNet/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/LuPNet/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/LuPNet/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23125038.svg)](https://doi.org/10.5281/zenodo.23125038)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-lightgrey)](LICENSE)

---

## What this is

LuPNet treats the lung as two coupled networks, the alveolar units that exchange gas and the vascular tree that perfuses them, and follows idiopathic pulmonary fibrosis (IPF) as it removes units, stiffens tissue and remodels vessels. It asks one question about the pulmonary circulation: **why is pulmonary hypertension in IPF so loosely related to the extent of fibrosis?**

Two results answer it.

**1. One remodelling law sets the healthy pulmonary tree.** Every vessel remodels until the shear it senses, the mean plus the first harmonic of the pulse, equals a target that depends on its radius and length:

$$\tau^* = \tau_0\, r^{\,b-1}\,\ell^{\,(b-1)/2}, \qquad S=\bar\tau\,[1+w_1\,\varphi_1\,G(\alpha)]$$

The law, its exponent and the weight of the pulsatile term are derived and tested across species in [Pérez-Calixto et al., bioRxiv 2026](https://doi.org/10.64898/2026.09.29.755319) ([InFlow](https://github.com/Danpc11/InFlow)).

Fitted to the human pulmonary arterial tree (Huang et al. 1996, 16 orders from 20 µm to 3 cm) with one free scale and exponent b = 0.90, the law predicts the diameters within 14%; applied to the venous tree without refitting, within 16% with a venous pulsatility φ₁ = 0.4 (27% and 32% with φ₁ = 0 and 0.8). In the lung model the scale is fixed on the main pulmonary artery and b = 0.925; the explicit arterial generations then match the measured ones within 6%, and, with the capillary and venular resistances of the reference calibration, the healthy mean pulmonary artery pressure is 14.8 mmHg with 1,024 units and default parameters (14.3 mmHg with the 256 units and calibrated parameters of the cohorts) (the law sets about half of the pressure drop; [THEORY.md](THEORY.md) §2.4).

**2. Fibrosis alone does not produce pulmonary hypertension; a vascular lesion does.** In in silico IPF cohorts, tissue loss alone keeps mPAP below about 19 mmHg in mild-to-moderate disease and ties it tightly to FVC (r = −0.92), unlike patients. Adding a lesion of the small arteries in 70% of the simulated cases, mostly mild, with an associated capillary loss, reproduces all seven targets of mild-to-moderate IPF inside their 95% intervals (200 simulated cases):

| | Model | 95% CI | Patients |
|---|---|---|---|
| mPAP < 20 mmHg | 54% | 45–63 | 51% |
| mPAP 20–25 | 27% | 19–35 | 30% |
| mPAP ≥ 25 | 19% | 12–26 | 19% |
| mPAP ≥ 35 | 5% | 2–10 | 4% |
| r(mPAP, DLCO) | −0.38 | −0.51 to −0.20 | ≈ −0.30 |
| r(mPAP, FVC) | −0.08 | −0.25 to 0.10 | ≈ 0 |
| DLCO with / without PH | 0.81 | 0.73–0.90 | 0.73 |

These seven numbers were fitted with five vasculopathy parameters. Not fitted: the prevalence of PH with FVC below 50% (48% vs 46% in transplant candidates) and of severe PH (8% vs about 9%), the small FVC difference between patients with and without PH, and a distinct low-DLCO phenotype that contains 71% of all PH.

**3. The lesion that matters sits in preserved lung, and only part of it can follow the scar.** Fibrosis thickens the vessels it surrounds, so the vasculopathy cannot be assumed independent of the fibrosis. LuPNet separates the possibilities by where the lesion sits (same fibrosis and same lesion draws in every cohort):

| Lesion | Benchmarks met | r(mPAP, FVC) | PH, FVC < 50% | mPAP > 40, FVC < 50% | PH, FVC ≥ 70% |
|---|---|---|---|---|---|
| Diffuse (manuscript) | 7 of 7 | −0.08 | 48% | 8% | 23% |
| Mixed: diffuse at half severity plus spreading (0.5 A diffuse + 1.5 A × field weight) | 7 of 7 | −0.17 | 54% | 16% | 14% |
| Spreading from the scar into adjacent lung | 4 of 7 | −0.34 | 68% | 31% | 0% |
| Inside the scar only | 1 of 7 | −0.19 | 0% | 0% | 0% |
| Patients | | ≈ 0 | 46% | ≈ 9% | |

A lesion confined to the scar produces no PH in 200 simulated cases (blood bypasses the scar). A lesion that spreads from the scar produces PH that follows lung function and triples severe PH in advanced disease. A minor share of the lesion can therefore be a local consequence of fibrosis; the rest sits in preserved lung and is not proportional to the extent of fibrosis.

**4. The same holds for the UIP and the fibrotic NSIP pattern, and at altitude.** With a subpleural, basal, honeycombing pattern (UIP) or a homogeneous one without honeycombing (NSIP), PH in mild-to-moderate disease is 17% in both and the pressure of a given lesion differs by 1 mmHg or less. At the barometric pressure of Mexico City (585 mmHg) pressures change by less than 1 mmHg; blood gases do change.

The derivation, the model and all numbers are in **[THEORY.md](THEORY.md)**.

---

## Calculator

A browser-based research calculator is available at **https://danpc11.github.io/LuPNet/**. It runs locally in the browser and works on phones, tablets and desktops.

Enter FVC and DLCO (required) and, if available, the echocardiographic systolic pulmonary artery pressure, the alveolar–arterial O₂ gradient and resting SpO₂. The calculator weights 5,036 simulated visits by how likely they are to produce the patient's data and returns

- the estimated mPAP with its 90% interval, against the mPAP expected from fibrosis alone,
- the probabilities of mPAP ≥ 25 mmHg and of relevant vasculopathy,
- the probability of the low-DLCO phenotype and the DLCO expected for that FVC,

and places the patient on a map of the 260 simulated cases. In held-out simulated cases, the area under the ROC curve for PH rises from 0.77 with FVC and DLCO to 0.94 with all optional inputs (`src/validate_calculator.py`). It is a research tool built on simulated cases and does not replace right heart catheterization.

---

## Non-invasive individual estimate

`src/noninvasive_tool.py` estimates PH for one patient without right heart catheterisation. FVC and DLCO are required; echocardiographic systolic pressure, pulmonary artery diameter on CT, BNP, right ventricular fractional area change and right atrial area are optional, in any combination. The simulated visits are the prior; catheterised patients are used only to fit how each measurement relates to pressure and to score the tool. In 43 catheterised patients with fibrosing interstitial lung disease (cross-validated), the area under the ROC curve for mPAP ≥ 25 mmHg is 0.70 with FVC and DLCO, 0.82 with the pulmonary artery diameter and 0.86 with every optional input (logistic regression on the same inputs: 0.64, 0.86 and 0.90, the last on the 34 patients with complete data), and the estimate of mPAP has a mean error of 8.6 mmHg. The tool also returns the pressure of the aged lung, what the lost lung adds and what the lesion adds, from the surface mPAP = F(lung lost, lesion) of the simulations (R² 0.998). The clinical file is not part of the repository:

```bash
python src/noninvasive_tool.py --clinical patients.xlsx      # tables S22-S24
python src/lesion_location.py --clinical patients.xlsx       # tables S25-S26
```

It is a research tool, evaluated in one cohort. It estimates pre-capillary PH and does not use SpO₂ or blood gases.

## Quick start

```bash
conda env create -f environment.yml && conda activate lupnet      # or: pip install -r requirements.txt -r requirements-dev.txt
python -m pytest tests -q                       # fast checks
python src/pulmonary_shared_law.py              # the shared law against human arterial morphometry
python src/ph_cohort.py analyze --table data/virtual_ph_cohort.tsv --tag calibrated   # the shipped calibrated cohort
./run_all.sh                                    # every figure and table from the shipped simulations, about 6 minutes
```

```python
from lupnet import Params
from lupnet.lung import Lung
L = Lung(Params(law="shared", b=0.925))          # healthy lung, 1,024 units
print(L.base_state["vs"]["PPA"])                 # mean pulmonary artery pressure, mmHg
```

## Reproducing the manuscript

*A vascular contribution to pulmonary hypertension in idiopathic pulmonary fibrosis* (in preparation). The cohorts used in the paper are included as `results/ph_cohorts/big2.tsv.gz` (calibrated, 200 cases), `adv2.tsv.gz` (advanced disease, 60) and `altF.tsv.gz` (vasculopathy linked to fibrosis, 100), together with the two cohorts without vasculopathy (`hom.tsv.gz`, 60 cases; `het.tsv.gz`, 38 cases with variable vascular responses), so figures and tables can be rebuilt without re-running the simulations. All cohorts use a 256-unit tree:

```bash
python src/make_figures_story.py      # Figures 1-6 of the article (add --clinical FILE for Figure 5 and Figure 6D-F; legends in docs/FIGURES.md)
python src/make_figures_erj.py        # figures 1-4 of manuscript v21 (V21_Figure1-4) and Table 1 -> results/figures/
python src/supplement_tables.py       # Supplementary tables S2, S4, S5 and the fibrosis-only cohorts -> results/figures/
python src/noninvasive.py             # V21_Figure5 and tables S13-S15 (echocardiography, blood gases, emulator)
python src/manuscript_numbers.py      # numbers quoted in the text and table S11
python src/validate_calculator.py     # calculator validation
python src/exercise_response.py analyze   # table S16 and V21_Figure6 (pressure-flow response of 70 cases; Figure 6A-B of the article)
python src/tree_size.py analyze           # table S6 (256 against 1,024 units, 12 cases)
python src/sensitivity.py --ref big2 --tags tau1.0 tau3.0 pv0.50 pv0.85 --labels 'T_les 1 y' 'T_les 3 y' '70% to 50%' '70% to 85%'   # table S17, figure S2
```

Supplementary tables carry the numbers of the online supplement. Up to version 0.3.14 the repository used its own numbering; the correspondence is

| Content | Supplement and v0.4.0 | Up to v0.3.14 |
|---|---|---|
| Tree of 256 and 1,024 units | S6 | not in the repository |
| Variation of the healthy baseline | S8 | S6 |
| Wedge pressure | S9 | S8 |
| Low-DLCO phenotype without labels | S10 | S9 |
| PH by FVC range | S11 | S10 |
| Pressure categories with noise | S12 | S11 |
| Emulator | S15 | S13a |
| Pressure–flow response | S16 | not in the repository |
| Sensitivity to the lesion time constant and prevalence | S17 | S16 |
| Audit checks (pressure decomposition, aged lung, lesion saturation, phenotype under baseline spread) | S18–S21 | S17–S20 |
| Non-invasive tool (cross-validation, measurement models, individual report) | S22–S24 | S21–S23 |
| Location of the lesion | S25–S26 | S24–S25 |

The model has more options than the manuscript uses. The manuscript cohorts use the symmetric tree, the mirrored venous tree, the shared law, gravity, surfactant mechanics with the time-domain breath, the right ventricle and vessel ageing (`data/calibrated_ipf_params.json`). The septal venous tree, the anatomical tree, epithelial senescence (`epi_bio`), acute exacerbations (`ae`), traction, the per-vessel law and the unbounded lesion (`k_vasc`) run and converge but were not used for any published number.

To regenerate the manuscript cohorts from scratch (about 5 hours on one core, resumable):

```bash
python src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85
python src/ph_cohort.py run --tag adv2 --n 60  --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45
python src/ph_cohort.py run --tag altF --n 100 --seed 7000 --pseed 40000 --vmode fibrosis
python src/fibrosis_only_cohorts.py --arm hom --n 60
python src/fibrosis_only_cohorts.py --arm het --n 38
```

Each simulated case is an independent run of the model, not a real patient.

## Cohorts in `results/ph_cohorts/`

All cohorts use a 256-unit tree and, within a family, the same seeds, so a difference between two cohorts comes from the one setting that differs. `big2*`: 200 cases, visit at 25–85% of healthy tissue; `adv2*`: 60 cases, 12–45%; `altF*`: 100 cases.

| Cohort | What differs from the manuscript cohort | Extra arguments of `ph_cohort.py run` |
|---|---|---|
| `big2`, `adv2` | Manuscript: calibrated pattern, diffuse lesion, sea level | none |
| `altF` | Lesion severity scaled by the extent of fibrosis | `--vmode fibrosis` |
| `hom`, `het` | No vasculopathy; variable vascular responses | `fibrosis_only_cohorts.py` |
| `big2_uip`, `adv2_uip` | UIP pattern | `--pattern uip` |
| `big2_nsip`, `adv2_nsip` | Fibrotic NSIP pattern | `--pattern nsip` |
| `big2_cdmx`, `adv2_cdmx`, `altF_cdmx` | Barometric pressure of Mexico City | `--PB 585` |
| `big2_scar` | Lesion in fibrotic and honeycomb units only | `--where scar` |
| `big2_field`, `adv2_field` | Lesion in the scar and the healthy units next to it | `--where field --a_scale 5.8` |
| `big2_mix50`, `adv2_mix50` | Mixed lesion: 0.5 A in every unit plus 1.5 A × field weight | `--where mixed --mix 0.25 --a_scale 2.0` |
| `tau1.0`, `tau3.0` | Time constant of the lesion 1 and 3 years | `--tau 1.0`, `--tau 3.0` |
| `pv0.50`, `pv0.85` | Lesion in 50% and 85% of the cases | `--p_vasc 0.50`, `--p_vasc 0.85` |

`scripts/cohorts_local.sh` regenerates the last six families on the cores of one machine without a scheduler (`WHAT=patterns`, `advanced`, `altitude`, `location` or `sensitivity`; about 25 to 60 s of one core per case, resumable).

## Repository

| Path | Content |
|---|---|
| `src/lupnet/` | The model: parameters, lung network and vessels, breathing, septal mechanics, anatomical tree, simulation loop |
| `src/run_experiment.py` | One disease simulation from a parameter file |
| `src/pulmonary_shared_law.py`, `src/pulmonary_veins_shared_law.py` | Tests of the shared law against human morphometry |
| `src/ph_cohort.py` | In silico PH cohorts (run and analyse; `--vmode fibrosis` for the alternative mechanism) |
| `src/make_figures_erj.py`, `src/supplement_tables.py`, `src/make_supplement_figures.py` | Figures and tables of the manuscript and its supplement |
| `docs/AUDIT.md` | Reproducibility audit |
| `results/ph_cohorts/*.tsv.gz` | Simulated cohorts, one compressed file each (see "Cohorts") |
| `src/cohort_io.py` | Reads cohorts stored per case or packed per cohort |
| `src/manuscript_numbers.py` | Every number quoted in the text of the manuscript, and supplementary table S11 |
| `src/validate_calculator.py` | Validation of the calculator on held-out simulated cases |
| `src/noninvasive.py` | Emulator, identification of vasculopathy and PH with echocardiographic pressure and arterial blood gases, individual estimation (V21_Figure5, tables S13–S15) |
| `scripts/sensitivity_parallel.sh`, `scripts/sensitivity.slurm` | Generate the four sensitivity cohorts across several cores, or as a SLURM array, then build Table S17 and Figure S2 |
| `src/sensitivity.py` | Benchmarks recomputed across cohorts generated at other values of the lesion time constant and of the proportion of cases carrying a lesion (Table S17, Figure S2) |
| `src/noninvasive_tool.py` | Individual estimate of PH without catheterisation: prior from the simulated visits, measurement models from catheterised patients, cross-validation (tables S22–S24) |
| `src/exercise_response.py` | Pressure–flow response at imposed cardiac output (table S16, V21_Figure6); `results/figures/pressure_flow_cases.tsv` holds every case |
| `src/tree_size.py` | The same cases with 256 and 1,024 units (table S6); `results/figures/tree_size_cases.tsv` |
| `data/README.md`, `CITATION.cff` | Provenance of the morphometric tables and the parameter file; how to cite |
| `src/make_figures_story.py` | Figures 1–6 of the article, on one grid and one palette; legends in `docs/FIGURES.md`. Figure 5 and Figure 6D–F need `--clinical` |
| `src/cohort_manifest.py`, `data/cohorts_manifest.tsv` | Command, seeds, code version and md5 of every cohort; `--check` verifies them |
| `src/table_checksums.py`, `data/table_checksums.tsv` | Checksums of the tables of `run_all.sh` in the locked environment; `--check` compares |
| `src/verify_cohorts.py` | Simulates the cases again from the manifest and compares them with the stored cohorts (`data/cohorts_verification.tsv`) |
| `src/synthetic_clinical.py`, `data/synthetic_patients.tsv` | Synthetic file with the columns of the clinical file, to run and test the scripts that need one; its rows are not patients |
| `src/lesion_location.py` | Benchmarks by location of the lesion and PH by FVC band (tables S25–S26) |
| `scripts/cohorts_local.sh` | Pattern, altitude and lesion-location cohorts on the cores of one machine |
| `environment.yml` | Conda environment with the locked versions |
| `tests/` | Mathematical invariants of the model and a regression test of every number published in the manuscript |
| `src/audit_extra.py` | Checks of the 0.3.7 audit: where the healthy pressure drop comes from, the aged lung without fibrosis, lesion saturation at the visit, the low-DLCO phenotype under baseline spread (tables S18–S21) |
| `src/ph_stats.py` | The single statistics routine behind Table 1, the supplementary tables and `ph_cohort.py analyze` |
| `src/fibrosis_only_cohorts.py` | Cohorts without vasculopathy (fibrosis alone and variable vascular responses) |
| `src/fit_calculator.py`, `src/build_app.py`, `src/app_template.html` | Calculator inputs and page |
| `data/calibrated_ipf_params.json` | Calibrated parameters of aged IPF |
| `data/huang1996_pulmonary_*.tsv` | Human pulmonary arterial and venous morphometry |
| `data/virtual_ph_cohort.tsv` | Final visit of every simulated case |
| `results/ph_calculator.json` | Coefficients used by the calculator |

## Data availability

The simulated cohorts, the code and every derived table are in this repository. The clinical data are not: tables S22–S24, Figure 5, Figure 6D–F and the patients row of table S26 need a clinical file passed with `--clinical` (columns in `noninvasive_tool.COLS`; `data/synthetic_patients.tsv` has the same layout with synthetic rows; what is computed from it goes to `results/synthetic/` and its figures carry a watermark). Tables S1 and S3 of the supplement are not produced by a script: S1 comes from the companion preprint and S3 lists the parameters of `data/calibrated_ipf_params.json`. `TableS_fibrosis_only.tsv` holds the numbers of the cohorts without lesion quoted in the text.

Two sets of figures are produced. `Figure1`–`Figure6` (`make_figures_story.py`) are the figures of the article. `V21_Figure1`–`V21_Figure6` are those of manuscript v21, kept because the regression test and the audit refer to them.

## Limitations

- The lesion is a resistance added to the arteriole, not a cellular mechanism; its parameters were chosen against the mild-to-moderate benchmarks. Advanced disease, the pressure–flow response and the patients are out of sample.
- Pre-capillary pulmonary hypertension only; alveolar ventilation is fixed, so blood gases are not used for the individual estimate.
- The arterial morphometry comes from the casts of one lung; the exponent b has no confidence interval. The pulsatility of the pulmonary artery is an effective value.
- The vessels of the calibrated lung are aged: without fibrosis it sits at 16.5 mmHg.
- The individual estimate was evaluated in one cohort, referred for suspected PH.

The full list is in [THEORY.md](THEORY.md), section 8.

## How to cite

Toscano-Marquez F, Cisneros J, Maldonado M, Cervera A, Tovar H, Vázquez-Victorio G, Pardo A, Selman M, Pérez-Calixto D. *LuPNet: Lung Perfusion Network model of the fibrotic lung and its pulmonary circulation* (v0.4.0). Zenodo, 2026. [concept DOI, resolves to the latest version: https://doi.org/10.5281/zenodo.23125038]

## License

LuPNet is released under the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0) (`LICENSE`): it may be used, modified and shared for noncommercial purposes, including academic research, teaching and personal use. **Commercial use requires a separate license from the authors**; contact Daniel Pérez-Calixto (dperez@inmegen.gob.mx). Transcribed morphometric tables keep the terms of their publishers. Versions up to 0.2.1 were released under the MIT License.
