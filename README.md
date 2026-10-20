# LuPNet
## Lung Perfusion Network model

### A physics-based model of the fibrotic lung and its pulmonary circulation

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
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

With one constant fixed on the main pulmonary artery, the law predicts the diameters of the human pulmonary arterial tree (Huang et al. 1996, 16 orders from 20 µm to 3 cm) within 14%, the venous tree within 14–18% without refitting, and a normal mean pulmonary artery pressure of 14.8 mmHg as an output rather than an input.

**2. Pulmonary hypertension needs a vasculopathy of its own.** In in silico IPF cohorts, tissue loss alone keeps mPAP below about 19 mmHg in mild-to-moderate disease and ties it tightly to FVC (r = −0.92), unlike patients. Adding a structural small-artery vasculopathy in 70% of the simulated cases, mostly mild, with its own course and an associated capillary loss, reproduces all seven targets of mild-to-moderate IPF inside their 95% intervals (200 simulated cases):

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

## Quick start

```bash
pip install -r requirements.txt -r requirements-dev.txt
python -m pytest tests -q                       # fast checks
python src/pulmonary_shared_law.py              # the shared law against human arterial morphometry
python src/ph_cohort.py analyze --table data/virtual_ph_cohort.tsv --tag calibrated   # the shipped calibrated cohort
./run_all.sh                                    # everything, about 2 hours on one core (resumable)
```

```python
from lupnet import Params
from lupnet.lung import Lung
L = Lung(Params(law="shared", b=0.925))          # healthy lung, 1,024 units
print(L.base_state["vs"]["PPA"])                 # mean pulmonary artery pressure, mmHg
```

## Reproducing the manuscript

*A vascular contribution to pulmonary hypertension in idiopathic pulmonary fibrosis* (in preparation). The three cohorts used in the paper are included as `results/ph_cohorts/big2.tsv.gz` (calibrated, 200 cases), `adv2.tsv.gz` (advanced disease, 60) and `altF.tsv.gz` (vasculopathy linked to fibrosis, 100), together with the two cohorts without vasculopathy (`hom.tsv.gz`, 60 cases; `het.tsv.gz`, 38 cases with variable vascular responses), so figures and tables can be rebuilt without re-running the simulations. All cohorts use a 256-unit tree:

```bash
python src/make_figures_erj.py        # Figures 1-4 (Figure 2a-e the calibrated cohort, 2f-h the two clocks) and Table 1 -> results/figures/
python src/supplement_tables.py       # Supplementary tables S2, S4, S5 and the fibrosis-only cohorts -> results/figures/
python src/noninvasive.py             # Figure 5 and tables S13-S15 (echocardiography, blood gases, emulator)
python src/manuscript_numbers.py      # numbers quoted in the text and table S10
python src/validate_calculator.py     # calculator validation
```

To regenerate the cohorts from scratch (about 3 hours on one core, resumable):

```bash
python src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85
python src/ph_cohort.py run --tag adv2 --n 60  --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45
python src/ph_cohort.py run --tag altF --n 100 --seed 7000 --pseed 40000 --vmode fibrosis
python src/fibrosis_only_cohorts.py --arm hom --n 60
python src/fibrosis_only_cohorts.py --arm het --n 38
```

Each simulated case is an independent run of the model, not a real patient.

## Repository

| Path | Content |
|---|---|
| `src/lupnet/` | The model: parameters, lung network and vessels, breathing, septal mechanics, anatomical tree, simulation loop |
| `src/run_experiment.py` | One disease simulation from a parameter file |
| `src/pulmonary_shared_law.py`, `src/pulmonary_veins_shared_law.py` | Tests of the shared law against human morphometry |
| `src/ph_cohort.py` | In silico PH cohorts (run and analyse; `--vmode fibrosis` for the alternative mechanism) |
| `src/make_figures_erj.py`, `src/supplement_tables.py`, `src/make_supplement_figures.py` | Figures and tables of the manuscript and its supplement |
| `docs/AUDIT.md` | Reproducibility audit |
| `results/ph_cohorts/*.tsv.gz` | Simulations of the three cohorts used in the manuscript (one compressed file per cohort) |
| `src/cohort_io.py` | Reads cohorts stored per case or packed per cohort |
| `src/manuscript_numbers.py` | Every number quoted in the text of the manuscript, and supplementary table S10 |
| `src/validate_calculator.py` | Validation of the calculator on held-out simulated cases |
| `src/noninvasive.py` | Emulator, identification of vasculopathy and PH with echocardiographic pressure and arterial blood gases, individual estimation (Figure 5, tables S13–S15) |
| `tests/` | Mathematical invariants of the model and a regression test of every number published in the manuscript |
| `src/ph_stats.py` | The single statistics routine behind Table 1, the supplementary tables and `ph_cohort.py analyze` |
| `src/fibrosis_only_cohorts.py` | Cohorts without vasculopathy (fibrosis alone and variable vascular responses) |
| `src/fit_calculator.py`, `src/build_app.py`, `src/app_template.html` | Calculator inputs and page |
| `data/calibrated_ipf_params.json` | Calibrated parameters of aged IPF |
| `data/huang1996_pulmonary_*.tsv` | Human pulmonary arterial and venous morphometry |
| `data/virtual_ph_cohort.tsv` | Final visit of every simulated case |
| `results/ph_calculator.json` | Coefficients used by the calculator |

## Limitations

- The five vasculopathy parameters were chosen against the mild-to-moderate targets; the advanced-disease results and the phenotype analysis are out of sample.
- The arterial morphometry comes from the casts of one lung, so the pulmonary exponent b has no confidence interval yet.
- The pulsatility of the pulmonary artery is an effective value, not a measured waveform.

## How to cite

Toscano-Marquez F, Cisneros J, Maldonado M, Cervera A, Tovar H, Vázquez-Victorio G, Pardo A, Selman M, Pérez-Calixto D. *LuPNet: Lung Perfusion Network model of the fibrotic lung and its pulmonary circulation* (v0.3.3). Zenodo, 2026. [DOI of v0.3.3 to be added on release; earlier versions: https://doi.org/10.5281/zenodo.23125038]

## License

LuPNet is released under the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0) (`LICENSE`): it may be used, modified and shared for noncommercial purposes, including academic research, teaching and personal use. **Commercial use requires a separate license from the authors**; contact Daniel Pérez-Calixto (dperez@inmegen.gob.mx). Transcribed morphometric tables keep the terms of their publishers. Versions up to 0.2.1 were released under the MIT License.
