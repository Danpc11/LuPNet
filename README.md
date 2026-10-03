# LuPNet
## Lung Perfusion Network model

### A physics-based model of the fibrotic lung and its pulmonary circulation

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/LuPNet/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/LuPNet/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/LuPNet/)
![Version](https://img.shields.io/badge/version-0.2.0-1f6feb)

---

## What this is

LuPNet treats the lung as two coupled networks, the alveolar units that exchange gas and the vascular tree that perfuses them, and follows idiopathic pulmonary fibrosis (IPF) as it removes units, stiffens tissue and remodels vessels. It asks one question about the pulmonary circulation: **why is pulmonary hypertension in IPF so loosely related to the extent of fibrosis?**

Two results answer it.

**1. One remodelling law sets the healthy pulmonary tree.** Every vessel remodels until the shear it senses, the mean plus the first harmonic of the pulse, equals a target that depends on its radius and length:

$$\tau^* = \tau_0\, r^{\,b-1}\,\ell^{\,(b-1)/2}, \qquad S=\bar\tau\,[1+w_1\,\varphi_1\,G(\alpha)]$$

With one constant fixed on the main pulmonary artery, the law predicts the diameters of the human pulmonary arterial tree (Huang et al. 1996, 16 orders from 20 µm to 3 cm) within 14%, the venous tree within 14–18% without refitting, and a normal mean pulmonary artery pressure of 14.8 mmHg as an output rather than an input.

**2. Pulmonary hypertension needs a vasculopathy of its own.** In in silico IPF cohorts, tissue loss alone keeps mPAP below about 21 mmHg and ties it tightly to FVC (r = −0.88), unlike patients. Adding a structural small-artery vasculopathy in half of the simulated cases, with its own course and an associated capillary loss, reproduces all seven targets of mild-to-moderate IPF inside their 95% intervals (200 simulated cases):

| | Model | 95% CI | Patients |
|---|---|---|---|
| mPAP < 20 mmHg | 54% | 46–61 | 51% |
| mPAP 20–25 | 29% | 22–36 | 30% |
| mPAP ≥ 25 | 17% | 12–23 | 19% |
| mPAP ≥ 35 | 4% | 1–7 | 4% |
| r(mPAP, DLCO) | −0.35 | −0.46 to −0.23 | ≈ −0.30 |
| r(mPAP, FVC) | −0.13 | −0.27 to 0.01 | ≈ 0 |
| DLCO with / without PH | 0.78 | 0.68–0.88 | 0.73 |

These seven numbers were fitted with five vasculopathy parameters. Not fitted: the prevalence of PH in advanced disease (47% vs 46% in transplant candidates in the calibrated cohort), the share of severe PH there, the absence of an FVC difference between patients with and without PH, and a distinct low-DLCO phenotype that concentrates almost all PH.

The derivation, the model and all numbers are in **[THEORY.md](THEORY.md)**.

---

## Calculator

A browser-based research calculator is available at **https://danpc11.github.io/LuPNet/**. It runs locally in the browser.

Enter FVC and DLCO (and optionally resting SpO₂) of a patient with IPF. The calculator returns

- the probability of mPAP ≥ 25 mmHg,
- the probability of the vascular phenotype (a DLCO lower than fibrosis alone explains),
- the mPAP expected from fibrosis alone and the DLCO expected for that FVC,
- the mPAP of the 25 most similar simulated cases,
- the Zisman 2007 formula, when SpO₂ is given, for comparison,

and places the patient on a map of the 260 simulated cases. It is a research tool built on simulated cases and does not replace right heart catheterization.

---

## Quick start

```bash
pip install -r requirements.txt
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

*A vascular contribution to pulmonary hypertension in idiopathic pulmonary fibrosis* (submitted). The three cohorts used in the paper are included as `results/ph_cohorts/big2.tsv.gz` (calibrated, 200 cases), `adv2.tsv.gz` (advanced disease, 60) and `altF.tsv.gz` (vasculopathy linked to fibrosis, 100), so figures and tables can be rebuilt without re-running the simulations:

```bash
python src/make_figures_erj.py        # Figures 1-4 and Table 1 -> results/figures/
python src/supplement_tables.py       # Supplementary tables S2, S4 and S5 -> results/figures/
```

To regenerate the cohorts from scratch (about 3 hours on one core, resumable):

```bash
python src/ph_cohort.py run --tag big2 --n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85
python src/ph_cohort.py run --tag adv2 --n 60  --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45
python src/ph_cohort.py run --tag altF --n 100 --seed 7000 --pseed 40000 --vmode fibrosis
```

Each simulated case is an independent run of the model, not a real patient.

## Repository

| Path | Content |
|---|---|
| `src/lupnet/` | The model: parameters, lung network and vessels, breathing, septal mechanics, anatomical tree, simulation loop |
| `src/run_experiment.py` | One disease simulation from a parameter file |
| `src/pulmonary_shared_law.py`, `src/pulmonary_veins_shared_law.py` | Tests of the shared law against human morphometry |
| `src/ph_cohort.py` | In silico PH cohorts (run and analyse; `--vmode fibrosis` for the alternative mechanism) |
| `src/make_figures_erj.py`, `src/supplement_tables.py` | Figures and tables of the manuscript |
| `results/ph_cohorts/*.tsv.gz` | Simulations of the three cohorts used in the manuscript (one compressed file per cohort) |
| `src/cohort_io.py` | Reads cohorts stored per case or packed per cohort |
| `src/fit_calculator.py`, `src/build_app.py`, `src/app_template.html` | Calculator inputs and page |
| `data/calibrated_ipf_params.json` | Calibrated parameters of aged IPF |
| `data/huang1996_pulmonary_*.tsv` | Human pulmonary arterial and venous morphometry |
| `data/virtual_ph_cohort.tsv` | Final visit of every simulated case |
| `results/ph_calculator.json` | Coefficients used by the calculator |

## Limitations

- In very advanced disease the simulated FVC does not fall below about 42%, so at a given low FVC the model carries more tissue loss than patients do; the overall PH prevalence in a dedicated advanced cohort is overestimated (67–71% vs 46%).
- The arterial morphometry comes from the casts of one lung, so the pulmonary exponent b has no confidence interval yet.
- The pulsatility of the pulmonary artery is an effective value, not a measured waveform.

## License

MIT for the code (`LICENSE`). Transcribed morphometric tables keep the terms of their publishers.
