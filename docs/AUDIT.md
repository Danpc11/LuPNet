# Reproducibility audit

Manuscript: *A vascular contribution to pulmonary hypertension in idiopathic pulmonary fibrosis* (in preparation; not part of this repository).
Code: this repository, version 0.3.0.

## Method

1. The repository was cloned into a clean environment. The tests (4/4) and `run_all.sh` were run in full: arterial and venous morphometry, Figures 1–4, Table 1, supplementary tables, supplementary figure S1 and the analysis of the five cohorts.
2. An independent script, separate from the repository code, recomputed every number quoted in the manuscript from the stored simulations (`results/ph_cohorts/*.tsv.gz`).
3. Model physics and units were checked against the baseline values of individual simulated cases.

## Reproducibility

Every number in the manuscript is reproduced exactly from the stored simulations.

| Quantity | Manuscript | Recomputed |
|---|---|---|
| Calibrated cases / with FVC 50–90% | 200 / 113 | 200 / 113 |
| Realised prevalence of vasculopathy (calibrated / advanced / alternative) | 70% / – / 64% | 70% / 72% / 64% |
| mPAP < 20 mmHg | 54% (45–63) | 54.0% (45.1–62.8) |
| mPAP 20–25 mmHg | 27% (19–35) | 27.4% (19.5–35.4) |
| mPAP ≥ 25 mmHg | 19% (12–26) | 18.6% (11.5–25.7) |
| mPAP ≥ 35 mmHg | 5% (2–10) | 5.3% (1.8–9.7) |
| r(mPAP, DLCO) | −0.38 (−0.51 to −0.20) | −0.38 (−0.51 to −0.20) |
| r(mPAP, FVC) | −0.08 (−0.25 to 0.10) | −0.08 (−0.25 to 0.10) |
| DLCO ratio, PH / no PH | 0.81 (0.73–0.90) | 0.81 (0.73–0.90) |
| DLCO ratio at the 21 mmHg threshold | 0.84 (0.76–0.92) | 0.84 (0.76–0.92) |
| No vasculopathy: mPAP | 17.7 ± 0.6, max 19.2 | 17.7 ± 0.6, 19.2 |
| Vasculopathy: mPAP, % PH | 23.7 ± 5.8, 28% | 23.7 ± 5.8, 28% |
| PH with FVC < 50%, pooled (n = 147) | 48% (39–56) | 48% (39–56) |
| Severe PH (mPAP > 40 mmHg), pooled | 8% (4–13) | 8% (4–13) |
| FVC with / without PH, FVC < 50% | 35.7 / 39.1% | 35.7 / 39.1% |
| 2022 precapillary PH, mild-to-moderate | 46% (37–55) | 46% (37–55) |
| 2022 precapillary PH without vasculopathy | none | 0/39 |
| PVR > 5 WU | 7% (3–12), only with vasculopathy | 7%; 8 with / 0 without |
| ΔBIC, low-DLCO phenotype | 78.8 | 78.8 |
| Low-DLCO component | 50; all with vasculopathy; 84% PH | same |
| Sensitivity / specificity for PH | 71% / 94% | 71% (42/59) / 94% |
| Mean AUC (FVC/DLCO, DLCO, FVC) | 0.75 / 0.74 / 0.65 | 0.745 / 0.738 / 0.654 |
| Fibrosis alone | 17.6 ± 0.5; r −0.92 and −0.92 | same |
| Variable vascular responses | SD 0.8, max 20.3; Spearman 0.47 and −0.15 | same |
| Alternative mechanism | PH 4% (0–10); r(FVC) −0.27 (−0.50 to −0.02) | same |
| Tree size (1,024 − 256 units), mPAP | −0.08 ± 0.41 mmHg | same |
| Morphometry, arterial / venous error factor | 1.14 / 1.16 | 1.135 / 1.156 |
| Healthy mPAP | 14.8 mmHg | 14.78 |

## Model checks

- **Physiological baseline:** FVC 3.8 L, cardiac output 5 L/min, mPAP 14.3 mmHg, PVR 1.26 WU, PaO₂ 97.5 mmHg, SaO₂ 97%.
- **PVR units:** (mPAP − left atrial pressure) / cardiac output, in Wood units; (14.3 − 8) / 5 = 1.26.
- **Womersley factor:** 1 for α → 0 and 2.78 at α = 10, guarded by a test (it was half its exact value before version 0.1.0).
- **Remodelling law:** the healthy tree is a rest point of the rule (residual 3×10⁻⁶).
- **Determinism:** every cohort is reproduced case by case from its seeds (checked to four decimals).
- **Statistics:** one routine (`src/ph_stats.py`) produces Table 1, the supplementary tables and `ph_cohort.py analyze`; percentile bootstrap over simulated cases; correlations averaged over 1,000 realisations of measurement noise.

## Assumptions examined

**Range of lung function.** Each simulated case starts from a healthy lung (FVC and DLCO at 100% of its own baseline). The single cross-sectional visit is drawn between 25% and 85% of healthy tissue remaining, so the cohorts describe established disease: FVC 42–80% and DLCO 11–65% at the visit, whereas earlier points of the same trajectories reach values close to 100% (supplementary figure S1, table S7). DLCO falls faster than FVC, as in IPF.

**Baseline in % predicted.** Expressing lung function relative to each case's own baseline ignores the normal spread of predicted values between individuals. Adding that spread (SD 13% for FVC and 15% for DLCO, correlation 0.5; 20 replicates) widens lung function to clinical ranges and leaves all seven calibration targets within the intervals of Table 1 (supplementary table S6).

**Altitude.** The model assumes sea-level conditions (baseline PaO₂ 97.5 mmHg).

## Open items

- Confirm references 20 (conference abstract) and 11, and the DLCO AUC of Joseph et al. against the full text.
- The length of the main pulmonary artery (9.05 cm) comes from a compilation of morphometric data; its diameter is cited from Singhal et al. (1973).

## Mathematical and methodological audit

### Invariants (now automated in `tests/test_invariants.py`)

| Check | Result |
|---|---|
| Flow conservation: unit flows add up to cardiac output; every bifurcation conserves flow | Pass (relative error < 1e-9) |
| mPAP = LAP + Q × equivalent resistance; PVR = (mPAP − LAP) / Q in Wood units | Pass |
| Poiseuille resistance and unit conversion to Wood units | Pass |
| Shared law: the mismatch between sensed shear and target decreases monotonically with radius, so the radius is unique | Pass |
| Healthy tree is a rest point of the remodelling rule; DLCO normalised to 1 at baseline | Pass |
| Womersley factor: 1 for α → 0, 2.78 at α = 10 | Pass (`tests/test_lupnet.py`) |

### Findings and their impact

1. **FVC plateau in advanced disease (resolved in 0.3.0).** FVC is computed from summed unit compliance. With honeycomb units at 60% of normal compliance, FVC could not fall much below 40% and PH prevalence at FVC < 50% was overestimated. Version 0.3.0 sets honeycomb compliance to that of fibrotic units (15%) and recalibrates the vasculopathy prevalence (70%); FVC now reaches 24% and PH prevalence at FVC < 50% is 48% (transplant candidates 46%).
2. **Phenotype definition used the vasculopathy labels (methodological).** The expected DLCO was fitted to cases without vasculopathy, information not available in patients. Fitting it to all cases (ordinary or trimmed least squares) keeps the mixture bimodal (ΔBIC 78–79) with sensitivity 70–73% and specificity 94% for PH (table S9). **Resolved:** the phenotype can be identified from FVC and DLCO alone.
3. **Fixed left atrial pressure (8 mmHg).** Varying wedge pressure leaves PH ≥ 25 mmHg, the correlations and the DLCO ratio essentially unchanged, but a 1 mmHg shift moves many cases across 20 mmHg, because one third of cases lie between 18 and 20 mmHg (table S8). Relevant for the 2022 definition.
4. **Measurement noise applied to correlations but not to pressure categories (methodological inconsistency).** Adding noise to mPAP moves the categories slightly closer to the benchmarks (table S11); conclusions unchanged.
5. **Diffusing capacity.** The membrane and capillary components are given equal weight in a healthy unit (Roughton–Forster), and θ does not depend on haemoglobin or inspired oxygen. Adequate for relative changes at sea level; altitude or anaemia would require θ(PO₂, Hb).
6. **Vital capacity as summed compliance.** A linear simplification of the pressure–volume curve; it reproduces relative FVC decline but not absolute volumes in severe restriction.
7. **Calibration and validation.** The five vasculopathy parameters were chosen against the mild-to-moderate targets; bootstrap intervals reflect sampling of simulated cases, not parameter uncertainty.
8. **Morphometry.** Flow per order assumes equal division among elements (Q/N); supernumerary branches and the single-lung source are not represented.
