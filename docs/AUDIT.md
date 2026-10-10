# Reproducibility audit

Manuscript: *A vascular contribution to pulmonary hypertension in idiopathic pulmonary fibrosis* (in preparation; not part of this repository).
Code: this repository, version 0.2.2.

## Method

1. The repository was cloned into a clean environment. The tests (4/4) and `run_all.sh` were run in full: arterial and venous morphometry, Figures 1–4, Table 1, supplementary tables, supplementary figure S1 and the analysis of the five cohorts.
2. An independent script, separate from the repository code, recomputed every number quoted in the manuscript from the stored simulations (`results/ph_cohorts/*.tsv.gz`).
3. Model physics and units were checked against the baseline values of individual simulated cases.

## Reproducibility

Every number in the manuscript is reproduced exactly from the stored simulations.

| Quantity | Manuscript | Recomputed |
|---|---|---|
| Calibrated cases / with FVC 50–90% | 200 / 162 | 200 / 162 |
| Realised prevalence of vasculopathy (calibrated / alternative) | 49% / 38% | 49% / 38% |
| mPAP < 20 mmHg | 54% (46–62) | 53.7% (45.7–61.7) |
| mPAP 20–25 mmHg | 29% (22–36) | 29.0% (22.2–36.4) |
| mPAP ≥ 25 mmHg | 17% (12–23) | 17.3% (11.7–23.5) |
| mPAP ≥ 35 mmHg | 4% (1–7) | 4.3% (1.2–7.4) |
| r(mPAP, DLCO) | −0.35 (−0.47 to −0.22) | −0.353 (−0.468 to −0.225) |
| r(mPAP, FVC) | −0.13 (−0.27 to 0.02) | −0.126 (−0.267 to 0.019) |
| DLCO ratio, PH / no PH | 0.78 (0.68–0.88) | 0.779 (0.677–0.882) |
| DLCO ratio at the 21 mmHg threshold | 0.76 | 0.764 |
| No vasculopathy: mPAP | 18.6 ± 1.3, max 23.2 | 18.6 ± 1.3, 23.2 |
| Vasculopathy: mPAP, % PH | 25.3 ± 7.0, 37% | 25.3 ± 7.0, 37% |
| PH with FVC < 50%: calibrated / advanced / pooled | 47% / 71% / 61% | 47% / 71% / 61% |
| Severe PH (mPAP > 40 mmHg), pooled | 6% (1–11) | 6% (1–11) |
| FVC with / without PH, FVC < 50% | 47 / 47% | 47.0 / 47.4% |
| 2022 precapillary PH without vasculopathy | 13%, PVR ≤ 3.0 WU | 11/86, 3.03 WU |
| PVR > 5 WU, with / without vasculopathy | only with | 9 / 0 |
| ΔBIC, low-DLCO phenotype | 17.6 | 17.6 |
| Low-DLCO component | 27; 85% vasculopathy; 89% PH | same |
| Sensitivity / specificity for PH | 52% / 98% | 52% / 98% |
| Mean AUC (FVC/DLCO, DLCO, FVC) | 0.75 / 0.73 / 0.63 | 0.752 / 0.733 / 0.627 |
| Fibrosis alone | 18.2 ± 1.1; r −0.88 and −0.93 | same |
| Variable vascular responses | SD 1.1; Spearman 0.06 and 0.00 | 1.14; 0.06 and 0.00 |
| Alternative mechanism (linked to fibrosis) | −0.26, −0.47, 0.51, 66% | same |
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

1. **FVC plateau in advanced disease (model assumption).** FVC is computed from summed unit compliance; honeycomb units keep 60% of normal compliance (`hc_C = 0.60`), so FVC cannot fall much below 40%. This explains the overestimate of PH prevalence at FVC < 50%. A first-order correction to `hc_C = 0.15` extends FVC to 27% and brings PH prevalence at FVC < 50% to 47% (transplant candidates 46%), with little change in mild-to-moderate disease (table S10). **Recommended fix:** set `hc_C` to the value of fibrotic units, recalibrate the FVC decline and re-run the cohorts.
2. **Phenotype definition used the vasculopathy labels (methodological).** The expected DLCO was fitted to cases without vasculopathy, information not available in patients. Fitting it to all cases (ordinary or trimmed least squares) keeps the mixture bimodal (ΔBIC 26–28) with sensitivity 59–61% and specificity 94–96% for PH (table S9). **Resolved:** the phenotype can be identified from FVC and DLCO alone.
3. **Fixed left atrial pressure (8 mmHg).** Varying wedge pressure leaves PH ≥ 25 mmHg, the correlations and the DLCO ratio essentially unchanged, but a 1 mmHg shift moves many cases across 20 mmHg, because one third of cases lie between 18 and 20 mmHg (table S8). Relevant for the 2022 definition.
4. **Measurement noise applied to correlations but not to pressure categories (methodological inconsistency).** Adding noise to mPAP moves the categories slightly closer to the benchmarks (table S11); conclusions unchanged.
5. **Diffusing capacity.** The membrane and capillary components are given equal weight in a healthy unit (Roughton–Forster), and θ does not depend on haemoglobin or inspired oxygen. Adequate for relative changes at sea level; altitude or anaemia would require θ(PO₂, Hb).
6. **Vital capacity as summed compliance.** A linear simplification of the pressure–volume curve; it reproduces relative FVC decline but not absolute volumes in severe restriction.
7. **Calibration and validation.** The five vasculopathy parameters were chosen against the mild-to-moderate targets; bootstrap intervals reflect sampling of simulated cases, not parameter uncertainty.
8. **Morphometry.** Flow per order assumes equal division among elements (Q/N); supernumerary branches and the single-lung source are not represented.
