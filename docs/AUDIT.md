# Reproducibility audit

Manuscript: *A vascular contribution to pulmonary hypertension in idiopathic pulmonary fibrosis* (in preparation; not part of this repository).
Code: this repository, version 0.4.0 (audit of 0.3.7 and of the additions up to 0.4.0; the numbers below are unchanged, except the pressure–flow response, finding 21).

## Method

1. The repository was cloned into a clean environment. The tests (16/16, including the regression test of the published numbers) and `run_all.sh` were run in full: arterial and venous morphometry, Figures 1–4 (Figure 2 including panels f–h), Table 1, supplementary tables, supplementary figure S1 and the analysis of the five cohorts.
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

## Non-invasive markers (Figure 5)

| Quantity | Manuscript | Recomputed |
|---|---|---|
| AUC, relevant vasculopathy: FVC + DLCO / + SpO2 / + PaO2 / + A–a / + echo / + echo and A–a | 0.69 / 0.72 / 0.74 / 0.77 / 0.84 / 0.85 | 0.687 / 0.721 / 0.738 / 0.768 / 0.837 / 0.853 |
| AUC, PH (same order) | 0.78 / 0.84 / 0.87 / 0.87 / 0.92 / 0.94 | 0.776 / 0.836 / 0.869 / 0.874 / 0.917 / 0.937 |
| Partial r(PaO2, lesion severity given FVC and DLCO) | −0.38 | −0.38 |
| Individual mPAP error: FVC + DLCO / + A–a / + echo | 4.2 / 3.0 / 2.2 mmHg | 4.21 / 3.03 / 2.21 |
| Coverage of the 90% intervals | 89–92% | 89–92% |
| PH by FVC range, 34–66% (table S11) | 31% (25–38) | 31% (25–38) |

All figures and tables (Figures 1–5, S1; Table 1; tables S2, S4–S12, S13–S15) regenerated from a clean clone are value-identical to those in the manuscript (Figure 2f and 2g differ visually since 0.3.8, findings 9 and 10 below); `src/manuscript_numbers.py` recomputes every other number quoted in the text.

## Model checks

- **Physiological baseline:** FVC 3.8 L, cardiac output 5 L/min, mPAP 14.3 mmHg, PVR 1.26 WU, PaO₂ 97.5 mmHg, SaO₂ 97%.
- **PVR units:** (mPAP − left atrial pressure) / cardiac output, in Wood units; (14.3 − 8) / 5 = 1.26.
- **Womersley factor:** 1 for α → 0 and 2.78 at α = 10, guarded by a test (it was half its exact value before version 0.1.0).
- **Remodelling law:** the healthy tree is a rest point of the rule (residual 3×10⁻⁶).
- **Determinism:** every cohort is reproduced case by case from its seeds (checked to four decimals).
- **Statistics:** one routine (`src/ph_stats.py`) produces Table 1, the supplementary tables and `ph_cohort.py analyze`; percentile bootstrap over simulated cases; correlations averaged over 1,000 realisations of measurement noise.

## Assumptions examined

**Range of lung function.** Each simulated case starts from a healthy lung (FVC and DLCO at 100% of its own baseline). The single cross-sectional visit is drawn between 25% and 85% of healthy tissue remaining, so the cohorts describe established disease: FVC 42–80% and DLCO 11–65% at the visit, whereas earlier points of the same trajectories reach values close to 100% (supplementary figure S1, table S7). DLCO falls faster than FVC, as in IPF.

**Baseline in % predicted.** Expressing lung function relative to each case's own baseline ignores the normal spread of predicted values between individuals. Adding that spread (SD 13% for FVC and 15% for DLCO, correlation 0.5; 20 replicates) widens lung function to clinical ranges and leaves all seven calibration targets within the intervals of Table 1 (supplementary table S8).

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
2. **Phenotype definition used the vasculopathy labels (methodological).** The expected DLCO was fitted to cases without vasculopathy, information not available in patients. Fitting it to all cases (ordinary or trimmed least squares) keeps the mixture bimodal (ΔBIC 78–79) with sensitivity 70–73% and specificity 94% for PH (table S10). **Resolved:** the phenotype can be identified from FVC and DLCO alone.
3. **Fixed left atrial pressure (8 mmHg).** Varying wedge pressure leaves PH ≥ 25 mmHg, the correlations and the DLCO ratio essentially unchanged, but a 1 mmHg shift moves many cases across 20 mmHg, because one third of cases lie between 18 and 20 mmHg (table S9). Relevant for the 2022 definition.
4. **Measurement noise applied to correlations but not to pressure categories (methodological inconsistency).** Adding noise to mPAP moves the categories slightly closer to the benchmarks (table S12); conclusions unchanged.
5. **Diffusing capacity.** The cohorts use the product of membrane conductance and capillary volume (`rf = False`); the Roughton–Forster series form exists in the code (`rf = True`) but was not used. Capillary loss from the lesion (`cap_frac`) therefore enters DLCO linearly. θ does not depend on haemoglobin or inspired oxygen; altitude or anaemia would require θ(PO₂, Hb). (Corrected in 0.3.8: earlier versions of this audit described the cohorts as Roughton–Forster.)
6. **Vital capacity as summed compliance.** A linear simplification of the pressure–volume curve; it reproduces relative FVC decline but not absolute volumes in severe restriction.
7. **Calibration and validation.** The five vasculopathy parameters were chosen against the mild-to-moderate targets; bootstrap intervals reflect sampling of simulated cases, not parameter uncertainty.
8. **Morphometry.** Flow per order assumes equal division among elements (Q/N); supernumerary branches and the single-lung source are not represented. The morphometric fits use a free τ₀ (least squares over 16 orders, viscosity per order); the lung model fixes τ₀ on the main pulmonary artery with μ = 3.5 mPa s.

### Findings of the 0.3.7 audit (v0.3.8; `src/audit_extra.py`, tables S18–S21, `tests/test_audit_extra.py`)

9. **Figure 2f used the wrong clock (bug, fixed).** The lesion severity along the trajectories was computed from `t × ypu`, which omits the 6.71 years per model unit; the curves ran 6.7 times too slowly. The stored `years` column is now used. No number quoted in the text changes; the panel does: the independent lesion is within a few percent of its final value before 20% of the lung is scarred.
10. **Ageing inside the "fibrosis" contribution.** The calibrated parameters age the vessels after the healthy baseline (`age_M = 0.3`). The aged lung with no fibrosis has an mPAP of 16.5 mmHg (healthy 14.3, G = 8; table S19). Of the published contributions of fibrosis in Figure 2g (3.5, 2.7, 2.2 mmHg above 14.8), 1.7 mmHg is ageing; fibrosis alone adds 1.8, 1.0 and 0.4 mmHg. Figure 2g now draws the ageing offset as a reference line and the script prints both decompositions. The manuscript text (abstract and "separate clocks") should state the contributions above the aged lung.
11. **Lesion saturated at the visit.** With τ = 2 years and t₀ ~ U(0, 10), the lesion is at its final severity in nearly all cases at the visit (table S20), so τ and t₀ are not identified by the cohorts; the sensitivity rows for τ (table S17) move little for this reason.
12. **Phenotype fragile to baseline spread.** The low-DLCO phenotype (ΔBIC 78.8, sensitivity 71%) depends on the model's DLCO–FVC correlation of 0.97. Under the baseline spread of table S8 plus measurement noise, ΔBIC has a median of 7, a second component is preferred in 40% of replicates, and sensitivity is 34% (table S21).
13. **Healthy pressure as output.** 47% of the healthy pressure drop is set by the law; the capillary and venular shares (3.5 of 6.6 mmHg) come from `dP0 × split` (table S18).
14. **Fixed ventilation.** VA = 4.2 L/min in every state; 49 of 200 calibrated cases have PaCO₂ > 45 mmHg at the visit. Relevant for PaO₂ and the A–a gradient of Figure 5.
15. **Runs ending before the drawn stage.** 7 of 60 cases of the advanced cohort reach `t_max = 3.0` before the drawn fraction of healthy tissue; `ph_cohort.py run` now reports such cases and takes `--t_max`.
16. **Shared `Params` object.** `Lung._baseline` and the exacerbation branch of `sim.run` alter fields of the `Params` instance and restore them; both restores are now in `finally` blocks.
17. **Tag prefixes (bug, fixed in 0.3.11).** `cohort_io.cases("big2")` matched `big2_*.tsv` and would have read the per-case files of `big2_uip` or `big2_cdmx` as the manuscript cohort. It now reads only `tag_<digits>.tsv`; guarded by `tests/test_cohort_io.py`.
18. **Non-viable visits (0.3.14).** Five runs of the field cohorts end with the oxygen demand not met (SaO₂ down to 0.27). `cohort_io.viable` and `final_visits` use the last viable visit and count the run as a death; the manuscript cohorts have none, so their numbers are unchanged.
19. **Severity as applied (0.3.14).** `vasc_A` is the drawn severity; with `--a_scale` the applied one is `vasc_A_scaled`, and each visit now records the lesion in place (`lesion_mean`, `lesion_healthy`). The tool uses the applied severity.
20. **Two-threshold rule of the tool.** With probabilities ≥ 0.7 and < 0.3 the tool classifies about two thirds of the 43 patients with about 90% correct, but the counts move between cross-validation partitions and references (25 of 27, 26 of 30). The rule is not stable at this sample size and is not proposed for use.
21. **Pressure–flow response reported at the wrong output (error in the analysis of manuscript v21, corrected in 0.4.0).** Figure 6 and table S16 were described as a cardiac output of 10 L/min. The numbers of v21 (mPAP 21.3 against 18.1 mmHg without lesion, slope 0.63 mmHg per L/min) are reproduced only when 10 L/min are demanded of the resting right ventricle, which delivers about 7.3 L/min; the slope had been divided by 5 L/min. With the flow imposed, the same cases reach about 23.6 mmHg and the slope is about 1.26 without lesion and 2.5 to 3.5 with a moderate lesion. `src/exercise_response.py` imposes the flow, records 7.5 and 10 L/min, and keeps the earlier quantities (`mPAP_rv`, `Q_rv`) so that the difference can be traced. The sentences of the manuscript on the proportion of cases unmasked by flow and on the low absolute slope have to be rewritten from the new table.
22. **Results of the manuscript without code (closed in 0.4.0).** Figure 6, table S16 and table S6 (256 against 1,024 units) had no script in the repository; the supplement cited `src/exercise_response.py` and `results/figures/pressure_flow_cases.tsv`, which were not versioned. Both are now scripts with archived per-case results.
23. **Table numbers (closed in 0.4.0).** The repository numbered the supplementary tables on its own; they now carry the numbers of the supplement (correspondence in the README).
24. **Provenance of the cohorts (closed in 0.4.0).** Cohort files did not record how they were generated, and the lesion time constant was written by hand in five scripts. `ph_cohort.py run` now stores tau, cap, p_vasc, the lesion distribution, seeds and the code version in every file, and readers take the time constant from the file (`cohort_io.lesion_tau`; 2.0 for cohorts written before 0.4.0).
25. **Partial re-runs (closed in 0.4.0).** A single per-case file next to a packed cohort made `cohort_io.cases` read that one case as the whole cohort. It now uses the packed cohort unless the per-case files are at least as many, and warns.
26. **Web calculator.** `index.html` is the calculator of v0.3.2 (FVC, DLCO, echocardiographic pressure, A–a gradient, SpO₂; simulated visits as reference, validated on simulated cases only). It does not include the measurement models fitted on patients nor the pulmonary artery diameter, and its blood-gas inputs depend on the fixed ventilation of finding 14. The page states that it is a research tool; it has not been replaced.
