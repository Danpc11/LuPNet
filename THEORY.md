# LuPNet: theory

This document describes the model, the simulated cohorts and the individual estimate, and gives every number reported in the README. Parameter names refer to `src/lupnet/params.py`; calibrated values for aged IPF are in `data/calibrated_ipf_params.json`.

## 1. Structure

The lung is a set of N gas-exchange units (1,024 by default, `G = 10` generations of a symmetric tree; the in silico cohorts of the manuscript use 256 units, `G = 8`, which give the same results within 0.1 mmHg on average; or an anatomical 3D tree) perfused by an arterial tree and drained by a venous tree. Each unit lumps a region of acini and its capillary sheet. A unit is healthy, fibrotic or honeycombed; fibrotic and honeycomb units have lower compliance, volume, capillary bed and membrane conductance.

**Breathing and mechanics.** Ventilation is distributed by unit compliance; healthy units can lose surfactant function, collapse and recruit cyclically. The forced vital capacity is computed from the volume–pressure behaviour of all units.

**Gas exchange.** Each unit has a membrane conductance and a capillary volume; its diffusing capacity is their product (`rf = False`, the setting of every cohort of the manuscript; `rf = True` gives the Roughton–Forster series form, which the cohorts do not use), DLCO is the sum over ventilated units, and arterial oxygen saturation follows from ventilation–perfusion matching. Alveolar ventilation is held at 4.2 L/min; PaCO₂ is a readout (0.863 VCO₂ / effective VA), not a regulated variable, so cases with large alveolar dead space are hypercapnic (49 of 200 calibrated cases above 45 mmHg, 14 above 50).

**Perfusion.** Arterial and venous conduits follow Poiseuille's law with viscosity μ; the capillary sheet follows sheet-flow mechanics with gravity, alveolar and vascular pressures; hypoxic pulmonary vasoconstriction acts on the terminal arterioles.

**Disease.** Fibrosis spreads between units with a hazard that depends on surfactant loss, collapse and neighbourhood (`beta = 0` in the calibrated set: no strain feedback); the vessels age after the healthy baseline (`age_M = 0.3`, `age_sigma = 0.3`: rest point of the remodelling rule at 1.3 times the target shear), which by itself raises mPAP from 14.3 to 16.5 mmHg in the 256-unit lung with no fibrosis (table S19); time is in model units (6.71 years per unit in the calibrated set, `years_per_unit`). The left atrial pressure is fixed at 8 mmHg, so the model represents pre-capillary pulmonary hypertension only. The right ventricle adapts to afterload. Acute exacerbations, epithelial senescence, the septal venous tree and the anatomical tree are options of the code that no cohort uses.

## 2. The shared target law for the pulmonary vessels

### 2.1 Target

For a tree that minimizes viscous dissipation at a fixed maintenance cost, with the cost of a vessel proportional to its mass to the power b, the shear stress of every vessel at the optimum is

$$\tau^*=\tau_0\,r^{\,b-1}\,\ell^{\,(b-1)/2}$$

where r is the lumen radius and ℓ the segment length. b = 1 is Murray's law (uniform shear). The law, the exponent b and the weight w₁ are derived and tested across species in Pérez-Calixto et al., *An invariant endothelial flow signal across species*, bioRxiv 2026 (https://doi.org/10.64898/2026.09.29.755319; code at https://github.com/Danpc11/InFlow). The local rule

$$\frac{d\ln r}{dt}=\kappa\left(\frac{S}{\tau^*}-1\right)$$

lowers the total cost for steady flow. With `law="shared"`, one τ₀ is used for every vessel; it is fixed so that the main pulmonary artery has the radius `r_pa`.

### 2.2 Sensed shear

$$S=\bar\tau\left[1+w_1\,\varphi_1\,G(\alpha)\right],\qquad \bar\tau=\frac{4\mu Q}{\pi r^3},\qquad \alpha=r\sqrt{\omega\rho/\mu},\ \ \omega=2\pi\,\mathrm{HR}/60$$

G is the exact Womersley ratio of oscillatory to quasi-steady wall shear in a rigid tube (1 for α → 0, ≈ α/4 for large α), w₁ = 2.865 is the weight of the first harmonic fitted to endothelial NF-κB responses (Feaver et al. 2013), and φ₁ the first-harmonic amplitude of flow relative to its mean (effective values 0.8 in arteries, 0.4 in veins).

### 2.3 The terminal element

The explicit tree stops at generation G. The terminal arteriole of each unit stands for the sub-tree below it; its resistance is that of the same law continued generation by generation down to 10 µm radius.

### 2.4 Tests against human morphometry

Huang et al. (1996) give, for 15 orders of human pulmonary arteries, the number of elements, diameter and length; we add the main pulmonary artery as order 16 (diameter 3.0 cm, Singhal et al. 1973; length 9.05 cm). Each order carries Q/N of a 5 L/min cardiac output.

| Model | Best b | Diameter error (16 orders) |
|---|---|---|
| Pulsatile sensing (φ₁ = 0.8) | 0.90 | factor 1.14 |
| Mean shear only | 0.825 | factor 1.14 |
| Murray (b = 1, mean only) | – | factor 1.26 |

With b = 0.675, the exponent fitted to systemic arteries across mammals, the predicted main pulmonary artery is twice its real size. The arterial law (b = 0.90) applied unchanged to the 15 venous orders predicts their diameters within a factor 1.16 (venous φ₁ = 0.4). In these morphometric tests τ₀ is a free scale fitted by least squares in log D over the 16 orders, with the apparent viscosity of each order. Integrated in the lung model (`law="shared"`, b = 0.925, inside the flat region of the error profile: factor 1.136 vs 1.135 at b = 0.90) τ₀ is instead fixed by the radius of the main pulmonary artery and μ = 3.5 mPa s; the explicit arterial generations then match Huang with a geometric error factor of 1.06 (largest deviation 8.2%), and the healthy mean pulmonary artery pressure is 14.8 mmHg with 1,024 units and default parameters, 14.6 with 256 units and default parameters, and 14.3 with the 256 units and calibrated parameters of the cohorts, which is the healthy lung every simulated case starts from (per-vessel law, calibrated: 14.0). The law sets the arterial conduits, the venous conduits and the terminal sub-tree, 47% of the healthy pressure drop (0.32 + 0.67 + 2.13 mmHg of 6.6, 256 units, default parameters); the capillary sheet (2.10) and the venules (1.39) keep the resistances of the reference calibration, `dP0` × `split` (table S18). The healthy mPAP is therefore an output of the law for its arterial half and an input for the rest.

## 3. Vascular lesion

A non-remodellable resistance in series with the terminal arteriole of unit i, with capillary loss:

$$R_{\mathrm{les},i}=R_{ta,0}\,A\,w_i\left(1-e^{-(t_\mathrm{yr}+t_0)/\tau}\right),\qquad c_{\mathrm{cap},i}=\frac{1}{1+f_\mathrm{cap}\,A\,w_i\left(1-e^{-(t_\mathrm{yr}+t_0)/\tau}\right)}$$

A is the severity (`vasc_A`, in units of the normal arteriolar resistance), t₀ the years of lesion before the start (`t_vasc0`), τ its time constant (`vasc_tau`), f_cap the capillary loss per unit of severity (`cap_frac`) and w_i the weight of the unit, set by where the lesion sits (`vasc_where`, `sim.lesion_weights`):

| `vasc_where` | w_i |
|---|---|
| `diffuse` | 1 in every unit |
| `scar` | 1 in fibrotic and honeycomb units, 0 in preserved ones |
| `field` | 1 in the scar; in a preserved unit, the fraction of its neighbours that are scar |
| `mixed` | m + (1 − m) × field, with m = `vasc_mix` |

With A = 0 the model is unchanged.

**Simulated cases** (`src/ph_cohort.py`). With probability 0.7 a case carries a lesion, A ~ lognormal(median 1.5, σ = 1.1), t₀ ~ U(0, 10) years, τ = 2 years, f_cap = 0.05, diffuse; the speed of fibrosis varies between cases (years per model unit × lognormal, σ = 0.5); one visit at a random stage (25–85% of preserved units in the calibrated cohort, 12–45% in the advanced one). Every cohort file stores these values, the seeds and the code version.

At the visit the lesion has reached its final severity in almost every case (1 − e^{−(t+t₀)/τ}: minimum 0.60, 5th percentile 0.96, median 1.00 in the calibrated cohort; table S20). The cohorts identify the prevalence of the lesion and the distribution of A, not τ or t₀.

## 4. Pulmonary hypertension in the simulated cohorts

### 4.1 Fibrosis without lesion

In 60 simulated cases without lesion, mPAP in mild-to-moderate disease (FVC 50–90%) is 17.6 ± 0.5 mmHg; the healthy lung of the cohorts is at 14.3 mmHg and the aged lung before any fibrosis at 16.5 (table S19); r(mPAP, FVC) = −0.92 and r(mPAP, DLCO) = −0.92; no case reaches 25 mmHg. Varying hypoxic vasoconstriction or pressure-driven maladaptive remodelling between cases spreads mPAP by 0.8 mmHg (maximum 20.3 mmHg).

### 4.2 Pressure surface

In every cohort mPAP is a function of two variables, the fraction of lung lost and the effective lesion: mPAP = F(lost, A), fitted to 11,022 simulated visits with R² = 0.998 (0.996–0.999 by cohort).

| Lung lost | A = 0 | A = 1 | A = 2 | A = 4 | A = 8 |
|---|---|---|---|---|---|
| 5% | 16.0 | 18.2 | 20.2 | 24.3 | 33.3 |
| 30% | 17.6 | 20.1 | 23.0 | 28.0 | 36.7 |
| 60% | 19.7 | 22.4 | 25.8 | 30.7 | 42.3 |

mPAP in mmHg. Without lesion, losing 60% of the lung adds 3.7 mmHg; a lesion A = 2 adds 4.2 mmHg at 5% of lung lost and 6.1 mmHg at 60%. The effect of fibrosis is a partial dependence: small alone, and a multiplier of the lesion, because the same cardiac output crosses fewer units.

### 4.3 Diffuse lesion: benchmarks

200 simulated cases, 113 with FVC 50–90% (README table). The seven benchmarks (ARTEMIS-IPF, Raghu et al. 2015; Zisman et al. 2007; the DLCO ratio of a Japanese cohort) lie inside the 95% bootstrap intervals; correlations include measurement noise (mPAP 2.5 mmHg, FVC and DLCO 5 points).

Not fitted:

| | Model | Reference |
|---|---|---|
| PH with FVC < 50%, pooled calibrated and advanced cohorts (n = 147) | 48% (39–56) | 46% (Shorr et al. 2007, 2,525 transplant candidates) |
| mPAP > 40 mmHg, FVC < 50%, pooled (n = 147) | 8% (4–13) | ≈ 9% |
| FVC, PH vs no PH, FVC < 50% | 35.7 vs 39.1% | 48.4 vs 51.4% |
| AUC of FVC/DLCO for mPAP ≥ 25 (mean of 1,000 noise realisations) | 0.75 | 0.69–0.74 in ILD cohorts |
| AUC of DLCO alone / FVC alone | 0.74 / 0.65 | ≈ 0.80 / no discrimination |

**Lesion scaled by the extent of fibrosis** (`vasc_mode="fibrosis"`, 100 cases): PH 4% and r(mPAP, FVC) = −0.27 with FVC 50–90%, both outside the intervals of the benchmarks.

**Sensitivity** (`src/sensitivity.py`, table S17, figure S2; 200 cases each, FVC 50–90%, n = 113):

| | Cases with lesion | PH | mPAP ≥ 35 | r(mPAP, FVC) | r(mPAP, DLCO) | DLCO ratio | PH, FVC < 50% |
|---|---|---|---|---|---|---|---|
| Calibrated (τ = 2 years, 70%) | 65% | 18.6% | 5.3% | −0.08 | −0.38 | 0.81 | 44% |
| τ = 1 year | 65% | 18.6% | 6.2% | −0.07 | −0.37 | 0.81 | 44% |
| τ = 3 years | 65% | 16.8% | 5.3% | −0.10 | −0.38 | 0.81 | 44% |
| Lesion in 50% of cases | 48% | 14.2% | 4.4% | −0.05 | −0.35 | 0.81 | 29% |
| Lesion in 85% of cases | 81% | 21.2% | 5.3% | −0.05 | −0.35 | 0.85 | 49% |

**Tree size** (`src/tree_size.py`, table S6). The same 12 cases with 1,024 and 256 units differ by −0.08 ± 0.41 mmHg in mPAP (r = 0.997), −0.6 ± 1.7 points in FVC and −0.2 ± 0.8 points in DLCO; 11 of 12 keep the same PH status.

### 4.4 Low-DLCO phenotype

The deficit of DLCO relative to the value expected for the FVC is bimodal (two-component mixture, ΔBIC = 78.8). The low-DLCO component has 50 of 200 cases, all with lesion, and 84% have PH; in the rest, 60% carry a lesion and 11% have PH. The separation rests on the DLCO–FVC relation of the model (r = 0.97 at the visit; about 0.5–0.7 in patients): with the spread of the healthy baseline of table S8 and measurement noise, ΔBIC falls to a median of 7 (range −12 to 46 over 20 replicates), a second component is preferred in 40% of replicates, and the sensitivity for PH falls from 71% to 34% (table S21).

## 5. Where the lesion sits, histological pattern and altitude

### 5.1 Where the lesion sits

Cohorts with the fibrosis course and the lesion draws of the calibrated and advanced cohorts, differing only in `vasc_where` (`src/lesion_location.py`, tables S25–S26). `--a_scale` multiplies A so that each location has the PH prevalence of the diffuse cohort in mild-to-moderate disease (field: 5.8; mixed: `--where mixed --mix 0.25 --a_scale 2.0`, that is 0.5 A in every unit plus 1.5 A × field weight).

| Lesion | Benchmarks met | r(mPAP, FVC) | PH, FVC < 50% | mPAP > 40, FVC < 50% | PH, FVC ≥ 70% |
|---|---|---|---|---|---|
| Diffuse | 7 of 7 | −0.08 | 48% (39–56) | 8% (4–13) | 23% |
| Mixed (0.5 A diffuse + 1.5 A × field weight) | 7 of 7 | −0.17 | 54% (46–63) | 16% (11–22) | 14% |
| Field | 4 of 7 | −0.34 | 68% (60–76) | 31% (23–38) | 0% |
| Scar only | 1 of 7 | −0.19 | 0% | 0% | 0% |
| Patients | | ≈ 0 | 46% | ≈ 9% | |

The scar-only lesion gives no PH in 200 cases (maximum 23.9 mmHg): the capillary bed of the scar is largely gone and flow goes elsewhere. Five runs of the field cohorts end with the oxygen demand not met; their last viable visit is used (`cohort_io.viable`).

### 5.2 Histological pattern

`ph_cohort.py run --pattern`. The calibrated set is patchy without subpleural predominance (units lost: 0.46 subpleural, 0.44 core; honeycomb 22% of lost units).

| | `uip` | `nsip` |
|---|---|---|
| Parameters | `sub_gain = 3`, `kn = 1.5`, `kn_s = 3`, `lam_hc = 0.5` | `kn = 0`, `kn_s = 0`, `at2_sigma = 0.3`, `lam_hc = 0.05`, `fib_C = 0.30`, `nbr_D = 1` |
| Units lost, subpleural / core | 0.60 / 0.25 | 0.47 / 0.38 |
| Units lost, base / apex | 0.61 / 0.41 | 0.68 / 0.30 |
| Honeycomb, share of lost units | 13% | 2% |
| Years to the visit | 3.9 | 6.9 |
| PH, FVC 50–90% | 17% | 17% |
| Benchmarks met | 7 of 7 | 6 of 7 (DLCO ratio 0.85) |
| PH, FVC < 50% | 45% | 35% |

`sub_gain` multiplies surfactant decline and the direct fibrosis hazard by 1 + sub_gain (1.5 − y) in the two outer rings of units, y = 0 at the base (`sim.subpleural_factor`). The pressure surfaces of the two patterns differ by 1 mmHg or less.

### 5.3 Altitude

`--PB 585` (Mexico City). Hypoxic vasoconstriction is normalised to the alveolar PO₂ of each lung at baseline, so a healthy lung has the same mPAP at any barometric pressure. In disease, pressures differ by less than 1 mmHg from sea level (PH 16% against 19%); SaO₂ at the visit is 85% against 94%.

## 6. Pressure–flow response

`src/exercise_response.py`, Figure 6, table S16. The first 70 cases of the calibrated cohort (47 with lesion, 23 without) are re-equilibrated at imposed cardiac outputs of 7.5 and 10 L/min with the vessels held at the calibre they had reached; ventilation and oxygen uptake scale with the flow. Slope = (mPAP at 10 L/min − resting mPAP) / (10 − resting output).

| | No lesion (n = 23) | Lesion (n = 47) |
|---|---|---|
| Resting mPAP, mmHg | 18.1 (17.5–19.3) | 22.3 (20.7–25.8) |
| mPAP at 7.5 L/min | 21.7 (21.0–23.3) | 27.5 (25.4–32.4) |
| mPAP at 10 L/min | 24.8 (23.9–26.6) | 32.2 (29.4–38.4) |
| Slope, mmHg per L/min | 1.33 (1.28–1.44) | 1.96 (1.72–2.53) |
| Resting mPAP < 25 that reach 25 at 7.5 L/min | 1 of 23 (4%; 0–13) | 23 of 32 (72%; 56–84) |
| Resting mPAP < 25 that reach 25 at 10 L/min | 10 of 23 (43%; 22–65) | 32 of 32 |
| Slope > 3 | 0 of 23 | 9 of 47 (19%; 9–30) |

Median (interquartile range) and proportion (95% bootstrap interval). The resting right ventricle of the model delivers 7.0–7.3 L/min when 10 are demanded; the flow is therefore imposed.

## 7. Non-invasive estimate

### 7.1 Markers in the simulated cohorts

`src/noninvasive.py`, tables S13–S15 (figure 5 of manuscript v21). The trajectories give resting SpO₂, arterial blood gases and the systolic pressure that echocardiography would estimate (mPAP = 0.61 sPAP + 2 mmHg, noise 8 mmHg). With measurement noise, the area under the ROC curve for a relevant lesion (effective severity ≥ 1) is 0.69 with FVC and DLCO, 0.72 with SpO₂, 0.74 with PaO₂, 0.77 with the alveolar–arterial gradient, 0.84 with echocardiographic pressure and 0.85 with both; for PH, 0.78, 0.84, 0.87, 0.87, 0.92 and 0.94. With an emulator of the simulator, the error of the estimated mPAP in held-out simulated cases falls from 4.2 to 2.2 mmHg. The web calculator (`index.html`) uses the simulated visits as reference in the same way.

### 7.2 Individual estimate in patients

`src/noninvasive_tool.py`, tables S22–S24. For one patient, the posterior over the simulated visits is

$$p(\text{visit}\mid\text{data})\;\propto\;\pi(\text{visit})\;p(\mathrm{FVC},\mathrm{DLCO}\mid\text{visit})\prod_k p(x_k\mid \mathrm{mPAP})$$

- FVC and DLCO are required; their likelihood has the measurement noise (5 points) and the spread of the healthy baseline (13% and 15%).
- Optional inputs x_k, in any combination: echocardiographic sPAP, pulmonary artery diameter on CT, log BNP, right ventricular fractional area change, right atrial area. Each is linear in mPAP with Gaussian noise, fitted on catheterised patients.
- The prior π gives the simulated visits the distribution of mPAP of the population in five classes (`--prior severity`). The reference is the cohort of the patient's tomographic pattern (`--reference pattern`).
- Output: probability of PH, mPAP with a 90% interval, PVR, and the parts of the pressure due to the aged lung, the lung lost and the lesion, from the surface F.
- SpO₂ and blood gases are not inputs. Catheterisation is used only to fit the measurement models and to score the estimate.

Patients: 47 adults with fibrosing interstitial lung disease and right heart catheterisation at one centre in Mexico City, referred for suspected PH; 43 with FVC, DLCO and mPAP (22 NSIP, 9 UIP, 12 other patterns; 31 with mPAP ≥ 25 mmHg). r(mPAP, FVC) = 0.05 (n = 47) and r(mPAP, DLCO) = −0.28 (n = 43); 6 of 9 patients with FVC ≥ 70% have PH.

Repeated 5-fold cross-validation (20 repeats; measurement models and prior from the training folds):

| Inputs | AUC, mPAP ≥ 25 (range over repeats) | Error of mPAP, mmHg | Logistic regression, AUC (complete cases) |
|---|---|---|---|
| FVC + DLCO | 0.70 (0.66–0.72) | 9.9 | 0.64 (43) |
| + echocardiographic sPAP | 0.72 (0.69–0.75) | 9.8 | 0.70 (43) |
| + BNP | 0.72 (0.68–0.76) | 9.9 | 0.68 (42) |
| + pulmonary artery diameter | 0.82 (0.78–0.84) | 8.9 | 0.86 (42) |
| + pulmonary artery diameter + BNP | 0.83 (0.80–0.85) | 8.9 | 0.82 (41) |
| + all optional inputs | 0.86 (0.83–0.87) | 8.6 | 0.90 (34) |

The 90% interval covers 72–77% of the patients. With the calibrated cohort or the mixed-lesion cohort as reference the AUC changes by 0.01 or less; with the field cohort it is 0.82 with all inputs. Measurement models: pulmonary artery diameter 0.25 mm per mmHg (SD 4.0 mm), echocardiographic sPAP 0.42 mmHg per mmHg (SD 16.5 mmHg).

## 8. Limitations

- The lesion is a resistance added to the arteriole, not a cellular mechanism; its five parameters were chosen against the mild-to-moderate benchmarks. Advanced disease, the pressure–flow response and the patients are out of sample.
- Pre-capillary pulmonary hypertension only: the left atrial pressure is fixed.
- Alveolar ventilation is fixed (4.2 L/min; 16 breaths/min, 0.41 L); a quarter of the calibrated cases are hypercapnic, which lowers their PaO₂ and alveolar–arterial gradient (section 7.1). Blood gases at altitude have no ventilatory response or acclimatisation.
- The vessel ageing offset (`age_M = 0.3`) has no measured anchor; it adds 2.2 mmHg to every case.
- Conduit arteries inside scar remodel freely to their shear target (`remodel_conduits = True`).
- Fibrotic and honeycomb units keep 15% of normal compliance; FVC is a linear approximation of the pressure–volume curve. The DLCO–FVC relation of the model is tighter than in patients.
- The arterial morphometry is from one lung; b has no interval. The pulsatility of the pulmonary artery is an effective value.
- The individual estimate was evaluated in one cohort, referred for suspected PH (12 patients without PH) and mostly connective-tissue disease; no thresholds are proposed.

## References

Feaver RE, Gelfand BD, Blackman BR (2013) Nat Commun 4:1525. · Huang W, Yen RT, McLaurine M, Bledsoe G (1996) J Appl Physiol 81:2123–2133. · Pérez-Calixto D, et al. (2026) bioRxiv, doi:10.64898/2026.09.29.755319. · Raghu G, et al. (2015) Eur Respir J 46:1370–1377. · Shorr AF, et al. (2007) Eur Respir J 30:715–721. · Singhal S, et al. (1973) Circ Res 33:190–197. · Womersley JR (1955) J Physiol 127:553–563. · Zisman DA, et al. (2007) Respir Med 101:2153–2159.
