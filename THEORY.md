# LuPNet: theory

This document describes the model and gives every number reported in the README. Parameter names refer to `src/lupnet/params.py`; calibrated values for aged IPF are in `data/calibrated_ipf_params.json`.

## 1. Structure

The lung is a set of N gas-exchange units (1,024 by default, `G = 10` generations of a symmetric tree; the in silico cohorts of the manuscript use 256 units, `G = 8`, which give the same results within 0.1 mmHg on average; or an anatomical 3D tree) perfused by an arterial tree and drained by a venous tree. Each unit lumps a region of acini and its capillary sheet. A unit is healthy, fibrotic or honeycombed; fibrotic and honeycomb units have lower compliance, volume, capillary bed and membrane conductance.

**Breathing and mechanics.** Ventilation is distributed by unit compliance; healthy units can lose surfactant function, collapse and recruit cyclically. The forced vital capacity is computed from the volume–pressure behaviour of all units.

**Gas exchange.** Each unit has a membrane conductance and a capillary volume; its diffusing capacity is their product (`rf = False`, the setting of every cohort of the manuscript; `rf = True` gives the Roughton–Forster series form, which the cohorts do not use), DLCO is the sum over ventilated units, and arterial oxygen saturation follows from ventilation–perfusion matching. Alveolar ventilation is held at 4.2 L/min; PaCO₂ is a readout (0.863 VCO₂ / effective VA), not a regulated variable, so cases with large alveolar dead space are hypercapnic (49 of 200 calibrated cases above 45 mmHg, 14 above 50).

**Perfusion.** Arterial and venous conduits follow Poiseuille's law with viscosity μ; the capillary sheet follows sheet-flow mechanics with gravity, alveolar and vascular pressures; hypoxic pulmonary vasoconstriction acts on the terminal arterioles.

**Disease.** Fibrosis spreads between units with a hazard that depends on surfactant loss, collapse and neighbourhood (`beta = 0` in the calibrated set: no strain feedback); the vessels age after the healthy baseline (`age_M = 0.3`, `age_sigma = 0.3`: rest point of the remodelling rule at 1.3 times the target shear), which by itself raises mPAP from 14.3 to 16.5 mmHg in the 256-unit lung with no fibrosis (table S19); time is in model units (6.71 years per unit in the calibrated set, `years_per_unit`). The right ventricle adapts to afterload; exacerbations flood units transiently.

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

With b = 0.675, the exponent fitted to systemic arteries across mammals, the predicted main pulmonary artery is twice its real size. The arterial law (b = 0.90) applied unchanged to the 15 venous orders predicts their diameters within a factor 1.16 (venous φ₁ = 0.4). In these morphometric tests τ₀ is a free scale fitted by least squares in log D over the 16 orders, with the apparent viscosity of each order. Integrated in the lung model (`law="shared"`, b = 0.925, inside the flat region of the error profile: factor 1.136 vs 1.135 at b = 0.90) τ₀ is instead fixed by the radius of the main pulmonary artery and μ = 3.5 mPa s; the explicit arterial generations then match Huang with a geometric error factor of 1.06 (largest deviation 8.2%), and the healthy mean pulmonary artery pressure is 14.8 mmHg (14.6 with the 256-unit tree of the cohorts; calibrated value of the previous, per-vessel law: 14.0). The law sets the arterial conduits, the venous conduits and the terminal sub-tree, 47% of the healthy pressure drop (0.32 + 0.67 + 2.13 mmHg of 6.6); the capillary sheet (2.10) and the venules (1.39) keep the resistances of the reference calibration, `dP0` × `split` (table S18). The healthy mPAP is therefore an output of the law for its arterial half and an input for the rest.

## 3. Structural vasculopathy

A non-remodellable lesion in series with each terminal arteriole, independent of fibrosis:

$$R_\mathrm{les}=R_{ta,0}\,A\left(1-e^{-(t_\mathrm{yr}+t_0)/\tau}\right),\qquad c_\mathrm{cap}=\frac{1}{1+f_\mathrm{cap}\,A\left(1-e^{-(t_\mathrm{yr}+t_0)/\tau}\right)}$$

A is the maximum severity of the patient (`vasc_A`), t₀ the years of vasculopathy before the start (`t_vasc0`), τ its time constant (`vasc_tau`), and c_cap multiplies the capillary bed of every unit (`cap_frac` = f_cap). With A = 0 the model is unchanged.

**Simulated cases** (`src/ph_cohort.py`): with probability 0.7 a simulated case has vasculopathy, A ~ lognormal(median 1.5, σ = 1.1), t₀ ~ U(0, 10) years, τ = 2 years, f_cap = 0.05; the speed of fibrosis varies between simulated cases (years per model unit × lognormal, σ = 0.5); one visit at a random stage. With these values the lesion has reached its final severity A at the visit in almost every case (saturation 1 − e^{−(t+t₀)/τ}: minimum 0.60, 5th percentile 0.96, median 1.00 in the calibrated cohort; minimum 0.96 in advanced disease; table S20). The cohort results therefore identify the prevalence and the distribution of A, not τ or t₀, and the independence of the lesion from fibrosis is an assumption of the design, tested against the alternative `vasc_mode="fibrosis"`, not a quantity inferred from the cohorts.

## 4. Pulmonary hypertension in IPF

**Without vasculopathy** (60 simulated cases), mPAP in mild-to-moderate disease is 17.6 ± 0.5 mmHg (of which 16.5 is the aged lung before any fibrosis, table S19), r(mPAP, FVC) = −0.92 and r(mPAP, DLCO) = −0.92, and no case reaches 25 mmHg. Varying hypoxic vasoconstriction or pressure-driven maladaptive remodelling between cases spreads mPAP by only 0.8 mmHg (maximum 20.3 mmHg).

**With vasculopathy** (200 simulated cases, 113 with FVC 50–90%): see the README table. All seven targets (ARTEMIS-IPF, Raghu et al. 2015; Zisman et al. 2007; the DLCO ratio of a Japanese cohort) lie inside the 95% bootstrap intervals; correlations include measurement noise (mPAP 2.5 mmHg, FVC and DLCO 5 points).

**Not fitted:**

| | Model | Reference |
|---|---|---|
| PH with FVC < 50%, pooled calibrated and advanced cohorts (n = 147) | 48% (39–56) | 46% (Shorr et al. 2007, 2,525 transplant candidates) |
| mPAP > 40 mmHg, FVC < 50%, pooled (n = 147) | 8% (4–13) | ≈ 9% |
| FVC, PH vs no PH, FVC < 50% | 35.7 vs 39.1% | 48.4 vs 51.4% |
| AUC of FVC/DLCO for mPAP ≥ 25 (mean of 1,000 noise realisations) | 0.75 | 0.69–0.74 in ILD cohorts |
| AUC of DLCO alone / FVC alone | 0.74 / 0.65 | ≈ 0.80 / no discrimination |

**Vascular phenotype.** The deficit of DLCO relative to the value expected for the FVC is bimodal (two-component mixture preferred by ΔBIC = 78.8). The low-DLCO component has 50 of 200 simulated cases, all with vasculopathy, and 84% have PH; in the rest, 60% carry vasculopathy and 11% have PH. Most vasculopathy is therefore subclinical in lung function, and only its severe forms show as a disproportionately low DLCO. The separation rests on the near-deterministic DLCO–FVC relation of the model (r = 0.97 at the visit, against about 0.5–0.7 in patients): with the inter-individual spread of the healthy baseline of table S8 plus measurement noise, ΔBIC falls to a median of 7 (range −12 to 46 over 20 replicates), a second component is preferred in 40% of replicates, and the sensitivity for PH falls from 71% to 34% (table S21). The phenotype is a prediction to be tested in patients, not a result that the simulated cohorts establish.

**Robustness.** The same 12 simulated cases with 1,024 and 256 units differ by −0.08 ± 0.41 mmHg in mPAP (r = 0.997), −0.6 ± 1.7 points in FVC and −0.2 ± 0.8 points in DLCO.

## 5. Histological pattern, altitude and where the lesion sits

**Pattern** (`ph_cohort.py run --pattern`). The calibrated set is patchy without subpleural predominance (units lost: 0.46 subpleural, 0.44 core; honeycomb 22% of lost units). `uip` adds a subpleural seed weighted to the base (`sub_gain = 3`: surfactant decline and direct fibrosis hazard × (1 + sub_gain (1.5 − y)) in the two outer rings, y = 0 base), stronger spread by contiguity (`kn = 1.5`, `kn_s = 3`) and honeycombing (`lam_hc = 0.5`): units lost 0.60 subpleural against 0.25 core, 0.61 base against 0.41 apex, honeycomb 13%, 3.9 years to the visit. `nsip` removes the seed and the contiguity, makes vulnerability nearly uniform (`at2_sigma = 0.3`), honeycombing rare (`lam_hc = 0.05`) and the scar softer (`fib_C = 0.30`): 0.47 against 0.38, 0.68 against 0.30, honeycomb 2%, 6.9 years. PH with FVC 50–90% is 17% in both (calibrated 19%); UIP meets the seven benchmarks and NSIP six (DLCO ratio 0.85, a benchmark measured in IPF). In advanced disease PH is 45% with UIP and 35% with NSIP (transplant candidates with IPF: 46%).

**Pressure surface.** In every cohort mPAP is a function of two variables, the fraction of lung lost and the effective lesion: mPAP = F(lost, A), R² 0.996–0.999. Without lesion, F rises from 16 mmHg (the aged lung) to 20 mmHg at 60% of lung lost; a lesion A = 2 adds 4 mmHg at 5% lost and 6 mmHg at 60% lost. The effect of fibrosis is therefore a partial dependence: small alone, and an amplifier of the lesion, because the same cardiac output crosses fewer units. The surfaces of UIP and NSIP differ by 1 mmHg or less.

**Altitude** (`--PB 585`). Hypoxic vasoconstriction is normalised to the alveolar PO₂ of each lung at baseline, so a healthy lung has the same mPAP at any barometric pressure; in disease, pressures at 585 mmHg differ by less than 1 mmHg from sea level (PH 16% against 19%). Blood gases do change (SaO₂ 85% against 94% at the visit), and without a ventilatory response or acclimatisation they are lower than in residents; SpO₂ is therefore not used as an input of the tool.

**Where the lesion sits** (`--where`, `sim.lesion_weights`). The series lesion of unit i is R_ta0 · A · w_i. `diffuse`: w = 1 (manuscript). `scar`: w = 1 in fibrotic and honeycomb units, 0 in healthy ones. `field`: w = 1 in the scar and, in a healthy unit, the fraction of its neighbours that are scar. `mixed`: w = m + (1 − m) · field. `--a_scale` multiplies A so that every mechanism has the PH prevalence of the diffuse cohort in mild-to-moderate disease. Results (`src/lesion_location.py`, tables S25–S26): the scar-only lesion gives no PH in 200 cases (maximum 23.9 mmHg); the field lesion meets 4 of 7 benchmarks, with r(mPAP, FVC) = −0.34, PH 68% and mPAP > 40 mmHg 31% with FVC < 50%, and no PH with FVC ≥ 70%; half diffuse and half field meets the seven benchmarks but already gives 16% (11–22) of mPAP > 40 mmHg with FVC < 50% (patients about 9%). A lesion inside the scar does not raise pressure because the capillary bed of those units is largely gone and flow goes elsewhere; the lesion that raises pressure is the one in preserved lung, and at most a minor share of it can follow the scar.

## 6. Non-invasive individual estimate

`src/noninvasive_tool.py`. Posterior over the simulated visits: p(visit | data) ∝ prior × p(FVC, DLCO | visit) × Π p(x_k | mPAP). FVC and DLCO are required and carry measurement noise (5 points) plus the spread of the healthy baseline (13% and 15%). Optional inputs x_k (echocardiographic sPAP, pulmonary artery diameter, log BNP, RV fractional area change, right atrial area) are linear in mPAP with Gaussian noise, fitted on catheterised patients in the training folds. The prior gives the simulated visits the distribution of mPAP of the population in five classes (`--prior severity`). Output: mPAP with a 90% interval, probability of PH, PVR, and the split of the pressure into aged lung, lung lost and lesion from the surface F. In 43 catheterised patients with fibrosing ILD (22 NSIP, 9 UIP, 12 other patterns; repeated 5-fold cross-validation): AUC 0.71 with FVC and DLCO, 0.81 with the pulmonary artery diameter, 0.87 with every optional input; error of mPAP 8.5 mmHg; the 90% interval covers 74%. The reference can be the cohort of the patient's pattern or the calibrated cohort, diffuse or half field, with the same result within 0.01 of AUC; a field-only reference lowers it to 0.82. Measured slopes: pulmonary artery diameter 0.25 mm per mmHg (SD 4.0 mm), echocardiographic sPAP 0.42 mmHg per mmHg (SD 16.5), so the estimate of the pressure itself is limited by the measurements, not by the model.

**Pressure–flow response** (`src/exercise_response.py`, Figure 6, table S16). The first 70 cases of the calibrated cohort (47 with lesion, 23 without) are re-equilibrated at imposed cardiac outputs of 7.5 and 10 L/min with the vessels held at the calibre they had reached; ventilation and oxygen uptake scale with the flow. The slope is (mPAP at 10 L/min − resting mPAP) / (10 − resting output). The resting right ventricle of the model delivers about 7 L/min when 10 are demanded, so the flow has to be imposed; versions of the manuscript up to v21 reported the ventricle-limited state as 10 L/min and its slope is half the correct one (docs/AUDIT.md, finding 21).

## 7. Limitations

- The clinical cohort used for the tool was referred for suspected PH (69% with PH, 12 patients without) and is mostly connective-tissue ILD; thresholds need a second, unselected cohort. Post-capillary PH is not represented (LAP is fixed at 8 mmHg).
- Fibrotic and honeycomb units keep 15% of normal compliance; FVC is a linear approximation of the pressure–volume curve.
- Five vasculopathy parameters were chosen against the mild-to-moderate targets; only the advanced-disease results and the phenotype analysis are out of sample.
- The arterial morphometry is from one lung; b has no interval yet.
- The pulmonary-artery pulsatility is an effective value; measured waveforms should replace it.
- Alveolar ventilation is fixed; the breathing pattern (16 breaths/min, 0.41 L) does not change with restriction and a quarter of the calibrated cases are hypercapnic, which lowers their PaO₂ and A–a gradient (Figure 5).
- The vessel ageing offset (`age_M = 0.3`) has no measured anchor; it adds 2.2 mmHg to every case.
- Conduit arteries inside scar remodel freely to their shear target (`remodel_conduits = True`); the opposite setting has not been run as a sensitivity arm.

## References

Feaver RE, Gelfand BD, Blackman BR (2013) Nat Commun 4:1525. · Huang W, Yen RT, McLaurine M, Bledsoe G (1996) J Appl Physiol 81:2123–2133. · Raghu G, et al. (2015) Eur Respir J 46:1370–1377 (ARTEMIS-IPF haemodynamics). · Shorr AF, et al. (2007) Eur Respir J 30:715–721. · Zisman DA, et al. (2007) Respir Med 101:2153–2159. · Womersley JR (1955) J Physiol 127:553–563.


## 8. Non-invasive markers in the simulated cohorts (Figure 5)

The trajectories also give resting SpO2, the arterial blood gases and the systolic pulmonary artery pressure that echocardiography would estimate (mPAP = 0.61 sPAP + 2 mmHg). At a given FVC and DLCO, PaO2 falls with the severity of the vasculopathy, because capillary loss limits oxygen diffusion. With measurement noise, the area under the ROC curve for relevant vasculopathy (effective lesion at least equal to the baseline arteriolar resistance) is 0.69 with FVC and DLCO, 0.72 with SpO2, 0.74 with PaO2, 0.77 with the alveolar-arterial gradient, 0.84 with echocardiographic pressure and 0.85 with the gradient and echocardiographic pressure; for PH, 0.78, 0.84, 0.87, 0.87, 0.92 and 0.94. An emulator of the simulator (`src/noninvasive.py`) allows a posterior over each case's disease state; in held-out simulated cases the error of the estimated mPAP falls from 4.2 to 2.2 mmHg, with calibrated 90% intervals. The calculator uses the same principle with the simulated visits as reference.
