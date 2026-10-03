# LuPNet: theory

This document describes the model and gives every number reported in the README. Parameter names refer to `src/lupnet/params.py`; calibrated values for aged IPF are in `data/calibrated_ipf_params.json`.

## 1. Structure

The lung is a set of N gas-exchange units (1,024 by default; `G = 10` generations of a symmetric tree, or an anatomical 3D tree) perfused by an arterial tree and drained by a venous tree. Each unit lumps a region of acini and its capillary sheet. A unit is healthy, fibrotic or honeycombed; fibrotic and honeycomb units have lower compliance, volume, capillary bed and membrane conductance.

**Breathing and mechanics.** Ventilation is distributed by unit compliance; healthy units can lose surfactant function, collapse and recruit cyclically. The forced vital capacity is computed from the volume–pressure behaviour of all units.

**Gas exchange.** Each unit has a diffusing capacity composed in series of a membrane conductance and a capillary component (Roughton–Forster), DLCO and DLNO follow from the sum over units, and arterial oxygen saturation from ventilation–perfusion matching.

**Perfusion.** Arterial and venous conduits follow Poiseuille's law with viscosity μ; the capillary sheet follows sheet-flow mechanics with gravity, alveolar and vascular pressures; hypoxic pulmonary vasoconstriction acts on the terminal arterioles.

**Disease.** Fibrosis spreads between units with a hazard that depends on injury, neighbourhood and age; time is in model units (6.71 years per unit in the calibrated set, `years_per_unit`). The right ventricle adapts to afterload; exacerbations flood units transiently.

## 2. The shared target law for the pulmonary vessels

### 2.1 Target

For a tree that minimizes viscous dissipation at a fixed maintenance cost, with the cost of a vessel proportional to its mass to the power b, the shear stress of every vessel at the optimum is

$$\tau^*=\tau_0\,r^{\,b-1}\,\ell^{\,(b-1)/2}$$

where r is the lumen radius and ℓ the segment length. b = 1 is Murray's law (uniform shear). The local rule

$$\frac{d\ln r}{dt}=\kappa\left(\frac{S}{\tau^*}-1\right)$$

lowers the total cost for steady flow. With `law="shared"`, one τ₀ is used for every vessel; it is fixed so that the main pulmonary artery has the radius `r_pa`.

### 2.2 Sensed shear

$$S=\bar\tau\left[1+w_1\,\varphi_1\,G(\alpha)\right],\qquad \bar\tau=\frac{4\mu Q}{\pi r^3},\qquad \alpha=r\sqrt{\omega\rho/\mu},\ \ \omega=2\pi\,\mathrm{HR}/60$$

G is the exact Womersley ratio of oscillatory to quasi-steady wall shear in a rigid tube (1 for α → 0, ≈ α/4 for large α), w₁ = 2.865 is the weight of the first harmonic fitted to endothelial NF-κB responses (Feaver et al. 2013), and φ₁ the first-harmonic amplitude of flow relative to its mean (effective values 0.8 in arteries, 0.4 in veins).

### 2.3 The terminal element

The explicit tree stops at generation G. The terminal arteriole of each unit stands for the sub-tree below it; its resistance is that of the same law continued generation by generation down to 10 µm radius.

### 2.4 Tests against human morphometry

Huang et al. (1996) give, for 16 orders of human pulmonary arteries, the number of elements, diameter and length. Each order carries Q/N of a 5 L/min cardiac output.

| Model | Best b | Diameter error (16 orders) |
|---|---|---|
| Pulsatile sensing | 0.90–0.925 | factor 1.14 |
| Mean shear only | 0.825 | factor 1.14 |
| Murray (b = 1, mean only) | – | factor 1.26 |

With b = 0.675, the exponent fitted to systemic arteries across mammals, the predicted main pulmonary artery is twice its real size. The arterial law applied unchanged to the 15 venous orders predicts their diameters within a factor 1.14–1.18 (venous φ₁ 0.4–0.8). Integrated in the lung model (`law="shared"`, b = 0.925), the explicit arterial generations match Huang within 6%, and the healthy mean pulmonary artery pressure is 14.8 mmHg (calibrated value of the previous, per-vessel law: 14.0).

## 3. Structural vasculopathy

A non-remodellable lesion in series with each terminal arteriole, independent of fibrosis:

$$R_\mathrm{les}=R_{ta,0}\,A\left(1-e^{-(t_\mathrm{yr}+t_0)/\tau}\right),\qquad c_\mathrm{cap}=\frac{1}{1+f_\mathrm{cap}\,A\left(1-e^{-(t_\mathrm{yr}+t_0)/\tau}\right)}$$

A is the maximum severity of the patient (`vasc_A`), t₀ the years of vasculopathy before the start (`t_vasc0`), τ its time constant (`vasc_tau`), and c_cap multiplies the capillary bed of every unit (`cap_frac` = f_cap). With A = 0 the model is unchanged.

**Simulated cases** (`src/ph_cohort.py`): with probability 0.5 a patient has vasculopathy, A ~ lognormal(median 1.5, σ = 1.1), t₀ ~ U(0, 10) years, τ = 2 years, f_cap = 0.05; the speed of fibrosis varies between patients (years per model unit × lognormal, σ = 0.5); one visit at a random stage.

## 4. Pulmonary hypertension in IPF

**Without vasculopathy** (60 patients), mPAP in mild-to-moderate disease is 18.2 ± 1.1 mmHg, r(mPAP, FVC) = −0.88 and r(mPAP, DLCO) = −0.93, and no patient reaches 25 mmHg. Varying hypoxic vasoconstriction or pressure-driven maladaptive remodelling between patients does not change this.

**With vasculopathy** (200 patients, 162 with FVC 50–90%): see the README table. All seven targets (ARTEMIS-IPF, Raghu et al. 2015; Zisman et al. 2007; the DLCO ratio of a Japanese cohort) lie inside the 95% bootstrap intervals; correlations include measurement noise (mPAP 2.5 mmHg, FVC and DLCO 5 points).

**Not fitted:**

| | Model | Reference |
|---|---|---|
| PH with FVC < 50% (calibrated cohort, n = 38) | 47% | 46% (Shorr et al. 2007, 2,525 transplant candidates) |
| mPAP > 40 mmHg, dedicated advanced cohort (FVC 34–66%, n = 60) | 10% (3–18) | ≈ 9% |
| FVC, PH vs no PH, advanced cohort | 47 vs 48% | 48.4 vs 51.4% |
| PH, dedicated advanced cohort | 67% (55–78) | 46% (overestimated; see limitations) |
| AUC of FVC/DLCO for mPAP ≥ 25 | 0.78 with noise | 0.69–0.74 in ILD cohorts |
| AUC of DLCO alone / FVC alone | 0.75 / 0.65 | ≈ 0.80 / no discrimination |

**Vascular phenotype.** The deficit of DLCO relative to the value expected for the FVC is bimodal (two-component mixture preferred by ΔBIC = 17.6). The low-DLCO component has 27 of 200 patients, of whom 85% carry vasculopathy and 89% have PH; in the rest, 43% carry vasculopathy and 13% have PH. Most vasculopathy is therefore subclinical in lung function, and only its severe forms show as a disproportionately low DLCO.

**Robustness.** The same 12 patients with 256 and 1,024 units differ by +0.06 ± 0.41 mmHg in mPAP (r = 0.997), −0.03 ± 2.2 points in FVC and −0.3 ± 0.7 points in DLCO.

## 5. Limitations

- The simulated FVC flattens near 42% in advanced disease, which inflates the tissue component at low FVC.
- Five vasculopathy parameters were chosen against the mild-to-moderate targets; only the advanced-disease results and the phenotype analysis are out of sample.
- The arterial morphometry is from one lung; b has no interval yet.
- The pulmonary-artery pulsatility is an effective value; measured waveforms should replace it.

## References

Feaver RE, Gelfand BD, Blackman BR (2013) Nat Commun 4:1525. · Huang W, Yen RT, McLaurine M, Bledsoe G (1996) J Appl Physiol 81:2123–2133. · Raghu G, et al. (2015) Eur Respir J 46:1370–1377 (ARTEMIS-IPF haemodynamics). · Shorr AF, et al. (2007) Eur Respir J 30:715–721. · Zisman DA, et al. (2007) Respir Med 101:2153–2159. · Womersley JR (1955) J Physiol 127:553–563.
