# Data

No patient data are stored in this repository.

| File | Content | Source |
|---|---|---|
| `huang1996_pulmonary_arteries.tsv` | Human pulmonary arterial tree by diameter-defined Strahler order: number of elements `N`, diameter `D_cm`, length `L_cm`, apparent blood viscosity `mu_cP`. Orders 1 to 15; order 16 is the main pulmonary artery, added here. | Orders 1–15 transcribed from Huang W, Yen RT, McLaurine M, Bledsoe G. Morphometry of the human pulmonary vasculature. J Appl Physiol 1996;81:2123–2133. Order 16: diameter 3.0 cm from Singhal S, et al. Circ Res 1973;33:190–197; length 9.05 cm from a compilation of morphometric data (to be confirmed against the primary source, docs/AUDIT.md, open items). |
| `huang1996_pulmonary_veins.tsv` | The same quantities for the 15 orders of the pulmonary venous tree. | Huang et al. 1996, as above. |
| `calibrated_ipf_params.json` | Parameters of the lung model for aged idiopathic pulmonary fibrosis used by every cohort; names are those of `src/lupnet/params.py`. Fields that are absent take the defaults of that file (for example `rf = False`). The cohort scripts set `law = "shared"`, `b = 0.925` and `G = 8` on top of it. | Calibrated in this work. |
| `virtual_ph_cohort.tsv` | Final visit of every simulated case of the calibrated and advanced cohorts (mPAP, FVC and DLCO in % of the case's own baseline, lesion). | Written by `src/fit_calculator.py` from `results/ph_cohorts/big2.tsv.gz` and `adv2.tsv.gz`. |

The transcribed morphometric tables keep the terms of their publishers; they are reproduced here as numerical data with attribution.
