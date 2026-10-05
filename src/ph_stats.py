"""One statistics routine shared by every script that reports the calibration targets (Table 1, Tables S4-S5,
ph_cohort.py analyze), so that identical quantities get identical point estimates and intervals."""
import numpy as np, pandas as pd
from scipy.stats import pearsonr

TARGETS = dict(lt20=51, b20_25=30, ge25=19, ge35=4, rDLCO=-0.30, rFVC=0.0, ratio=0.73)
SD_MPAP, SD_PFT = 2.5, 5.0          # measurement noise: mPAP (mmHg), FVC and DLCO (percentage points)


def noise_avg_r(D, col, n=1000, seed=0):
    """Pearson r between mPAP and FVC or DLCO, averaged over n realisations of measurement noise."""
    rng = np.random.default_rng(seed)
    return float(np.mean([pearsonr(D.mPAP + rng.normal(0, SD_MPAP, len(D)), D[col] + rng.normal(0, SD_PFT, len(D)))[0] for _ in range(n)]))


def _one(D, rng):
    mp = D.mPAP.values + rng.normal(0, SD_MPAP, len(D)); ph = D.mPAP.values >= 25
    return dict(lt20=100 * (D.mPAP < 20).mean(), b20_25=100 * ((D.mPAP >= 20) & (D.mPAP < 25)).mean(), ge25=100 * ph.mean(),
                ge35=100 * (D.mPAP >= 35).mean(), rDLCO=pearsonr(mp, D.DLCO + rng.normal(0, SD_PFT, len(D)))[0],
                rFVC=pearsonr(mp, D.FVC + rng.normal(0, SD_PFT, len(D)))[0],
                ratio=D.DLCO[ph].mean() / D.DLCO[~ph].mean() if ph.any() and (~ph).any() else np.nan)


def point(D):
    """Point estimates: proportions and ratio without noise; correlations averaged over 1,000 noise realisations."""
    p = _one(D, np.random.default_rng(2)); p["rDLCO"] = noise_avg_r(D, "DLCO"); p["rFVC"] = noise_avg_r(D, "FVC")
    return p


def bootstrap(D, n=1000, seed=1):
    """Percentile 95% intervals from bootstrap resamples of simulated cases (noise redrawn in each resample)."""
    rng = np.random.default_rng(seed); D = D.reset_index(drop=True)
    B = pd.DataFrame([_one(D.iloc[rng.integers(0, len(D), len(D))].reset_index(drop=True), rng) for _ in range(n)])
    return {k: tuple(np.nanpercentile(B[k], [2.5, 97.5])) for k in B}


def prevalence_ci(D, f, n=1000, seed=1):
    """Point estimate and percentile 95% interval of a proportion defined by f(D) -> boolean array."""
    rng = np.random.default_rng(seed); D = D.reset_index(drop=True)
    b = [100 * f(D.iloc[rng.integers(0, len(D), len(D))]).mean() for _ in range(n)]
    return 100 * f(D).mean(), *np.percentile(b, [2.5, 97.5])
