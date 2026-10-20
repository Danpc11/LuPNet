"""Regression test for the numbers published in the manuscript.

The invariant tests check that the model is self-consistent; these check that the
*results* have not moved. Any change to the model, to the cohorts or to the
statistics routine that shifts a published quantity makes this test fail, which
is the point: a number in the paper should never change silently.

Values are those of LuPNet 0.3.3-0.3.4, recomputed from the cohorts shipped in
results/ph_cohorts/. Tolerances are the last digit the manuscript quotes.
"""
import os
import sys

import pandas as pd
import pytest
from scipy.stats import pearsonr

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
COHORTS = os.path.join(ROOT, "results", "ph_cohorts")

from cohort_io import cases  # noqa: E402
import ph_stats  # noqa: E402


def final_visits(tag):
    """One row per simulated case: its last visit, lung function as % of its own healthy lung."""
    rows = []
    for _, d in cases(COHORTS, tag):
        r0, x = d.iloc[0], d.iloc[-1]
        rows.append(dict(mPAP=x.mPAP, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO,
                         A=float(x.vasc_A)))
    R = pd.DataFrame(rows)
    assert len(R), f"cohort {tag} not found in {COHORTS}"
    R["vasc"] = R.A > 0
    R["PH"] = R.mPAP >= 25
    return R


@pytest.fixture(scope="module")
def mild():
    CAL = final_visits("big2")
    assert len(CAL) == 200
    return CAL[(CAL.FVC >= 50) & (CAL.FVC <= 90)].reset_index(drop=True)


def test_mild_cohort_size(mild):
    assert len(mild) == 113


def test_calibration_targets(mild):
    """Table 1: the seven benchmarks of mild-to-moderate IPF."""
    p = ph_stats.point(mild)
    assert p["lt20"] == pytest.approx(54, abs=1)
    assert p["b20_25"] == pytest.approx(27, abs=1)
    assert p["ge25"] == pytest.approx(19, abs=1)
    assert p["ge35"] == pytest.approx(5, abs=1)
    assert p["rDLCO"] == pytest.approx(-0.38, abs=0.02)
    assert p["rFVC"] == pytest.approx(-0.08, abs=0.02)
    assert p["ratio"] == pytest.approx(0.81, abs=0.02)


def test_vasculopathy_prevalence_and_pressure(mild):
    """70% of cases carry a lesion; without one, pressure stays low."""
    assert mild.vasc.mean() == pytest.approx(0.65, abs=0.05)
    no, yes = mild.mPAP[~mild.vasc], mild.mPAP[mild.vasc]
    assert no.mean() == pytest.approx(17.7, abs=0.1)
    assert yes.mean() == pytest.approx(23.7, abs=0.2)
    assert not mild.PH[~mild.vasc].any(), "PH without vasculopathy in mild-to-moderate disease"


def test_fibrosis_alone_cannot_produce_ph():
    """Cohort 1: fibrosis and hypoxia only, FVC 50-90% predicted (n = 35 of 60)."""
    H = final_visits("hom")
    assert len(H) == 60
    M = H[(H.FVC >= 50) & (H.FVC <= 90)]
    assert len(M) == 35
    assert M.mPAP.mean() == pytest.approx(17.6, abs=0.1)
    assert M.mPAP.std() == pytest.approx(0.5, abs=0.1)
    assert M.mPAP.max() < 25
    # The manuscript quotes these two without measurement noise: with a spread of
    # 0.5 mmHg, 2.5 mmHg of noise would erase the correlation. Table 1 instead
    # quotes noise-averaged correlations, because those are compared with patients.
    assert pearsonr(M.mPAP, M.FVC)[0] == pytest.approx(-0.92, abs=0.02)
    assert pearsonr(M.mPAP, M.DLCO)[0] == pytest.approx(-0.92, abs=0.02)


def test_variable_vascular_responses():
    """Cohort 2: hypoxic vasoconstriction and remodelling vary between cases (n = 22 of 38)."""
    H = final_visits("het")
    assert len(H) == 38
    M = H[(H.FVC >= 50) & (H.FVC <= 90)]
    assert len(M) == 22
    assert M.mPAP.max() == pytest.approx(20.3, abs=0.2)
    assert M.mPAP.max() < 25


def test_advanced_disease_out_of_sample():
    """Nothing refitted: PH and severe PH with FVC below 50%, pooled cohorts."""
    POOL = pd.concat([final_visits("big2"), final_visits("adv2")], ignore_index=True)
    A = POOL[POOL.FVC < 50]
    assert len(A) == 147
    assert 100 * A.PH.mean() == pytest.approx(48, abs=1)
    assert 100 * (A.mPAP > 40).mean() == pytest.approx(8, abs=1)


def test_lesion_tied_to_fibrosis_fails():
    """The alternative mechanism: PH becomes rare and pressure follows FVC again."""
    M = final_visits("altF")
    M = M[(M.FVC >= 50) & (M.FVC <= 90)]
    assert 100 * M.PH.mean() == pytest.approx(4, abs=2)
    assert ph_stats.noise_avg_r(M, "FVC") == pytest.approx(-0.27, abs=0.03)
