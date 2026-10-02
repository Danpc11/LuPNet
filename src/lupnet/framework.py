"""The model as one variational core plus a stochastic biological layer, in three time levels.

    F[Q, V, f ; r ; s] = Sum Q^2 R(r) + lambda Sum (pi r^2 L)^b          vascular network (set-point theory)
                       + W_breath(V, f)                                   breathing power
                       + Sum_i [ E_elastic,i + gamma(s_i) A_i ]           parenchyma (tissue + surfactant)

    fast   (seconds): Q, P, V, f = argmin F given r, s     Kirchhoff (Thomson), min power, mechanics, Fick/Roughton-Forster
    medium (weeks):   dr/dt = -Gamma dF/dr (+ aging, maladaptive: non-variational terms)
                      = the local shear set-point rule d ln r/dt = kappa (tau/tau_set - 1), tau_set = tau0 r^(b-1) L^((b-1)/2)
    slow   (years):   parenchymal state s (surfactant quality, recruitment, scar) jumps with rates set by cycle
                      averages of the fast level; acute exacerbations are Poisson shocks on s
    readout:          z_P = z_F z_R ,  z_dP = z_V z_E        (LiPNet loads)

Check C19 verifies that the medium level is a descent of the Murray/Hu-Cai functional when the
non-variational terms are off, and that an aging offset climbs away from its minimum.
These functions only name the levels; the implementation lives in lung.py, breath.py and sim.py.
"""
import numpy as np
from .lung import Lung
from . import sim


def fast(L: Lung, **kw):
    """Fast level: flows, pressures, breath, gas exchange at fixed vessel structure and parenchymal state."""
    return L.equilibrate(remodel=False, **kw)


def medium(L: Lung, weeks: float = None):
    """Medium level: set-point remodelling for `weeks` of physical time (None = to its rest point)."""
    L._rem_budget = None if weeks is None or L.p.tau_rem_weeks <= 0 else weeks / L.p.tau_rem_weeks
    st = L.equilibrate(remodel=True)
    L._rem_budget = None
    return st


def slow(L: Lung, st, dt: float, rng):
    """Slow level: one stochastic step of the parenchymal state (model time units)."""
    lf, lh = sim.hazards(L, st)
    L.state[rng.random(L.N) < 1 - np.exp(-lf * dt)] = 1
    L.state[(rng.random(L.N) < 1 - np.exp(-lh * dt)) & (L.state == 1)] = 2
    if L.p.parenchyma == "surfactant":
        sim.surfactant_step(L, st, dt)


def functional(L: Lung, st):
    """Value of the Murray/Hu-Cai functional (W) for the current vascular state."""
    return L.functional_E(st["vs"])
