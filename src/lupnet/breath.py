"""Fast scale: one breathing cycle in the time domain (seconds).

The slow scale (years) never integrates seconds. Instead, at every slow step the breath is solved as a
periodic steady state and the slow processes receive cycle averages: tidal volume and strain of each
unit, fraction of the cycle it is open, and the number of reopening events per breath. Damage and
remodelling rates are expressed per unit time and multiplied by the breathing rate where they are
per-event (atelectrauma). This is the standard two-timescale (averaging) reduction: the fast variable
is slaved to the slow state, the slow variable sees only its cycle average.

Unit model (quasi-static tissue, dynamic recruitment, after Bates & Irvin 2002):
  V_i(t) = V_col,i + x_i(t) (V0_i - V_col,i + C_i (P_L,i(t) - P_exp,i))
  dx/dt = s_o (P_L - P_o)  if P_L > P_o ;  s_c (P_L - P_c)  if P_L < P_c ;  0 otherwise ;  x in [0, 1]
  P_L,i(t) = P_exp,i + u(t)   (common pleural swing; P_exp,i falls with height, pleural gradient)
u(t) is found at each instant so that the total volume follows V_T (1 - cos wt)/2 (volume-targeted
breath, the tidal volume and frequency chosen by minimum mechanical power). Airway resistance is
neglected inside the cycle (unit RC ~0.1 s against a ~4 s breath). x carries over between slow steps
(recruitment memory / hysteresis).
"""
import numpy as np


def breath_cycle(L, VT, f, n_breaths=3, nt=120):
    p = L.p
    N = L.N
    h = L.state == 0
    C, V0, _, _ = L.mech_tissue()
    Pc, Po = L.critical_pressures()
    Pexp = p.PL_FRC + p.PEEP + p.pl_grad * L.z_mech
    Vcol = V0 * p.col_V
    x = L.x.copy()
    x[~h] = 1.0                                    # scar and honeycomb: no recruitment dynamics
    T = 60.0 / f
    dt = T / nt
    t_all = np.arange(n_breaths * nt) * dt
    Vbase = None
    lastV, lastx = [], []
    ev = np.zeros(N)
    for k, t in enumerate(t_all):
        target = VT * 0.5 * (1 - np.cos(2 * np.pi * t / T))
        lo, hi = -20.0, 40.0
        for _ in range(40):
            u = 0.5 * (lo + hi)
            V = Vcol + x * (V0 - Vcol + C * u)
            if Vbase is None:
                break
            if V.sum() - Vbase < target:
                lo = u
            else:
                hi = u
        if Vbase is None:
            u = 0.0
            Vbase = (Vcol + x * (V0 - Vcol)).sum()
            V = Vcol + x * (V0 - Vcol)
        PL = Pexp + u
        dx = np.where(PL > Po, p.s_open * (PL - Po), np.where(PL < Pc, p.s_close * (PL - Pc), 0.0))
        xn = np.clip(x + dx * dt, 0.0, 1.0)
        xn[~h] = 1.0
        if k >= (n_breaths - 1) * nt:
            ev += (x < 0.5) & (xn >= 0.5)
            lastV.append(V); lastx.append(xn)
        x = xn
    lastV = np.array(lastV); lastx = np.array(lastx)
    VTi = lastV.max(0) - lastV.min(0)
    rd = np.where(lastx.max(0) < 0.5, 2, np.where(lastx.min(0) > 0.5, 0, 1))
    rd[~h] = 0
    return dict(x=x, VTi=VTi, open_frac=lastx.mean(0), events=ev, rd=rd,
                Vmean=lastV.mean(0))
