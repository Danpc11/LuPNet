"""Two coupled symmetric trees sharing N = 2^G terminal lung units.

Tree layout (both trees): heap indexing, node k has children 2k+1, 2k+2 and one
incoming edge that carries the same index k; node 0's edge is the trunk (main
pulmonary artery / trachea). Leaves are nodes 2^G-1 ... 2^(G+1)-2.

Leaves are placed on a 2^(G/2) x 2^(G/2) grid by recursive bisection
(odd depths split x, even depths split y), so the tree is an H-tree and grid
neighbours are anatomical neighbours. y = 0 is the base, y = 1 the apex.

Circulation: arterial conduit edge k -> ... -> leaf: terminal arteriole (HPV
acts here) -> capillary sheet -> venule -> venous conduit edge k -> ... -> LA.
Because the venous tree mirrors the arterial one the whole network is
series-parallel and is solved exactly by one bottom-up pass (equivalent
conductance) and one top-down pass (flow split, pressures).
"""
import numpy as np
from .params import Params

MMHG = 133.322                  # Pa
LPM = 1.0e-3 / 60.0             # m^3/s per L/min
CMH2O = 98.0665                 # Pa


def pois(mu, L, r):
    """Poiseuille resistance, Pa s / m^3."""
    return 8.0 * mu * L / (np.pi * r ** 4)



def womersley_G(alpha):
    """Exact |oscillatory / quasi-steady| wall-shear ratio for a rigid tube (linear continuation for alpha > 40)."""
    from scipy.special import jv
    a = np.atleast_1d(np.asarray(alpha, float)); out = np.ones_like(a)
    m = (a > 1e-3) & (a <= 40)
    lam = (1j ** 1.5) * a[m]
    j0, j1 = jv(0, lam), jv(1, lam)
    out[m] = np.abs(lam * j1 / j0 / (1 - 2 * j1 / (lam * j0))) / 4.0
    big = a > 40
    if big.any():
        lam39, lam40 = (1j ** 1.5) * 39.0, (1j ** 1.5) * 40.0
        g = lambda L: abs(L * jv(1, L) / jv(0, L) / (1 - 2 * jv(1, L) / (L * jv(0, L)))) / 4.0
        g39, g40 = g(lam39), g(lam40)
        out[big] = g40 + (a[big] - 40.0) * (g40 - g39)
    return out

class Lung:
    STATES = ("healthy", "fibrotic", "honeycomb")

    def __init__(self, p: Params = None):
        self.p = p = p or Params()
        self._topology()
        self._finish_topology()
        self._grid()
        self.state = np.zeros(self.N, int)
        self._heights()
        self._init_parenchyma()
        self._build_vessels()
        self._build_airways()
        # baseline equilibrium: defines unstressed radii and reference strain
        self.h = np.ones(self.N)
        self._baseline()

    # ------------------------------------------------------------ geometry
    def _grid(self):
        G, N = self.p.G, self.N
        nx_bits = (G + 1) // 2
        ny_bits = G // 2
        ix = np.zeros(N, int); iy = np.zeros(N, int)
        for j, node in enumerate(self.leaf):
            k = node; bits = []
            while k > 0:
                bits.append((self.depth[k], int(k % 2 == 0)))
                k = self.parent[k]
            for d, bit in bits:
                if d % 2 == 1:
                    pos = nx_bits - 1 - (d - 1) // 2
                    ix[j] |= bit << pos
                else:
                    pos = ny_bits - 1 - (d - 2) // 2
                    iy[j] |= bit << pos
        self.ix, self.iy = ix, iy
        self.nx, self.ny = 2 ** nx_bits, 2 ** ny_bits
        self.y = (iy + 0.5) / self.ny                  # 0 base -> 1 apex
        self.x = (ix + 0.5) / self.nx
        cell = -np.ones((self.nx, self.ny), int)
        cell[ix, iy] = np.arange(N)
        self.cell = cell
        nb = []
        for j in range(N):
            a, b = ix[j], iy[j]
            nb.append([cell[a + da, b + db] for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1))
                       if 0 <= a + da < self.nx and 0 <= b + db < self.ny])
        self.nbrs = nb
        edge = np.array([min(ix[j], self.nx - 1 - ix[j], iy[j], self.ny - 1 - iy[j]) for j in range(N)])
        self.pleural_dist = edge                       # 0 = subpleural ring

    # ------------------------------------------------------------ topology (any binary tree)
    def _topology(self):
        """Symmetric binary tree in heap order (subclasses may build any binary tree)."""
        G = self.p.G
        self.N = 2 ** G
        self.n = n = 2 ** (G + 1) - 1
        k = np.arange(n)
        self.parent = np.where(k > 0, (k - 1) // 2, -1)
        self.c1 = np.where(k < 2 ** G - 1, 2 * k + 1, -1)
        self.c2 = np.where(k < 2 ** G - 1, 2 * k + 2, -1)

    def _finish_topology(self):
        """Depth, leaves, levels and leaf counts from parent / child arrays (node 0 = root)."""
        n = self.n
        depth = np.zeros(n, int)
        order = [0]
        for k in order:
            for c in (self.c1[k], self.c2[k]):
                if c >= 0:
                    depth[c] = depth[k] + 1; order.append(c)
        self.depth = depth
        self.leaf = np.where(self.c1 < 0)[0]
        self.leaf_index = -np.ones(n, int); self.leaf_index[self.leaf] = np.arange(len(self.leaf))
        D = int(depth.max())
        self.levels = [np.where(depth == d)[0] for d in range(D + 1)]
        self.inner_levels = [l[self.c1[l] >= 0] for l in self.levels]
        nl = np.zeros(n); nl[self.leaf] = 1
        for l in reversed(self.inner_levels):
            nl[l] = nl[self.c1[l]] + nl[self.c2[l]]
        self.nleaf = nl
        self.Dmax = D

    # ------------------------------------------------------------ vessels
    def _vessel_lengths(self):
        return self.p.L_pa * 2.0 ** (-self.depth / 3.0)

    def _build_vessels(self):
        p, d = self.p, self.depth
        b = p.b
        self.Lv = self._vessel_lengths()
        f_ref = p.Q * LPM * self.nleaf / self.N
        if p.law == "shared":
            return self._build_vessels_shared(f_ref)
        # healthy tree is a rest point of the set-point rule with exponent b:
        # tau = 4 mu f / (pi r^3) = tau0 r^(b-1) L^((b-1)/2)
        tau0 = 4 * p.mu_b * p.Q * LPM / (np.pi * p.r_pa ** (2 + b) * p.L_pa ** ((b - 1) / 2))
        self.r_ref = (4 * p.mu_b * f_ref / (np.pi * tau0 * self.Lv ** ((b - 1) / 2))) ** (1 / (2 + b))
        R_edge = pois(p.mu_b, self.Lv, self.r_ref) * LPM / MMHG      # mmHg/(L/min)
        # dissipation-equivalent resistance of one conduit tree (= sum over levels of R/count if symmetric)
        R_cond = float((R_edge * (self.nleaf / self.N) ** 2).sum())
        R0 = p.dP0 / p.Q
        N = self.N
        R_ta = N * (p.split[0] * R0 - R_cond)
        R_c = N * p.split[1] * R0
        R_vn = N * (p.split[2] * R0 - R_cond)
        assert R_ta > 0 and R_vn > 0
        self.R_vn = R_vn
        self.R_cap0 = R_c
        self.r_ta_ref = (8 * p.mu_b * p.L_term / (np.pi * R_ta * MMHG / LPM)) ** 0.25
        self.R_cond = R_cond
        # venous tree: mirror of the arterial tree in this class (see septal.SeptalLung)
        self.r_ref_v, self.Lvv, self.vdepth = self.r_ref, self.Lv, self.depth

    def _law_radius(self, f, L, tau0, arterial=True):
        """Radius at which the sensed shear equals the shared target: 4 mu f/(pi r^3) F(r) = tau0 r^(b-1) L^((b-1)/2)."""
        from scipy.optimize import brentq
        p, b = self.p, self.p.b
        def g(x):
            r = np.exp(x)
            return np.log(4 * p.mu_b * f / (np.pi * r ** 3) * float(self.pulse_factor(np.array([r]), arterial)[0])) - \
                np.log(tau0 * r ** (b - 1) * L ** ((b - 1) / 2))
        return float(np.exp(brentq(g, np.log(1e-9), np.log(1.0))))

    def _build_vessels_shared(self, f_ref):
        """Shared target law: one tau0, calibrated on the main pulmonary artery radius, sets every conduit and the
        terminal arterioles. Capillary and venule resistances are kept from the reference calibration."""
        p, N = self.p, self.N
        b = p.b
        Fpa = float(self.pulse_factor(np.array([p.r_pa]), True)[0])
        self.tau0_shared = 4 * p.mu_b * p.Q * LPM / (np.pi * p.r_pa ** 3) * Fpa / (p.r_pa ** (b - 1) * p.L_pa ** ((b - 1) / 2))
        cache = {}
        r = np.empty(len(f_ref))
        for k, (f, L) in enumerate(zip(f_ref, self.Lv)):
            key = (round(float(f), 18), round(float(L), 12))
            if key not in cache:
                cache[key] = self._law_radius(float(f), float(L), self.tau0_shared, True)
            r[k] = cache[key]
        self.r_ref = r
        R_edge = pois(p.mu_b, self.Lv, self.r_ref) * LPM / MMHG
        self.R_cond = float((R_edge * (self.nleaf / self.N) ** 2).sum())
        # the lumped terminal element stands for the sub-tree below the last explicit generation: continue the
        # symmetric tree under the same law down to 10 um radius and use its equivalent resistance
        G0 = int(self.depth.max()); fl = p.Q * LPM / N; R_sub = 0.0; g = G0 + 1
        rows = []
        while True:
            f_g = p.Q * LPM / 2 ** g; L_g = p.L_pa * 2 ** (-g / 3)
            r_g = self._law_radius(f_g, L_g, self.tau0_shared, True)
            R_sub += (pois(p.mu_b, L_g, r_g) * LPM / MMHG) / 2 ** (g - G0)
            rows.append((g, r_g, L_g))
            if r_g < 10e-6 or g > G0 + 30:
                break
            g += 1
        self.subtree = rows
        self.R_sub = R_sub
        self.r_ta_ref = (8 * p.mu_b * p.L_term / (np.pi * R_sub * MMHG / LPM)) ** 0.25
        # capillary and venule resistances from the reference (local-law) calibration
        b0 = p.b_ref_calib
        t0 = 4 * p.mu_b * p.Q * LPM / (np.pi * p.r_pa ** (2 + b0) * p.L_pa ** ((b0 - 1) / 2))
        r0 = (4 * p.mu_b * f_ref / (np.pi * t0 * self.Lv ** ((b0 - 1) / 2))) ** (1 / (2 + b0))
        R_cond0 = float(((pois(p.mu_b, self.Lv, r0) * LPM / MMHG) * (self.nleaf / self.N) ** 2).sum())
        R0 = p.dP0 / p.Q
        self.R_cap0 = N * p.split[1] * R0
        self.R_vn = N * (p.split[2] * R0 - R_cond0)
        self.r_ref_v, self.Lvv, self.vdepth = self.r_ref, self.Lv, self.depth

    def _vres(self, r, L):
        return pois(self.p.mu_b, L, r) * LPM / MMHG

    def solve_vessels(self, ru_a, ru_v, ru_ta, c, h, Qtot, P=None, tol=1e-7, maxit=200, Pt=0.0):
        """Pressure-flow state with distensible vessels (r = ru (1 + alpha Ptm)).

        Returns dict with flows per edge Q, leaf flows, pressures, actual radii.
        Distension makes it a fixed point; it is iterated with damping.
        """
        p, N = self.p, self.N
        lv, leaf, par = self.levels, self.leaf, self.parent
        if P is None:
            P = self._P_base
        for it in range(maxit):
            # Pt: extra transmural pressure from tethering by stiff tissue
            r_a = ru_a * (1 + p.alpha * (self._ptm(P["Pm_a"], self.za_edge) + Pt))
            r_v = ru_v * (1 + p.alpha * (self._ptm(P["Pm_v"], self.zv_edge) + Pt))
            r_ta = ru_ta * (1 + p.alpha * self._ptm(P["Pm_ta"], self.z_leaf))
            Ra = self._vres(r_a, self.Lv); Rv = self._vres(r_v, self.Lv)
            Rta = self._vres(r_ta, p.L_term) * h
            Gcap = self._gcap(c, P)
            with np.errstate(divide="ignore"):
                Rterm = Rta + self.R_lesion + np.where(Gcap > 0, 1.0 / np.maximum(Gcap, 1e-300), np.inf) + self.R_vn
            Req = np.empty(self.n); Gs = np.zeros(self.n)
            Req[leaf] = Ra[leaf] + Rv[leaf] + Rterm
            for k in reversed(self.inner_levels):
                g1 = 1 / Req[self.c1[k]]; g2 = 1 / Req[self.c2[k]]
                Gs[k] = g1 + g2
                with np.errstate(divide="ignore"):
                    Req[k] = Ra[k] + Rv[k] + np.where(Gs[k] > 0, 1 / np.maximum(Gs[k], 1e-300), np.inf)
            Q = np.zeros(self.n); Q[0] = Qtot
            for k in self.inner_levels:
                for c_ in (self.c1[k], self.c2[k]):
                    Q[c_] = np.where(Gs[k] > 0, Q[k] * (1 / Req[c_]) / np.maximum(Gs[k], 1e-300), 0.0)
            PPA = p.LAP + Qtot * Req[0]
            Pa = np.empty(self.n); Pv = np.empty(self.n)
            Pa[0] = PPA - Q[0] * Ra[0]; Pv[0] = p.LAP + Q[0] * Rv[0]
            for k in self.levels[1:]:
                Pa[k] = Pa[par[k]] - Q[k] * Ra[k]
                Pv[k] = Pv[par[k]] + Q[k] * Rv[k]
            Pa_up = np.where(par >= 0, Pa[np.maximum(par, 0)], PPA)
            Pv_dn = np.where(par >= 0, Pv[np.maximum(par, 0)], p.LAP)
            Ql = Q[leaf]
            Pta_end = Pa[leaf] - Ql * (Rta + self.R_lesion)
            Pvs = Pv[leaf] + Ql * self.R_vn
            with np.errstate(divide="ignore", invalid="ignore"):
                Pcap_n = np.where(Gcap > 0, Pta_end - 0.5 * Ql / np.maximum(Gcap, 1e-300), Pa[leaf])
            new = dict(Pm_a=0.5 * (Pa_up + Pa), Pm_v=0.5 * (Pv + Pv_dn),
                       Pm_ta=0.5 * (Pa[leaf] + Pta_end), Pcap=Pcap_n, Pte=Pta_end, Pvs=Pvs)
            err = max(np.abs(new[key] - P[key]).max() for key in new)
            w = 0.6
            P = {key: w * new[key] + (1 - w) * P[key] for key in new}
            if err < tol:
                break
        Pcap = P["Pcap"]
        Vcap = self._vcap(c, P)
        return dict(Q=Q, Qv=Q, Ql=Ql, PPA=PPA, dP=PPA - p.LAP, Req=Req[0], r_a=r_a, r_v=r_v,
                    r_ta=r_ta, P=P, Vcap=Vcap, it=it, err=err, Pcap=Pcap)

    # ------------------------------------------------------------ gravity, zones, extravascular pressure
    def _heights(self):
        """Vertical position (cm, + above the heart) of units and vessel midpoints."""
        p = self.p
        self.z_leaf = (self.y - p.heart_frac) * p.lung_height_cm if p.perf_grav else np.zeros(self.N)
        zn = self.subtree_mean(self.z_leaf)
        zpar = np.where(self.parent >= 0, zn[np.maximum(self.parent, 0)], 0.0)
        self.za_edge = 0.5 * (zn + zpar)
        self.zv_edge = self.za_edge

    def _pev_ea(self):
        """Extravascular pressure of extra-alveolar vessels (mmHg): pleural pressure minus an
        interdependence term that grows with transpulmonary pressure (tethering)."""
        p = self.p
        PL = p.P_alv - p.P_pl
        return (p.P_pl - p.kappa_ev * PL) * 0.7356

    def _ptm(self, Phead, z):
        """Transmural pressure of an extra-alveolar vessel. Without perf_grav this is the
        intravascular pressure itself (earlier campaigns). With perf_grav: head minus hydrostatic
        column minus extravascular pressure."""
        p = self.p
        if not p.perf_grav:
            return Phead
        return Phead - p.rho_g * z - self._pev_ea()

    def _sheet(self, P, z):
        """Fung-Sobin sheet: H = 1 + a (P - P_alv) above alveolar pressure, collapsed below.
        Returns H and Phi = int H^3 dP from P_alv."""
        p = self.p
        # extravascular pressure of alveolar vessels: alveolar pressure plus an optional septal-stretch
        # term that grows with transpulmonary pressure above FRC (sigma_sep = 0: pure Fung-Sobin)
        x = P - p.rho_g * z - (p.P_alv + p.sigma_sep * (p.P_alv - p.P_pl - p.PL_FRC)) * 0.7356
        if p.parenchyma == "surfactant":
            # surface tension lowers pericapillary pressure below alveolar pressure (~3.5 mmHg at mid volume),
            # proportionally to the local surface tension
            gam = 1 + p.gamma_gain * (1 - self.sq)
            x = x + p.dst_mid * (gam if np.shape(P) == gam.shape else gam[:np.size(P)])
        a, w = p.sheet_a, p.sheet_w
        # smoothed collapse (softplus of width w mmHg) keeps the zone 2/3 switch differentiable;
        # w -> 0 recovers the sharp Fung-Sobin sheet
        xs = w * np.logaddexp(0.0, x / w)
        Phi = ((1 + a * xs) ** 4 - 1) / (4 * a)
        H = (1 + a * xs) / (1 + np.exp(-np.clip(x / w, -50, 50)))      # open-sheet fraction x height
        return H, Phi, x

    def _gcap(self, c, P):
        p = self.p
        if not p.perf_grav:
            gc = ((1 + p.alpha_c * P["Pcap"]) / (1 + p.alpha_c * self._Pcap0)) ** 3
            return c * gc / self.R_cap0
        Ha, Fa, xa = self._sheet(P["Pte"], self.z_leaf)
        Hv, Fv, xv = self._sheet(P["Pvs"], self.z_leaf)
        d = xa - xv
        xs_a = p.sheet_w * np.logaddexp(0.0, xa / p.sheet_w)
        dPhi = (1 + p.sheet_a * xs_a) ** 3 / (1 + np.exp(-np.clip(xa / p.sheet_w, -50, 50)))
        G = np.where(np.abs(d) > 1e-6, (Fa - Fv) / np.where(np.abs(d) > 1e-6, d, 1.0), dPhi)
        # a collapsed sheet keeps a residual conductance (1e-3 of an open one) so that a fixed
        # inflow can raise upstream pressure until the sheet reopens (and the solve stays regular)
        return c * np.maximum(G, 1e-3) * self.K0

    def _vcap(self, c, P):
        p = self.p
        if not p.perf_grav:
            return c * (p.Vcap / self.N) * (1 + p.alpha_c * P["Pcap"]) / (1 + p.alpha_c * self._Pcap0)
        Ha, _, _ = self._sheet(P["Pte"], self.z_leaf)
        Hv, _, _ = self._sheet(P["Pvs"], self.z_leaf)
        return c * (p.Vcap / self.N) * 0.5 * (Ha + Hv) / self.H_norm

    def zones(self, vs):
        """West zone of each unit: 1 (P_art < P_alv), 2 (P_ven < P_alv <= P_art), 3."""
        _, _, xa = self._sheet(vs["P"]["Pte"], self.z_leaf)
        _, _, xv = self._sheet(vs["P"]["Pvs"], self.z_leaf)
        return np.where(xa <= 0, 1, np.where(xv <= 0, 2, 3))

    def _calibrate_west(self, P0):
        """Sheet constant K0 so a mid-height unit matches R_cap0, then one common factor on the
        microvascular resistances so the healthy upright lung has mPAP - LAP = dP0."""
        p, N = self.p, self.N
        z0 = np.zeros(1)
        _, Fa, xa = self._sheet(P0["Pte"][:1], z0); Ha, _, _ = self._sheet(P0["Pte"][:1], z0)
        _, Fv, xv = self._sheet(P0["Pvs"][:1], z0); Hv, _, _ = self._sheet(P0["Pvs"][:1], z0)
        self.K0 = 1.0 / (self.R_cap0 * float(((Fa - Fv) / (xa - xv))[0]))
        self.H_norm = float((0.5 * (Ha + Hv))[0])
        dPc = 2 * self.R_cond * p.Q
        for _ in range(60):
            vs = self.solve_vessels(self.r_ref, self.r_ref_v, np.full(N, self.r_ta_ref), np.ones(N), np.ones(N), p.Q, P=P0)
            s = (p.dP0 - dPc) / max(vs["dP"] - dPc, 1e-9)
            self.R_vn *= s; self.r_ta_ref *= s ** -0.25; self.K0 /= s
            P0 = vs["P"]
            if abs(s - 1) < 1e-10 and vs["err"] < 1e-9:
                break

    def shear_z(self, vs, ru_a=None):
        """z = |tau| / tau_set for the remodelable vessels (conduits + terminal arterioles)."""
        p, b = self.p, self.p.b
        Qa = np.abs(vs["Q"]) * LPM
        Qt = np.abs(vs["Ql"]) * LPM
        za = (4 * p.mu_b * Qa / (np.pi * vs["r_a"] ** 3)) / (self.tau0_a * vs["r_a"] ** (b - 1) * self.Lv ** ((b - 1) / 2))
        Qvv = np.abs(vs["Qv"]) * LPM
        zv = (4 * p.mu_b * Qvv / (np.pi * vs["r_v"] ** 3)) / (self.tau0_v * vs["r_v"] ** (b - 1) * self.Lvv ** ((b - 1) / 2))
        zt = (4 * p.mu_b * Qt / (np.pi * vs["r_ta"] ** 3)) / (self.tau0_t * vs["r_ta"] ** (b - 1) * p.L_term ** ((b - 1) / 2))
        if p.pulsatile or p.law == "shared":
            za = za * self.pulse_factor(vs["r_a"], arterial=True, depth=self.depth)
            zv = zv * self.pulse_factor(vs["r_v"], arterial=False, depth=self.vdepth)
            zt = zt * self.pulse_factor(vs["r_ta"], arterial=True, depth=np.full(self.N, self.p.G + 1))
        return za, zv, zt

    def pulse_factor(self, r, arterial=True, depth=None):
        """Sensed shear / Poiseuille mean shear for pulsatile flow. Womersley number alpha = r sqrt(w rho/mu):
        at small alpha the oscillatory wall shear is Poiseuille-like, at large alpha it grows ~ alpha (the
        Stokes layer is thinner than the radius), so a vessel senses tau_mean [(1 - phi) + phi sqrt(1 +
        (alpha/wom_c)^2)], phi = pulsatile fraction (damped on the venous side)."""
        p = self.p
        w = 2 * np.pi * p.HR / 60.0
        a = r * np.sqrt(w * p.rho_b / p.mu_b)
        if p.law == "shared":
            phi1 = p.phi1_a if arterial else p.phi1_v
            return 1.0 + p.w1 * phi1 * womersley_G(a).reshape(np.shape(r))
        phi = p.pulse_frac if arterial else p.pulse_frac_v
        if p.pulse_damp > 0 and depth is not None:
            phi = phi * np.exp(-np.asarray(depth) / p.pulse_damp)     # pulsatility damped along the tree
        return (1 - phi) + phi * np.sqrt(1 + (a / p.wom_c) ** 2)

    # ------------------------------------------------------------ airways
    def _build_airways(self):
        p, d = self.p, self.depth
        self.r_aw = p.r_tr * 2.0 ** (-p.aw_r_exp * d)
        self.L_aw = p.L_tr * 2.0 ** (-p.aw_L_exp * d)
        self.R_aw = p.k_inertial * pois(p.mu_a, self.L_aw, self.r_aw) * 1e-3 / CMH2O   # cmH2O s/L
        self.V_cond = float((np.pi * self.r_aw ** 2 * self.L_aw).sum() * 1e3)          # L
        assert self.V_cond < p.VD
        g = 1 + p.grav * (0.5 - self.y); g /= g.mean()
        self.C_h = p.C_lung / self.N * g
        v = 1 + p.grav_V * (self.y - 0.5); v /= v.mean()
        self.V_h = p.FRC / self.N * v
        self.R_unit = self.N * (p.R_periph + p.R_tissue)

    def mech(self):
        p, s = self.p, self.state
        Cf = np.select([s == 0, s == 1, s == 2], [1.0, p.fib_C, p.hc_C])
        Vf = np.select([s == 0, s == 1, s == 2], [1.0, p.fib_V, p.hc_V])
        cf = np.select([s == 0, s == 1, s == 2], [1.0, p.fib_c, p.hc_c])
        cf = cf * getattr(self, 'cap_vasc', 1.0)   # microvascular loss from the vasculopathy (default none)
        Df = np.select([s == 0, s == 1, s == 2], [1.0, p.fib_D, p.hc_D])
        if p.nbr_D < 1.0:
            # optional field effect: healthy units next to scar lose membrane conductance
            Df = np.where(s == 0, Df * (1 - (1 - p.nbr_D) * self.stiff_nbr_frac()), Df)
        if p.parenchyma == "surfactant":
            h = s == 0
            rd = self.rd
            Cf = np.where(h & (rd == 2), p.col_C, np.where(h & (rd == 1), Cf * p.cyc_vent, Cf))
            Vf = np.where(h & (rd == 2), p.col_V, Vf)
            # higher surface tension shrinks the capillary bed (volume and conductance)
            cf = np.where(h, cf * (1 - p.surf_Vcap * (1 - self.sq)), cf)
            # the same epithelial state thickens the membrane (lower Dm) in not-yet-scarred units
            Df = np.where(h, Df * (1 - p.mem_gain * (1 - self.sq)), Df)
            # flooded units (exacerbation): unventilated, stiff, perfused, no membrane exchange
            Cf = np.where(self.flood, p.col_C, Cf); Df = np.where(self.flood, 0.0, Df)
            if p.mem_field > 0 or p.vc_field > 0:
                # diffuse field: membrane thickening and capillary rarefaction in not-yet-scarred tissue,
                # growing with the damaged fraction of the lung (saturating, square root)
                F = np.sqrt(float(((s > 0) | (rd == 2)).mean()))
                Df = np.where(h, Df * max(0.0, 1 - p.mem_field * F), Df)
                cf = np.where(h, cf * max(0.0, 1 - p.vc_field * F), cf)
            cf = np.where(h & (rd == 2), cf * p.col_c, cf)
        return self.C_h * Cf, self.V_h * Vf, cf, Df

    def stiff_nbr_frac(self):
        stiff = (self.state > 0).astype(float)
        return np.array([stiff[nb].mean() for nb in self.nbrs])

    def impedance(self, C, f):
        """Leaf flow fractions and input impedances for breathing frequencies f (1/min)."""
        p, lv, leaf = self.p, self.levels, self.leaf
        w = 2 * np.pi * np.atleast_1d(f)[:, None] / 60.0
        Z = np.empty((w.shape[0], self.n), complex); Ys = np.empty_like(Z)
        Z[:, leaf] = self.R_aw[leaf] + self.R_unit + 1 / (1j * w * C[None, :])
        for k in reversed(self.inner_levels):
            Ys[:, k] = 1 / Z[:, self.c1[k]] + 1 / Z[:, self.c2[k]]
            Z[:, k] = self.R_aw[k] + 1 / Ys[:, k]
        phi = np.ones_like(Z)
        for k in self.inner_levels:
            for c_ in (self.c1[k], self.c2[k]):
                phi[:, c_] = phi[:, k] * (1 / Z[:, c_]) / Ys[:, k]
        Zl = Z[:, 0]
        Zrs = Zl + p.R_cw + 1 / (1j * w[:, 0] * p.C_cw)
        return Zl, Zrs, phi[:, leaf], w[:, 0]

    def mech_tissue(self):
        """Mechanics of the open unit (recruitment state ignored)."""
        rd = self.rd; self.rd = np.zeros(self.N, int)
        out = self.mech(); self.rd = rd
        return out

    def breathe(self, VA):
        if self.p.parenchyma != "surfactant":
            return self._breathe(VA)
        if self.p.breath_model == "cycle":
            from .breath import breath_cycle
            for _ in range(2):
                air = self._breathe(VA)
                cyc = breath_cycle(self, max(air["VT"] - self.p.VD, 1e-4), air["f"])
                self.rd = cyc["rd"]
            self.x = cyc["x"]
            _, V0t, _, _ = self.mech_tissue()
            VTi = cyc["VTi"]
            air["VTi"] = VTi
            VTi = np.where(self.flood, 0.0, VTi)
            air["VAi"] = VA * VTi / max(VTi.sum(), 1e-12)
            air["strain"] = VTi / V0t
            air["events"] = np.where(self.state == 0, cyc["events"], 0.0)
            air["open_frac"] = cyc["open_frac"]
            return air
        # recruitment state and tidal distribution depend on each other: iterate to a fixed point
        for _ in range(12):
            air = self._breathe(VA)
            new = self.classify_rd(air)
            if (new == self.rd).all():
                break
            self.rd = new
        air = self._breathe(VA)
        vent = ((self.rd != 2) | (self.state != 0)) & ~self.flood
        air["VAi"] = np.where(vent, air["VAi"], 0.0)
        air["VAi"] = air["VAi"] * VA / max(air["VAi"].sum(), 1e-12)
        # cyclic units are strained by reopening every breath
        air["strain"] = np.where((self.rd == 1) & (self.state == 0), air["strain"] * self.p.cyc_strain, air["strain"])
        return air

    # ------------------------------------------------------------ surfactant, recruitment/derecruitment
    def _init_parenchyma(self):
        p = self.p
        self.sq = np.ones(self.N)                          # surfactant quality: 1 normal, 0 absent
        self.rd = np.zeros(self.N, int)                    # 0 open, 1 cyclic (closes and reopens each breath), 2 closed
        self.z_mech = (self.y - p.heart_frac) * p.lung_height_cm       # cm, for the pleural gradient
        rng = np.random.default_rng(p.seed + 12345)
        self.at2_rate = np.exp(p.at2_sigma * rng.standard_normal(self.N) - p.at2_sigma ** 2 / 2)
        self.x = np.ones(self.N)                          # Bates-Irvin recruitment variable (1 open)
        self.flood = np.zeros(self.N, bool)               # flooded units (diffuse alveolar damage)
        self.sen = np.zeros(self.N)                       # AT2 senescence burden (0-1)
        self.aberrant = np.zeros(self.N, bool)            # KRT5-/KRT17+ aberrant basaloid epithelium present
        self._senolytic_done = False
        self._rem_budget = None                           # remodelling budget per slow step (kappa x time)
        self.aging_on = False

    def critical_pressures(self):
        """Closing and opening transpulmonary pressures (cmH2O) from surfactant quality
        (Bates-Irvin units; opening above closing gives hysteresis)."""
        p = self.p
        Pc = p.Pc_off + p.Pc_slope * (1 - self.sq)
        Po = Pc + p.Po_off + p.Po_slope * (1 - self.sq)
        return Pc, Po

    def classify_rd(self, air):
        p = self.p
        Pc, Po = self.critical_pressures()
        Pexp = p.PL_FRC + p.PEEP + p.pl_grad * self.z_mech               # regional end-expiratory P_L
        Copen = self.C_h                                                  # tissue compliance of the open unit
        Pins = Pexp + np.where(self.rd == 2, 0.0, air["VTi"] / Copen)
        Pins_open = Pexp + air["VTi"] / Copen
        new = self.rd.copy()
        was_closed = self.rd == 2
        # closed units reopen only if the inspiratory pressure they would feel reaches Po
        reopen = was_closed & (Pexp + air["VT"] / self.N / Copen * 1.0 >= Po)
        closes = (~was_closed) & (Pexp <= Pc)
        new[(~was_closed) & (Pexp > Pc)] = 0
        new[closes & (Pins_open >= Po)] = 1
        new[closes & (Pins_open < Po)] = 2
        new[reopen] = np.where(Pexp[reopen] > Pc[reopen], 0, 1)
        new[self.state != 0] = 0
        return new

    def _breathe(self, VA):
        """Choose f minimising inspiratory mechanical power at alveolar ventilation VA (L/min)."""
        p = self.p
        C, V0, _, _ = self.mech()
        f = np.linspace(*p.f_grid[:2], int(p.f_grid[2]))
        Zl, Zrs, phi, w = self.impedance(C, f)
        VT = VA / f + p.VD
        Eeff = -w * Zrs.imag
        Qh = w * VT / 2
        power = f / 60 * Eeff * VT ** 2 / 2 + 0.5 * Zrs.real * Qh ** 2
        i = int(np.argmin(power))
        a = np.abs(phi[i]); share = a / a.sum()
        VTi = VT[i] * a
        EL = 1 / C.sum()
        return dict(f=f[i], VT=VT[i], VTi=VTi, VAi=VA * share, strain=VTi / V0,
                    EL=EL, dPL=VT[i] * EL, Zl=abs(Zl[i]), power=power[i],
                    FRC=V0.sum(), C=C.sum())

    # ------------------------------------------------------------ gas exchange
    def content(self, P):
        p = self.p
        S = P ** p.hill_n / (P ** p.hill_n + p.P50 ** p.hill_n)
        return 1.34 * p.Hb * S + 0.003 * P                         # mL O2 / dL

    def inv_content(self, C):
        lo = np.full(np.shape(C), 0.1); hi = np.full(np.shape(C), 800.0)
        for _ in range(60):
            m = 0.5 * (lo + hi); up = self.content(m) < C
            lo = np.where(up, m, lo); hi = np.where(up, hi, m)
        return 0.5 * (lo + hi)

    def rf_DL(self, D, Vcap):
        """Unit diffusing capacity (healthy = 1) from membrane (D) and capillary volume (Vcap)."""
        vc = Vcap / (self.p.Vcap / self.N)
        with np.errstate(divide="ignore"):
            return np.where((D > 0) & (vc > 0), 1.0 / (1.0 / np.maximum(2 * D, 1e-300) + 1.0 / np.maximum(2 * vc, 1e-300)), 0.0)

    def gas(self, Ql, VAi, Vcap, D, VO2, Cv=None, maxit=300):
        p = self.p
        PI = p.FIO2 * (p.PB - 47.0)
        Qt = Ql.sum()
        Cv = self.content(40.0) if Cv is None else Cv
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(Ql > 0, Vcap * 60.0 / np.maximum(Ql, 1e-300), 0.0)
        if p.rf:
            # Roughton-Forster per unit: 1/DL = 1/Dm + 1/(theta Vc), both normalised to 1/2 of a healthy
            # unit so that DL = 1 in health; the equilibration exponent is DL x (healthy transit) / t_eq
            DLu = self.rf_DL(D, Vcap)
            v0 = p.Vcap / self.N
            kt = np.where(Ql > 0, 3.0 / p.t_eq * DLu * v0 * 60.0 / np.maximum(Ql, 1e-300), 0.0)
        else:
            kt = 3.0 / p.t_eq * D * t
        ok = True; nfail = 0
        for it in range(maxit):
            Pv = float(self.inv_content(np.array(Cv)))
            lo = np.full(self.N, Pv); hi = np.full(self.N, PI)
            for _ in range(50):
                PA = 0.5 * (lo + hi)
                Pc = PA - (PA - Pv) * np.exp(-kt)
                F = VAi * 1000 * (PI - PA) / 863.0 - Ql * 10 * (self.content(Pc) - Cv)
                lo = np.where(F > 0, PA, lo); hi = np.where(F > 0, hi, PA)
            PA = 0.5 * (lo + hi)
            PA = np.where(Ql > 0, PA, PI)
            Pc = PA - (PA - Pv) * np.exp(-kt)
            Cc = self.content(Pc)
            Ca = float((Ql * Cc).sum() / Qt)
            Cv_new = Ca - VO2 / (10 * Qt)
            if Cv_new < 0.5:
                Cv_new, ok = 0.5, False; nfail += 1
                if nfail > 5:
                    Cv = Cv_new; break
            if abs(Cv_new - Cv) < 1e-7:
                Cv = Cv_new; break
            Cv = 0.5 * Cv + 0.5 * Cv_new
        Pa = float(self.inv_content(np.array(Ca)))
        Sa = Pa ** p.hill_n / (Pa ** p.hill_n + p.P50 ** p.hill_n)
        VO2_air = float((VAi * 1000 * (PI - PA) / 863.0).sum())
        return dict(PA=PA, Pc=Pc, Ca=Ca, Cv=Cv, PaO2=Pa, SaO2=Sa, Pv=Pv, ok=ok,
                    VO2_air=VO2_air, VO2_blood=Qt * 10 * (Ca - Cv), t_transit=t)

    # ------------------------------------------------------------ equilibrium
    def hpv_factor(self, PA):
        p = self.p
        if not p.hpv:
            return np.ones_like(PA)
        f = lambda x: 1 + p.hpv_max / (1 + (np.maximum(x, 1e-6) / p.hpv_P50) ** p.hpv_n)
        return f(PA) / f(self._PA_ref)

    def _baseline(self):
        p, N = self.p, self.N
        self.R_lesion = np.zeros(N)
        self._Pcap0 = np.full(N, p.LAP + p.dP0 * (p.split[1] / 2 + p.split[2]))
        zeros = dict(Pm_a=np.full(self.n, p.LAP + p.dP0), Pm_v=np.full(len(self.r_ref_v), p.LAP),
                     Pm_ta=np.full(N, p.LAP + p.dP0 * (1 - p.split[0] / 2)), Pcap=self._Pcap0.copy(),
                     Pte=np.full(N, p.LAP + p.dP0 * (p.split[1] + p.split[2])), Pvs=np.full(N, p.LAP + p.dP0 * p.split[2]))
        # solve with the reference radii held fixed (alpha applied to r_u = r_ref)
        a = p.alpha
        p.alpha = 0.0
        if p.perf_grav:
            self._calibrate_west(zeros)
        self._P_base = zeros
        vs = self.solve_vessels(self.r_ref, self.r_ref_v, np.full(N, self.r_ta_ref), np.ones(N), np.ones(N), p.Q)
        p.alpha = a
        self._Pcap0 = vs["Pcap"].copy()
        Pb = vs["P"]
        self._P_base = Pb
        self.ru_a = self.r_ref / (1 + a * self._ptm(Pb["Pm_a"], self.za_edge))
        self.ru_v = self.r_ref_v / (1 + a * self._ptm(Pb["Pm_v"], self.zv_edge))
        self.ru_ta = self.r_ta_ref / (1 + a * self._ptm(Pb["Pm_ta"], self.z_leaf))
        self.ru0 = (self.ru_a.copy(), self.ru_v.copy(), self.ru_ta.copy())
        b = p.b
        Qa = vs["Q"] * LPM; Qt = vs["Ql"] * LPM
        if p.law == "shared":
            self.tau0_a = np.full(self.n, self.tau0_shared)
            self.tau0_v = np.full(len(self.r_ref_v), self.tau0_shared)
            rt = self.r_ta_ref           # lumped sub-tree: at rest at baseline by construction (equivalent vessel)
            self.tau0_t = (4 * p.mu_b * Qt / (np.pi * rt ** 3)) / (rt ** (b - 1) * p.L_term ** ((b - 1) / 2)) * \
                self.pulse_factor(np.full(self.N, rt), True)
        else:
            self.tau0_a = (4 * p.mu_b * Qa / (np.pi * self.r_ref ** 3)) / (self.r_ref ** (b - 1) * self.Lv ** ((b - 1) / 2))
            if p.pulsatile:
                self.tau0_a = self.tau0_a * self.pulse_factor(self.r_ref, True, self.depth)
            rv_, Qv = self.r_ref_v, vs["Qv"] * LPM
            self.tau0_v = (4 * p.mu_b * Qv / (np.pi * rv_ ** 3)) / (rv_ ** (b - 1) * self.Lvv ** ((b - 1) / 2))
            if p.pulsatile:
                self.tau0_v = self.tau0_v * self.pulse_factor(rv_, False, self.vdepth)
            rt = self.r_ta_ref
            self.tau0_t = (4 * p.mu_b * Qt / (np.pi * rt ** 3)) / (rt ** (b - 1) * p.L_term ** ((b - 1) / 2))
            if p.pulsatile:
                self.tau0_t = self.tau0_t * self.pulse_factor(np.full(self.N, rt), True, np.full(self.N, p.G + 1))
        self._PA_ref = 100.0
        air = self.breathe(p.VA)
        C, V0, cf, Df = self.mech()
        g = self.gas(vs["Ql"], air["VAi"], vs["Vcap"], Df, p.VO2)
        self._PA_ref = float(np.median(g["PA"]))
        self.s_ref = float(np.median(air["strain"]))
        self.base = dict(vs=vs, air=air, gas=g)
        self.Ees, self.EDV = p.Ees0, p.EDV0
        self.base_state = self.equilibrate()
        self._Rta0 = float(self._vres(np.array([self.r_ta_ref]), p.L_term)[0])
        self.DL0 = self.base_state["DL"]
        self.Dm0, self.Vc0 = self._Dm, self._Vc
        # aging: imperfect set-point (inward offset age_M, per-vessel error age_sigma), active after baseline
        rng = np.random.default_rng(p.seed + 999)
        mk = lambda n: (1 + p.age_M) * np.exp(p.age_sigma * rng.standard_normal(n))
        self.age_a, self.age_v, self.age_t = mk(len(self.ru_a)), mk(len(self.ru_v)), mk(len(self.ru_ta))
        self.aging_on = bool(p.age_M > 0 or p.age_sigma > 0)
        self.PVV0 = self.pvv(self.base_state["vs"])

    def subtree_mean(self, leafval):
        """Mean of a leaf quantity over the subtree below every node (symmetric tree)."""
        v = np.empty(self.n); v[self.leaf] = leafval
        for k in reversed(self.inner_levels):
            a, b_ = self.c1[k], self.c2[k]
            v[k] = (self.nleaf[a] * v[a] + self.nleaf[b_] * v[b_]) / self.nleaf[k]
        return v

    def tether(self):
        """Traction on extra-alveolar conduits by stiff parenchyma.

        Each unit's elastic recoil at FRC scales with V0/C; the rise above its healthy
        value lowers perivascular pressure around the conduits that run through that
        region, i.e. adds to their transmural pressure. A conduit feels the mean over the
        units it supplies. P_el_FRC is the healthy recoil at FRC (cmH2O)."""
        p = self.p
        if p.traction <= 0:
            return 0.0
        C, V0, _, _ = self.mech()
        ratio = (V0 / C) / (self.V_h / self.C_h)
        dP = p.traction * p.P_el_FRC * (ratio - 1.0) / 1.36          # mmHg
        return self.subtree_mean(dP)

    def maladaptive(self, P):
        """Set-point shift M >= 0 on arterial vessels: rest point becomes z = 1 + M.

        Pressure drive: eta_p * max(0, Ptm/Ptm_healthy - 1)  (stretch-induced wall thickening)
        Scar drive:     eta_s * fraction of stiff units supplied (fibrotic micro-environment)
        Applied to terminal arterioles and to arterial conduits of depth >= mal_depth."""
        p = self.p
        if p.eta_p <= 0 and p.eta_s <= 0:
            return 0.0, 0.0
        Pb = self._P_base
        stiff = (self.state > 0).astype(float)
        Mt = p.eta_p * np.maximum(0, P["Pm_ta"] / Pb["Pm_ta"] - 1) + p.eta_s * stiff
        Ma = p.eta_p * np.maximum(0, P["Pm_a"] / Pb["Pm_a"] - 1) + p.eta_s * self.subtree_mean(stiff)
        Ma = np.where(self.depth >= p.mal_depth, Ma, 0.0)
        return Ma, Mt

    # ------------------------------------------------------------ right ventricle
    def rv_state(self, R, Qd, exercise=False):
        """Right ventricle as an end-systolic elastance coupled to the pulmonary network.

        ESPVR Pes = Ees (EDV - SV - V0), load Pes = LAP + HR R SV (R = network resistance in
        mmHg min/mL, Pes ~ mPAP), so SV = (Ees (EDV - V0) - LAP) / (Ees + HR R), exactly.
        Rest (chronic): homeometric adaptation first, Ees -> clip(ratio * Ea, Ees0, Ees_max),
        Ea = Pes/SV at the demanded SV; then heterometric: EDV takes the value that delivers the
        demanded SV, within [EDV_min, EDV_max]. Exercise (acute): no adaptation; inotropic reserve
        Ees x ex_inotropy for a normal RV, shrinking linearly to none as chronic Ees reaches
        Ees_max (exhausted contractile reserve); preload reserve EDV + ex_EDV (capped); HR_ex.
        Delivered output is min(demand, capacity). RAP follows an exponential EDPVR."""
        p = self.p
        if exercise:
            HR = p.HR_ex
            used = (self.Ees - p.Ees0) / (p.Ees_max - p.Ees0)       # fraction of hypertrophic reserve used
            Ees = self.Ees * (1 + (p.ex_inotropy - 1) * (1 - float(np.clip(used, 0, 1))))
            EDVcap = min(self.EDV + p.ex_EDV, p.EDV_max)
        else:
            HR = p.HR
            SVd = Qd * 1000 / HR
            Pes_d = p.LAP + HR * R * SVd
            Ees = float(np.clip(p.rv_ratio * Pes_d / SVd, p.Ees0, p.Ees_max))
            EDVcap = float(np.clip(p.V0_rv + Pes_d / Ees + SVd, p.EDV_min, p.EDV_max))
        SVcap = max(0.0, (Ees * (EDVcap - p.V0_rv) - p.LAP) / (Ees + HR * R))
        Q = min(Qd, HR * SVcap / 1000)
        SV = Q * 1000 / HR
        Pes = p.LAP + HR * R * SV
        EDV = min(EDVcap, p.V0_rv + Pes / Ees + SV)
        RAP = p.RAP0 * np.exp(p.k_edpvr * (EDV - p.EDV0))
        return dict(Q=Q, Qcap=HR * SVcap / 1000, SV=SV, Ees=Ees, EDV=EDV, Ea=Pes / SV if SV > 0 else np.inf,
                    coupling=Ees / (Pes / SV) if SV > 0 else 0.0, RAP=RAP, HR=HR,
                    limited=bool(HR * SVcap / 1000 < Qd * (1 - 1e-9)))

    def equilibrate(self, Qtot=None, VO2=None, VA=None, remodel=None, tol=1e-5, maxit=None, exercise=False):
        """Quasi-static equilibrium of perfusion, HPV, gas exchange, (optionally) remodelling and,
        with p.rv, right-ventricular output. Qtot is then the demanded output."""
        p = self.p
        maxit = p.eq_maxit if maxit is None else maxit
        Qtot = p.Q if Qtot is None else Qtot
        Qd = Qtot; rv = None; dq = 0.0
        VO2 = p.VO2 if VO2 is None else VO2
        VA = p.VA if VA is None else VA
        remodel = p.remodel if remodel is None else remodel
        air = self.breathe(VA)                       # breathing first: it updates the recruitment state
        C, V0, cf, Df = self.mech()
        h = self.h.copy()
        rem_left = self._rem_budget if self._rem_budget is not None else 0.0
        P = None; Cv = None
        ru_a, ru_v, ru_ta = self.ru_a.copy(), self.ru_v.copy(), self.ru_ta.copy()
        floor_a, floor_v, floor_t = 1e-3 * self.ru0[0], 1e-3 * self.ru0[1], 1e-3 * self.ru0[2]
        Pt = self.tether()
        for it in range(maxit):
            vs = self.solve_vessels(ru_a, ru_v, ru_ta, cf, h, Qtot, P=P, Pt=Pt)
            P = vs["P"]
            g = self.gas(vs["Ql"], air["VAi"], vs["Vcap"], Df, VO2, Cv=Cv)
            Cv = g["Cv"]
            h_new = self.hpv_factor(g["PA"])
            dh = np.abs(h_new - h).max()
            h = 0.5 * h + 0.5 * h_new
            if p.rv:
                rv = self.rv_state(vs["dP"] / Qtot / 1000.0, Qd, exercise=exercise)
                dq = abs(rv["Q"] - Qtot)
                Qtot = max(0.5 * Qtot + 0.5 * rv["Q"], 1e-3)
            zerr = 0.0
            if remodel and (p.remodel_steps <= 0 or it < p.remodel_steps):
                za, zv, zt = self.shear_z(vs)
                Ma, Mt = self.maladaptive(vs["P"])
                za = za / (1 + Ma); zt = zt / (1 + Mt)
                if self.aging_on:
                    za = za / self.age_a; zv = zv / self.age_v; zt = zt / self.age_t
                live_a = vs["Q"] > 1e-9; live_t = vs["Ql"] > 1e-9; live_v = vs["Qv"] > 1e-9
                if not p.remodel_conduits:
                    za = np.ones_like(za); zv = np.ones_like(zv)
                zerr = max(np.abs(za[live_a] - 1).max(initial=0), np.abs(zv[live_v] - 1).max(initial=0),
                           np.abs(zt[live_t] - 1).max(initial=0))
                kap = p.kappa
                if self._rem_budget is not None:
                    kap = min(p.kappa, rem_left); rem_left -= kap
                st = lambda z: np.clip(kap * (z - 1), -p.clip, p.clip)
                ru_a = np.maximum(ru_a * np.exp(st(za)), floor_a)
                ru_v = np.maximum(ru_v * np.exp(st(zv)), floor_v)
                ru_ta = np.maximum(ru_ta * np.exp(st(zt)), floor_t)
                if self._rem_budget is not None and rem_left <= 1e-12:
                    zerr = 0.0                 # remodelling budget of this slow step is spent
            if dh < tol and zerr < tol and dq < 1e-6 and vs["err"] < 1e-6:
                break
            if not g["ok"] and it > 200:
                break                      # O2 demand cannot be met: report the failed state
            if vs["PPA"] > 1.5 * p.stop_mPAP and it > 50:
                break                      # decompensated: stop refining a runaway state
        if remodel:
            self.ru_a, self.ru_v, self.ru_ta = ru_a, ru_v, ru_ta
        if rv is not None and not exercise:
            self.Ees, self.EDV = rv["Ees"], rv["EDV"]
        self.h = h
        vent = air["VAi"] > 0
        if p.rf:
            DL = float((self.rf_DL(Df, vs["Vcap"]) * vent).sum())
            self._Dm, self._Vc = float((Df * vent).sum()), float((vs["Vcap"] / (p.Vcap / self.N) * vent).sum())
        else:
            DL = float((Df * vs["Vcap"] * vent).sum())
            self._Dm = self._Vc = np.nan
        return dict(vs=vs, air=air, gas=g, h=h, it=it, zerr=zerr, dh=dh, DL=DL,
                    cf=cf, C=C, V0=V0, Pt=Pt, rv=rv, Q=Qtot)

    def functional_E(self, vs):
        """Murray/Hu-Cai functional of the remodelable vessels: dissipation of the whole circuit plus
        lambda (pi r^2 L)^b, with lambda per vessel from its set-point constant: tau0^2 = b lambda mu pi^(b-1)."""
        p, b = self.p, self.p.b
        Q = vs["Q"] * LPM; Ql = vs["Ql"] * LPM; Qv = vs["Qv"] * LPM
        mu = p.mu_b
        D = (vs["dP"] * MMHG) * (p.Q * LPM)                 # total power = inflow x pressure drop (W)
        cost = 0.0
        for tau0, r, L in ((self.tau0_a, vs["r_a"], self.Lv), (self.tau0_v, vs["r_v"], self.Lvv),
                           (self.tau0_t, vs["r_ta"], np.full(self.N, p.L_term))):
            lam = tau0 ** 2 / (b * mu * np.pi ** (b - 1))
            cost += float((lam * (np.pi * r ** 2 * L) ** b).sum())
        return D + cost

    def pvv(self, vs, cut=None):
        """Conduit vessel volume (arteries + veins), large = depth <= cut[0], small = depth >= cut[1].

        Default for G = 10: large = generations 0-4, small = 7-10 (the last four)."""
        d = self.depth
        cut = cut or (min(4, self.p.G // 2), self.p.G - 3)
        va = np.pi * vs["r_a"] ** 2 * self.Lv * 1e6                      # mL
        vv = np.pi * vs["r_v"] ** 2 * self.Lvv * 1e6
        dv = self.vdepth
        return dict(large=float(va[d <= cut[0]].sum() + vv[dv <= cut[0]].sum()),
                    small=float(va[d >= cut[1]].sum() + vv[dv >= cut[1]].sum()),
                    total=float(va.sum() + vv.sum()))
