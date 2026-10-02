"""Independent (septal) pulmonary venous tree.

Anatomy: pulmonary arteries run with the airways in the centre of each unit; veins run in the
interlobular septa, between units. Here the venous tree is built on the corners of the unit grid
((nx+1) x (ny+1) venules) by recursive median bisection, and every unit drains into its four
corner venules with equal share. A venule therefore collects from up to four neighbouring units,
which couples neighbours on the outflow side. The network is no longer series-parallel and is
solved as a graph Laplacian (Kirchhoff nodal analysis) with a sparse direct solver, as in circuit
and microvascular network codes.

  veins = "septal"  corner venules, interdigitated with the arterial H-tree (default of this class)
  veins = "mirror"  the same general solver on the mirrored tree: must reproduce Lung exactly (C13)

vein_occl in [0, 1): a corner venule loses conductance in proportion to the fraction of stiff
units around it, g -> g (1 - vein_occl f_stiff) (septal fibrosis compressing or obliterating veins).
"""
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from .lung import Lung, pois, MMHG, LPM


class SeptalLung(Lung):

    # ------------------------------------------------------------ venous topology
    def _venous_topology(self):
        p = self.p
        if p.veins in ("mirror", "mirror_general"):
            self.vpar = self.parent.copy()
            self.vdepth = self.depth.copy()
            self.vleaf = self.leaf.copy()
            self.inc_cell = np.arange(self.N)
            self.inc_vleaf = np.arange(self.N)                 # position in vleaf
            self.inc_share = np.ones(self.N)
            return
        cx, cy = np.meshgrid(np.arange(self.nx + 1), np.arange(self.ny + 1), indexing="ij")
        pts = np.c_[cx.ravel(), cy.ravel()].astype(float)
        par, dep, leaf_of_pt = [], [], np.empty(len(pts), int)

        def build(idx, parent, d):
            me = len(par); par.append(parent); dep.append(d)
            if len(idx) == 1:
                leaf_of_pt[idx[0]] = me
                return
            sub = pts[idx]
            ax = int(np.ptp(sub[:, 1]) > np.ptp(sub[:, 0]))
            order = idx[np.argsort(sub[:, ax], kind="stable")]
            h = len(order) // 2
            build(order[:h], me, d + 1); build(order[h:], me, d + 1)

        build(np.arange(len(pts)), -1, 0)
        self.vpar = np.array(par); self.vdepth = np.array(dep)
        is_leaf = np.ones(len(par), bool); is_leaf[self.vpar[self.vpar >= 0]] = False
        self.vleaf = np.where(is_leaf)[0]
        pos = -np.ones(len(par), int); pos[self.vleaf] = np.arange(len(self.vleaf))
        cid = lambda i, j: pos[leaf_of_pt[i * (self.ny + 1) + j]]
        cells, vls = [], []
        for c in range(self.N):
            i, j = self.ix[c], self.iy[c]
            for di, dj in ((0, 0), (1, 0), (0, 1), (1, 1)):
                cells.append(c); vls.append(cid(i + di, j + dj))
        self.inc_cell = np.array(cells); self.inc_vleaf = np.array(vls)
        self.inc_share = np.full(len(cells), 0.25)

    def _build_vessels(self):
        super()._build_vessels()
        p, b = self.p, self.p.b
        self._venous_topology()
        nv = self.nv = len(self.vpar)
        # children lists by depth for bottom-up sums
        self.vorder = np.argsort(-self.vdepth, kind="stable")
        self.Lvv = p.L_pa * 2.0 ** (-self.vdepth / 3.0)
        # healthy venous flows: every unit carries Q/N, split by share to its venules
        f = np.zeros(nv)
        np.add.at(f, self.vleaf[self.inc_vleaf], self.inc_share * p.Q / self.N)
        self.v_cells = np.zeros(nv)
        np.add.at(self.v_cells, self.vleaf[self.inc_vleaf], self.inc_share)
        for k in self.vorder:
            if self.vpar[k] >= 0:
                f[self.vpar[k]] += f[k]; self.v_cells[self.vpar[k]] += self.v_cells[k]
        tau0 = 4 * p.mu_b * p.Q * LPM / (np.pi * p.r_pa ** (2 + b) * p.L_pa ** ((b - 1) / 2))
        self.r_ref_v = (4 * p.mu_b * f * LPM / (np.pi * tau0 * self.Lvv ** ((b - 1) / 2))) ** (1 / (2 + b))
        zc = self.vsubtree_mean(self.z_leaf)
        zpar = np.where(self.vpar >= 0, zc[np.maximum(self.vpar, 0)], 0.0)
        self.zv_edge = 0.5 * (zc + zpar)
        self._nodes()

    def _nodes(self):
        n, N, nv = self.n, self.N, self.nv
        self.iPA = 0
        self.iA = 1 + np.arange(n)
        self.iO = 1 + n + np.arange(N)
        self.iV = 1 + n + N + np.arange(nv)
        self.nnodes = 1 + n + N + nv
        # static endpoints
        self.ea_up = np.where(self.parent >= 0, self.iA[np.maximum(self.parent, 0)], self.iPA)
        self.ea_dn = self.iA
        self.el_up = self.iA[self.leaf]; self.el_dn = self.iO
        self.ek_up = self.iO[self.inc_cell]; self.ek_dn = self.iV[self.vleaf[self.inc_vleaf]]
        self.ev_up = self.iV
        self.ev_dn = np.where(self.vpar >= 0, self.iV[np.maximum(self.vpar, 0)], -1)   # -1 = LA

    def corner_stiff(self):
        stiff = (self.state > 0).astype(float)
        fs = np.zeros(len(self.vleaf)); w = np.zeros(len(self.vleaf))
        np.add.at(fs, self.inc_vleaf, self.inc_share * stiff[self.inc_cell])
        np.add.at(w, self.inc_vleaf, self.inc_share)
        return fs / np.maximum(w, 1e-12)

    def vsubtree_mean(self, cellval):
        s = np.zeros(self.nv)
        np.add.at(s, self.vleaf[self.inc_vleaf], self.inc_share * cellval[self.inc_cell])
        for k in self.vorder:
            if self.vpar[k] >= 0:
                s[self.vpar[k]] += s[k]
        return s / np.maximum(self.v_cells, 1e-12)

    def tether(self):
        p = self.p
        if p.traction <= 0:
            return 0.0
        C, V0, _, _ = self.mech()
        ratio = (V0 / C) / (self.V_h / self.C_h)
        dP = p.traction * p.P_el_FRC * (ratio - 1.0) / 1.36
        return (self.subtree_mean(dP), self.vsubtree_mean(dP))

    # ------------------------------------------------------------ solver
    def _baseline(self):
        # calibrate the venule resistance so that the healthy mPAP - LAP is dP0 with this venous tree
        p, N = self.p, self.N
        a, ac = p.alpha, p.alpha_c
        p.alpha = p.alpha_c = 0.0
        self._Pcap0 = np.full(N, p.LAP + 3.0)
        self._P_base = dict(Pm_a=np.full(self.n, p.LAP + p.dP0), Pm_v=np.full(self.nv, p.LAP),
                            Pm_ta=np.full(N, p.LAP + p.dP0), Pcap=self._Pcap0.copy(),
                            Pte=np.full(N, p.LAP + 3.6), Pvs=np.full(N, p.LAP + 1.5))
        for _ in range(0 if p.perf_grav else 4):
            vs = self.solve_vessels(self.r_ref, self.r_ref_v, np.full(N, self.r_ta_ref), np.ones(N), np.ones(N), p.Q)
            self.R_vn += N * (p.dP0 - vs["dP"]) / p.Q
        p.alpha, p.alpha_c = a, ac
        assert self.R_vn > 0
        super()._baseline()

    def solve_vessels(self, ru_a, ru_v, ru_ta, c, h, Qtot, P=None, tol=1e-7, maxit=200, Pt=0.0):
        p, N = self.p, self.N
        if P is None:
            P = self._P_base
        else:
            maxit = min(maxit, 3)      # warm start: the outer equilibrium loop finishes the distension fixed point
        Pt_a, Pt_v = Pt if isinstance(Pt, tuple) else (Pt, Pt)
        occl = 1.0 - p.vein_occl * self.corner_stiff() if p.vein_occl > 0 else 1.0
        occl_k = occl[self.inc_vleaf] if np.ndim(occl) else 1.0
        rows_up = np.r_[self.ea_up, self.el_up, self.ek_up, self.ev_up]
        rows_dn = np.r_[self.ea_dn, self.el_dn, self.ek_dn, self.ev_dn]
        la = rows_dn < 0
        for it in range(maxit):
            r_a = ru_a * (1 + p.alpha * (self._ptm(P["Pm_a"], self.za_edge) + Pt_a))
            r_v = ru_v * (1 + p.alpha * (self._ptm(P["Pm_v"], self.zv_edge) + Pt_v))
            r_ta = ru_ta * (1 + p.alpha * self._ptm(P["Pm_ta"], self.z_leaf))
            Ra = self._vres(r_a, self.Lv); Rv = self._vres(r_v, self.Lvv)
            Rta = self._vres(r_ta, p.L_term) * h
            Gcap = self._gcap(c, P)
            gl = np.where(Gcap > 0, 1.0 / (Rta + 1.0 / np.maximum(Gcap, 1e-300)), 0.0)
            gk = self.inc_share * occl_k / self.R_vn
            g = np.r_[1 / Ra, gl, gk, 1 / Rv]
            up, dn = rows_up, np.where(la, 0, rows_dn)
            inner = ~la
            i = np.r_[up[inner], dn[inner], up[inner], dn[inner], up[la]]
            j = np.r_[up[inner], dn[inner], dn[inner], up[inner], up[la]]
            v = np.r_[g[inner], g[inner], -g[inner], -g[inner], g[la]]
            M = csc_matrix((v, (i, j)), shape=(self.nnodes, self.nnodes))
            M = M + csc_matrix((np.full(self.nnodes, 1e-15 * g.max()), (np.arange(self.nnodes), np.arange(self.nnodes))),
                               shape=M.shape)
            rhs = np.zeros(self.nnodes); rhs[self.iPA] += Qtot
            np.add.at(rhs, up[la], g[la] * p.LAP)
            phi = spsolve(M, rhs)
            Pn = lambda idx: np.where(idx >= 0, phi[np.maximum(idx, 0)], p.LAP)
            PPA = phi[self.iPA]
            Q = (Pn(self.ea_up) - Pn(self.ea_dn)) / Ra
            PA_leaf = phi[self.iA[self.leaf]]; PO = phi[self.iO]
            Ql = gl * (PA_leaf - PO)
            Qv = (Pn(self.ev_up) - Pn(self.ev_dn)) / Rv
            Pta_end = PA_leaf - Ql * Rta
            Pcap_n = np.where(Gcap > 0, Pta_end - 0.5 * Ql / np.maximum(Gcap, 1e-300), PA_leaf)
            new = dict(Pm_a=0.5 * (Pn(self.ea_up) + Pn(self.ea_dn)),
                       Pm_v=0.5 * (Pn(self.ev_up) + Pn(self.ev_dn)),
                       Pm_ta=0.5 * (PA_leaf + Pta_end), Pcap=Pcap_n, Pte=Pta_end, Pvs=PO)
            err = max(np.abs(new[key] - P[key]).max() for key in new)
            P = {key: 0.6 * new[key] + 0.4 * P[key] for key in new}
            if err < tol:
                break
        Vcap = self._vcap(c, P)
        return dict(Q=Q, Qv=Qv, Ql=Ql, PPA=PPA, dP=PPA - p.LAP, Req=(PPA - p.LAP) / Qtot, r_a=r_a, r_v=r_v,
                    r_ta=r_ta, P=P, Vcap=Vcap, it=it, err=err, Pcap=P["Pcap"], PO=PO)


def make_lung(p):
    if p.geometry == "anatomical":
        from .anatomy import AnatomicalLung
        return AnatomicalLung(p)
    return SeptalLung(p) if p.veins in ("septal", "mirror_general") else Lung(p)
