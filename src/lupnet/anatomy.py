"""Anatomical 3D lung: an asymmetric branching tree grown into a lung-shaped volume.

Volume-filling algorithm in the spirit of Tawhai et al. (2000) and Kitaoka et al. (1999): N terminal
units (acinus clusters) are seeded uniformly in the domain; starting from the trachea, each set of units
is split by the plane through its centroid normal to its longest principal axis, and each child branch
grows from the parent end 40% of the way towards the centroid of its half. The result is an asymmetric
binary tree whose depth varies from branch to branch; branch lengths and positions are geometric.
Pulmonary arteries run with the airways (same tree); veins mirror the arteries. Everything else in the
model (set-point remodelling, gravity, zones, surfactant, recruitment, gas exchange, RV) works on it
unchanged, through the general binary-tree solver.

Domain: one lung modelled as a half-ellipsoid dome (base flat), height H, lateral half-width a, antero-
posterior half-depth c (cm), with a medial concavity for the heart. y is height (0 base, 1 apex).
"""
import numpy as np
from scipy.spatial import cKDTree
from .lung import Lung, pois, CMH2O


class AnatomicalLung(Lung):

    def _inside(self, P):
        p = self.p
        a, c, H = p.anat_a, p.anat_c, p.anat_H
        x, y, z = P[:, 0], P[:, 1], P[:, 2]
        dome = (x / a) ** 2 + (z / c) ** 2 + (y / H) ** 2 <= 1.0
        heart = ((x + a) / (0.55 * a)) ** 2 + ((z - 0.2 * c) / (0.45 * c)) ** 2 + (y / (0.45 * H)) ** 2 <= 1.0
        return dome & (y >= 0) & ~heart

    def _topology(self):
        p = self.p
        rng = np.random.default_rng(p.anat_seed)
        N = p.anat_N
        pts = np.empty((0, 3))
        while len(pts) < N:
            c = rng.uniform([-p.anat_a, 0, -p.anat_c], [p.anat_a, p.anat_H, p.anat_c], size=(4 * N, 3))
            pts = np.vstack([pts, c[self._inside(c)]])
        pts = pts[:N]
        hilum = np.array([-0.35 * p.anat_a, 0.62 * p.anat_H, 0.0])
        top = hilum + np.array([0.0, 10.0, 0.0])
        parent, start, end, c1, c2, leaf_pt = [-1], [top], [hilum], [-1], [-1], [-1]
        stack = [(0, np.arange(N))]
        while stack:
            k, idx = stack.pop()
            S = pts[idx]
            cen = S.mean(0)
            if len(idx) == 2:
                halves = [idx[:1], idx[1:]]
            else:
                w, V = np.linalg.eigh(np.cov((S - cen).T))
                ax = V[:, -1]
                s = (S - cen) @ ax
                left = s < 0
                if left.all() or (~left).all():
                    left = s < np.median(s)
                halves = [idx[left], idx[~left]]
            kids = []
            for h in halves:
                j = len(parent)
                parent.append(k); start.append(end[k]); c1.append(-1); c2.append(-1)
                if len(h) == 1:
                    end.append(pts[h[0]]); leaf_pt.append(h[0])
                else:
                    tgt = pts[h].mean(0)
                    end.append(end[k] + 0.4 * (tgt - end[k])); leaf_pt.append(-1)
                    stack.append((j, h))
                kids.append(j)
            c1[k], c2[k] = kids
        self.n = len(parent)
        self.N = N
        self.parent = np.array(parent); self.c1 = np.array(c1); self.c2 = np.array(c2)
        self.pos_start = np.array(start); self.pos_end = np.array(end)
        self.L_geo = np.maximum(np.linalg.norm(self.pos_end - self.pos_start, axis=1), 0.1)   # cm
        self.leaf_pt = np.array(leaf_pt)
        self.unit_xyz = None

    def _grid(self):
        p = self.p
        P = self.pos_end[self.leaf]
        self.unit_xyz = P
        self.y = (P[:, 1] - P[:, 1].min()) / np.ptp(P[:, 1])
        self.x = (P[:, 0] - P[:, 0].min()) / np.ptp(P[:, 0])
        tree = cKDTree(P)
        d, _ = tree.query(P, k=2)
        r = 1.6 * np.median(d[:, 1])
        self.nbrs = [[j for j in tree.query_ball_point(P[i], r) if j != i] or [int(tree.query(P[i], k=2)[1][1])]
                     for i in range(self.N)]
        s = np.sqrt((P[:, 0] / p.anat_a) ** 2 + (P[:, 2] / p.anat_c) ** 2 + (P[:, 1] / p.anat_H) ** 2)
        self.pleural_dist = np.round((1 - np.clip(s, 0, 1)) * 10).astype(int)
        self.ix = np.round(self.x * 31).astype(int); self.iy = np.round(self.y * 31).astype(int)
        self.nx = self.ny = 32

    def _vessel_lengths(self):
        L = self.L_geo / 100.0
        L[0] = self.p.L_pa
        return L

    def _build_airways(self):
        p = self.p
        self.r_aw = p.r_tr * (self.nleaf / self.N) ** p.aw_r_exp
        self.L_aw = self.L_geo / 100.0
        self.L_aw[0] = p.L_tr
        self.R_aw = p.k_inertial * pois(p.mu_a, self.L_aw, self.r_aw) * 1e-3 / CMH2O
        self.V_cond = float((np.pi * self.r_aw ** 2 * self.L_aw).sum() * 1e3)
        g = 1 + p.grav * (0.5 - self.y); g /= g.mean()
        self.C_h = p.C_lung / self.N * g
        v = 1 + p.grav_V * (self.y - 0.5); v /= v.mean()
        self.V_h = p.FRC / self.N * v
        self.R_unit = self.N * (p.R_periph + p.R_tissue)

    def pvv(self, vs, cut=None):
        """Large = generations 0-4, small = the terminal quarter of the path (depth >= 75% of max)."""
        cut = cut or (4, int(0.75 * self.Dmax))
        return super().pvv(vs, cut)
