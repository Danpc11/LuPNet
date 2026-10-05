#!/usr/bin/env python3
"""Figures 1-4 and Table 1 of the manuscript on pulmonary hypertension in IPF.

    python3 src/make_figures_erj.py --cohorts results/ph_cohorts --out results/figures
Figure 1: the shared target law against human pulmonary morphometry. Figure 2: the calibrated virtual cohort.
Figure 3: advanced disease (out of sample). Figure 4: the vascular phenotype. Table 1: calibration targets.
"""
import argparse, glob, os, sys
import numpy as np, pandas as pd
import logging
import matplotlib
matplotlib.use("Agg")
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.mixture import GaussianMixture
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lupnet.lung import womersley_G

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ap = argparse.ArgumentParser(); ap.add_argument("--cohorts", default=os.path.join(ROOT, "results", "ph_cohorts"))
ap.add_argument("--out", default=os.path.join(ROOT, "results", "figures")); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
MM = 1 / 25.4
plt.rcParams.update({"font.family": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 600,
                     "lines.linewidth": 1.0, "pdf.fonttype": 42})
BLUE, RED, GREEN, GREY, LGREY, DARK = "#2166AC", "#B2182B", "#1B7837", "#8C8C8C", "#D9D9D9", "#252525"
rng0 = np.random.default_rng(0)


def tag(ax, s, x=-0.22, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    for ext in ("pdf", "png", "tif"):
        kw = dict(pil_kwargs={"compression": "tiff_lzw"}) if ext == "tif" else {}
        fig.savefig(os.path.join(a.out, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.04, **kw)
    plt.close(fig)


# ------------------------------------------------------------------ data
def cohort(tag_):
    from cohort_io import cases
    rows = []
    for _, d in cases(a.cohorts, tag_):
        r0, x = d.iloc[0], d.iloc[-1]
        rows.append(dict(mPAP=x.mPAP, FVC=100 * x.FVC / r0.FVC, DLCO=100 * x.DLCO / r0.DLCO, A=float(x.vasc_A)))
    R = pd.DataFrame(rows); R["vasc"] = R.A > 0; R["PH"] = R.mPAP >= 25
    return R


CAL, ADV = cohort("big2"), cohort("adv2")
MILD = CAL[(CAL.FVC >= 50) & (CAL.FVC <= 90)].reset_index(drop=True)


def noise_avg_r(D, col, n=1000, seed=0):
    """Pearson r between mPAP and FVC or DLCO, averaged over n realisations of measurement noise."""
    rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        mp = D.mPAP + rng.normal(0, 2.5, len(D)); x = D[col] + rng.normal(0, 5, len(D))
        out.append(pearsonr(mp, x)[0])
    return float(np.mean(out))


def noisy(D, rng):
    return (D.mPAP + rng.normal(0, 2.5, len(D)), D.FVC + rng.normal(0, 5, len(D)), (D.DLCO + rng.normal(0, 5, len(D))).clip(lower=3))


def network_schematic(ax):
    """Panel a of Figure 1: the networks modelled by LuPNet."""
    from matplotlib.patches import FancyBboxPatch, Polygon
    ax.set_xlim(0, 10); ax.set_ylim(0.15, 3.75); ax.axis("off")
    yl = np.linspace(0.35, 2.85, 8)                      # eight units drawn
    xr, dx = 0.55, 0.62
    def branch(x0, y0, level, ys, sign):
        if level == 3:
            return
        half = len(ys) // 2
        for sub in (ys[:half], ys[half:]):
            yc = float(np.mean(sub)); x1 = x0 + sign * dx
            ax.plot([x0, x0 + sign * dx * 0.35, x1], [y0, yc, yc], color=BLUE if sign > 0 else "#4D4D9F", lw=2.4 - 0.55 * level,
                    solid_capstyle="round", zorder=2)
            branch(x1, yc, level + 1, sub, sign)
    ax.plot([0.12, xr], [1.6, 1.6], color=BLUE, lw=3.0, solid_capstyle="round")
    branch(xr, 1.6, 0, yl, +1)
    xa = xr + 3 * dx                                      # arterial leaves
    states = ["h", "h", "f", "h", "c", "f", "h", "h"]
    col = {"h": "#9ECAE1", "f": "#FC9272", "c": "#BDBDBD"}
    xu0, xu1 = xa + 0.55, xa + 1.45
    for k, y in enumerate(yl):
        tri = Polygon([[xa, y], [xa + 0.42, y + 0.11], [xa + 0.42, y - 0.11]], closed=True, fc="white", ec=BLUE, lw=0.6, zorder=3)
        ax.add_patch(tri)
        ax.plot([xa + 0.42, xu0], [y, y], color=BLUE, lw=0.6, zorder=2)
        if k in (2, 5):                                   # vasculopathy: series lesion
            ax.plot([xa + 0.47, xa + 0.53], [y, y], color=RED, lw=3.2, solid_capstyle="butt", zorder=4)
        ax.add_patch(FancyBboxPatch((xu0, y - 0.11), xu1 - xu0, 0.22, boxstyle="round,pad=0.0,rounding_size=0.08",
                                    fc=col[states[k]], ec=DARK, lw=0.5, zorder=3))
        ax.plot([xu1, xu1 + 0.25], [y, y], color="#4D4D9F", lw=0.6, zorder=2)
    xv = xu1 + 0.25
    def merge(x0, ys, level):
        if level == 3:
            return
        half = len(ys) // 2
        for sub in (ys[:half], ys[half:]):
            merge(x0, sub, level + 1)
    # venous tree as the mirror image, converging to the left atrium
    def vbranch(x0, y0, level, ys):
        if level == 3:
            return
        half = len(ys) // 2
        for sub in (ys[:half], ys[half:]):
            yc = float(np.mean(sub)); x1 = x0 - dx
            ax.plot([x0, x0 - dx * 0.35, x1], [y0, yc, yc], color="#4D4D9F", lw=2.4 - 0.55 * level, solid_capstyle="round", zorder=2)
            vbranch(x1, yc, level + 1, sub)
    xla = xv + 3 * dx
    vbranch(xla, 1.6, 0, yl)
    ax.plot([xla, xla + 0.35], [1.6, 1.6], color="#4D4D9F", lw=3.0, solid_capstyle="round")
    # labels
    ax.text(0.12, 1.78, "Main pulmonary\nartery", fontsize=6, va="bottom")
    ax.text(xla + 0.05, 1.78, "Left atrium\n(8 mmHg)", fontsize=6, va="bottom")
    ax.text(xr + 1.5 * dx, 3.62, "Arterial tree\n10 explicit generations", fontsize=6, ha="center", va="top", color=BLUE)
    ax.text(xa + 0.21, 3.62, "Sub-tree to\n10 µm radius", fontsize=6, ha="center", va="top", color=BLUE)
    ax.text(0.5 * (xu0 + xu1), 3.62, "Gas-exchange\nunits", fontsize=6, ha="center", va="top")
    ax.text(xv + 1.5 * dx, 3.62, "Venous tree", fontsize=6, ha="center", va="top", color="#4D4D9F")
    # key on the right
    x0k = 7.45; yk = 3.3; dy = 0.36
    ax.plot([x0k, x0k + 0.25], [yk, yk], color=BLUE, lw=2.0); ax.text(x0k + 0.35, yk, "Shared law: τ* = τ₀ r$^{b-1}$ ℓ$^{(b-1)/2}$", va="center", fontsize=6)
    ax.text(x0k + 0.35, yk - 0.75 * dy, "sensed shear S = τ̄ [1 + w₁φ₁G(α)]", va="center", fontsize=6, color=GREY)
    items = [("#9ECAE1", "Healthy unit: ventilation, capillary sheet,"), (None, "membrane and capillary diffusion (DLCO)"),
             ("#FC9272", "Fibrotic unit: low compliance and capillary bed"), ("#BDBDBD", "Honeycomb unit")]
    y = yk - 2.0 * dy
    for c, t in items:
        if c:
            ax.add_patch(FancyBboxPatch((x0k, y - 0.07), 0.25, 0.14, boxstyle="round,pad=0.0,rounding_size=0.05", fc=c, ec=DARK, lw=0.5))
        ax.text(x0k + 0.35, y, t, va="center", fontsize=6); y -= dy * (0.9 if c is None else 1.0)
    ax.plot([x0k + 0.08, x0k + 0.17], [y, y], color=RED, lw=3.2, solid_capstyle="butt")
    ax.text(x0k + 0.35, y, "Vasculopathy: series lesion and capillary loss", va="center", fontsize=6); y -= dy
    ax.text(x0k + 0.35, y, "(half of the simulated cases; equation 3)", va="center", fontsize=6, color=GREY)
    tag(ax, "a", x=-0.05, y=0.98)


# ------------------------------------------------------------------ Figure 1: shared law vs morphometry
Q = 5.0 / 60 * 1e-3; RHO = 1060.0; HR = 70.0; W1 = 2.865; OM = 2 * np.pi * HR / 60


def tree(fname):
    H = pd.read_csv(os.path.join(ROOT, "data", fname), sep="\t")
    return H, Q / H.N.values, H.L_cm.values / 100, H.mu_cP.values * 1e-3, H.D_cm.values / 100


def radius(Qv, L, mu, b, tau0, phi):
    f = lambda x: np.log(4 * mu * Qv / (np.pi * np.exp(3 * x)) * (1 + W1 * phi * womersley_G(np.exp(x) * np.sqrt(OM * RHO / mu))[0])) - \
        np.log(tau0 * np.exp(x) ** (b - 1) * L ** ((b - 1) / 2))
    return float(np.exp(brentq(f, np.log(1e-9), np.log(5.0))))


def predict(T, b, tau0, phi):
    return np.array([2 * radius(q, l, m, b, tau0, phi) for q, l, m in zip(T[1], T[2], T[3])])


def err(T, D): return float(np.sqrt(np.mean((np.log10(D) - np.log10(T[4])) ** 2)))


def fit(T, b, phi):
    r = minimize_scalar(lambda lt: err(T, predict(T, b, 10 ** lt, phi)) ** 2, bounds=(-3, 2), method="bounded", options=dict(xatol=1e-4))
    return float(np.sqrt(r.fun)), 10 ** r.x


ART, VEN = tree("huang1996_pulmonary_arteries.tsv"), tree("huang1996_pulmonary_veins.tsv")
PHI_A, PHI_V = 0.8, 0.4
bgrid = np.round(np.arange(0.55, 1.101, 0.025), 3)
prof_p = np.array([fit(ART, b, PHI_A)[0] for b in bgrid]); prof_m = np.array([fit(ART, b, 0.0)[0] for b in bgrid])
b_best = float(bgrid[int(prof_p.argmin())])
_, t_best = fit(ART, b_best, PHI_A); _, t_675 = fit(ART, 0.675, PHI_A)
D_best, D_675 = predict(ART, b_best, t_best, PHI_A), predict(ART, 0.675, t_675, PHI_A)
D_ven = predict(VEN, b_best, t_best, PHI_V)
e_art, e_ven = 10 ** err(ART, D_best), 10 ** err(VEN, D_ven)

fig = plt.figure(figsize=(180 * MM, 178 * MM))
gs = fig.add_gridspec(3, 3, hspace=0.62, wspace=0.48, left=0.07, right=0.98, top=0.95, bottom=0.07, height_ratios=[0.95, 1, 1])
network_schematic(fig.add_subplot(gs[0, :]))
ax = fig.add_subplot(gs[1, 0])
o = ART[0].order.values
ax.semilogy(o, ART[4] * 1e3, "o", ms=3.2, mfc="white", mec=DARK, mew=0.8, label="Measured (Huang 1996)")
ax.semilogy(o, D_best * 1e3, "-", color=BLUE, label=f"Shared law, b = {b_best:.3g}")
ax.semilogy(o, D_675 * 1e3, "--", color=RED, label="Systemic exponent, b = 0.675")
ax.set(xlabel="Strahler order", ylabel="Diameter (mm)", xticks=[1, 4, 8, 12, 16], ylim=(0.01, 3000))
ax.legend(loc="upper left", handlelength=1.6)
tag(ax, "b")
ax = fig.add_subplot(gs[1, 1])
lim = [0.01, 50]
ax.loglog(lim, lim, color=LGREY, lw=0.8)
ax.loglog(ART[4] * 1e3, D_best * 1e3, "o", ms=3.2, color=BLUE, label=f"Arteries (fitted), ×{e_art:.2f}")
ax.loglog(VEN[4] * 1e3, D_ven * 1e3, "s", ms=3.0, mfc="white", mec=GREEN, mew=0.9, label=f"Veins (not refitted), ×{e_ven:.2f}")
ax.set(xlim=lim, ylim=lim, xlabel="Measured diameter (mm)", ylabel="Predicted diameter (mm)")
ax.legend(loc="upper left", handlelength=1.2)
tag(ax, "c")
ax = fig.add_subplot(gs[1, 2])
ax.plot(bgrid, 10 ** prof_p, color=BLUE, label="Pulsatile sensing")
ax.plot(bgrid, 10 ** prof_m, color=GREY, ls="--", label="Mean shear only")
ax.axvline(0.675, color=RED, lw=0.8, ls=":"); ax.text(0.668, 1.012, "systemic\nb = 0.675", color=RED, fontsize=6, va="bottom", ha="right")
ax.axvline(1.0, color=DARK, lw=0.6, ls=":"); ax.text(1.008, 1.012, "Murray", fontsize=6, va="bottom")
ax.set(xlabel="Maintenance-cost exponent b", ylabel="Diameter error (geometric factor)", ylim=(1.0, 1.45))
ax.legend(loc="upper center", bbox_to_anchor=(0.66, 1.02), handlelength=1.6)
tag(ax, "d")
# d: pressure partition of the healthy lung (shared law)
from lupnet.lung import Lung
from lupnet.params import Params
Lh = Lung(Params(law="shared", b=0.925)); vs = Lh.base_state["vs"]
Ll = Lung(Params())
ax = fig.add_subplot(gs[2, 0])
ax.axhspan(14.0 - 3.3, 14.0 + 3.3, color=LGREY, lw=0)
ax.text(1.45, 14.0 + 3.3 - 0.3, "Healthy adults\n14.0 ± 3.3 mmHg", fontsize=6, va="top", ha="right", color=DARK)
ax.plot([0], [Ll.base_state["vs"]["PPA"]], "o", ms=4.5, mfc="white", mec=GREY, mew=1.0)
ax.plot([1], [vs["PPA"]], "o", ms=4.5, color=BLUE)
ax.set(xlim=(-0.5, 1.5), ylim=(8, 20), xticks=[0, 1], ylabel="mPAP (mmHg)")
ax.set_xticklabels(["Per-vessel law\n(calibrated)", "Shared law\n(predicted)"])
tag(ax, "e")
# e: generation diameters of the integrated lung against Huang
ax = fig.add_subplot(gs[2, 1])
d = Lh.depth; ra = vs["r_a"]; LPM = 1e-3 / 60
fl = np.log(Q / ART[0].N.values); Dh = np.log(ART[4])
g = np.arange(0, d.max() + 1); Dm = np.array([2 * np.median(ra[d == k]) for k in g])
Dhu = np.exp(np.interp(np.log(Q / 2 ** g), fl, Dh))
ax.semilogy(g, Dhu * 1e3, "o", ms=3.2, mfc="white", mec=DARK, mew=0.8, label="Huang 1996, matched by flow")
ax.semilogy(g, Dm * 1e3, "-", color=BLUE, label="LuPNet tree, shared law")
ax.set(xlabel="Generation of the explicit tree", ylabel="Diameter (mm)")
ax.legend(loc="lower left", handlelength=1.6)
tag(ax, "f")
# f: where the pressure drops
ax = fig.add_subplot(gs[2, 2])
mmHg = 133.322
dp_ord = 8 * ART[3] * ART[2] * ART[1] / (np.pi * (D_best / 2) ** 4) / mmHg
ax.bar(o, dp_ord, color=BLUE, width=0.7)
ax.set(xlabel="Strahler order", ylabel="Pressure drop per order (mmHg)", xticks=[1, 4, 8, 12, 16])
ax.text(0.98, 0.95, f"arterial total {dp_ord.sum():.2f} mmHg", transform=ax.transAxes, ha="right", va="top", fontsize=6)
tag(ax, "g")
save(fig, "Figure1")

# ------------------------------------------------------------------ Figure 2: calibrated cohort
fig = plt.figure(figsize=(180 * MM, 118 * MM))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.48, left=0.07, right=0.98, top=0.92, bottom=0.10)
bins = [("<20", lambda m: m < 20, 51), ("20–25", lambda m: (m >= 20) & (m < 25), 30), ("25–35", lambda m: (m >= 25) & (m < 35), 15), ("≥35", lambda m: m >= 35, 4)]
B = np.array([[100 * f(MILD.sample(len(MILD), replace=True, random_state=int(rng0.integers(1e9))).mPAP).mean() for _, f, _ in bins] for _ in range(2000)])
pt = np.array([100 * f(MILD.mPAP).mean() for _, f, _ in bins]); lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
ax = fig.add_subplot(gs[0, 0]); x = np.arange(4)
ax.bar(x - 0.19, pt, 0.36, color=BLUE, label="Model (n = %d)" % len(MILD))
ax.errorbar(x - 0.19, pt, yerr=[pt - lo, hi - pt], fmt="none", ecolor=DARK, elinewidth=0.7, capsize=1.5)
ax.bar(x + 0.19, [t for *_, t in bins], 0.36, color=LGREY, edgecolor=DARK, lw=0.5, label="ARTEMIS-IPF (n = 488)")
ax.set(xticks=x, ylabel="Share (%)", xlabel="Mean pulmonary artery pressure (mmHg)", ylim=(0, 72))
ax.set_xticklabels([b for b, *_ in bins]); ax.legend(loc="upper right")
tag(ax, "a")
for k, (col, lab, letter) in enumerate((("FVC", "FVC (% predicted)", "b"), ("DLCO", "DLCO (% predicted)", "c"))):
    ax = fig.add_subplot(gs[0, 1 + k])
    for v, c, lb in ((False, GREY, "No vasculopathy"), (True, RED, "Vasculopathy")):
        m = MILD.vasc == v
        ax.plot(MILD[col][m], MILD.mPAP[m], "o", ms=2.6, color=c, alpha=0.75, mew=0, label=lb)
    ax.axhline(25, color=DARK, lw=0.5, ls=":")
    rmean = noise_avg_r(MILD, col)
    ref = "0" if col == "FVC" else "−0.30"
    ax.text(0.97, 0.95, f"model r = {rmean:.2f}".replace("-", "−") + f"\npatients r ≈ {ref}", transform=ax.transAxes, ha="right", va="top", fontsize=6)
    ax.set(xlabel=lab, ylabel="mPAP (mmHg)", ylim=(12, 50))
    if k == 0: ax.legend(loc="upper left", handlelength=0.8, markerscale=1.3)
    tag(ax, letter)
ax = fig.add_subplot(gs[1, 0])
parts = [MILD.mPAP[~MILD.vasc], MILD.mPAP[MILD.vasc]]
vp = ax.violinplot(parts, positions=[0, 1], widths=0.7, showextrema=False)
for body, c in zip(vp["bodies"], (GREY, RED)): body.set_facecolor(c); body.set_alpha(0.35); body.set_edgecolor("none")
for i, (p_, c) in enumerate(zip(parts, (GREY, RED))):
    ax.plot(np.full(len(p_), i) + rng0.uniform(-0.12, 0.12, len(p_)), p_, "o", ms=1.8, color=c, mew=0, alpha=0.8)
ax.axhline(25, color=DARK, lw=0.5, ls=":")
ax.set(xticks=[0, 1], xlim=(-0.6, 1.6), ylabel="mPAP (mmHg)", ylim=(12, 50))
ax.set_xticklabels([f"No vasculopathy\n(n = {len(parts[0])})", f"Vasculopathy\n(n = {len(parts[1])})"])
tag(ax, "d")
# e: Table-like dot plot of the seven targets relative to their CI
from ph_stats import bootstrap as ph_bootstrap, point as ph_point, prevalence_ci
CI_MILD = ph_bootstrap(MILD); PT_MILD = ph_point(MILD)
S = {"r_DLCO": (CI_MILD["rDLCO"], PT_MILD["rDLCO"]), "r_FVC": (CI_MILD["rFVC"], PT_MILD["rFVC"]), "ratio": (CI_MILD["ratio"], PT_MILD["ratio"])}
sub = gs[1, 1:].subgridspec(1, 2, width_ratios=[1.0, 0.75], wspace=0.08)
axL, axR = fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])
items = [("r(mPAP, DLCO)", "r_DLCO", -0.30), ("r(mPAP, FVC)", "r_FVC", 0.0), ("DLCO ratio,\nPH / no PH", "ratio", 0.73)]
for i, (lab, k, tgt) in enumerate(items):
    ax_ = axL if k != "ratio" else axR
    (l_, h_), m_ = S[k]
    ax_.plot([l_, h_], [i, i], color=BLUE, lw=1.6, solid_capstyle="butt"); ax_.plot(m_, i, "o", ms=4, color=BLUE)
    ax_.plot(tgt, i, "D", ms=4, mfc="white", mec=DARK, mew=0.9)
for ax_ in (axL, axR):
    ax_.set_ylim(2.5, -0.5); ax_.set_yticks(range(3))
axL.set_yticklabels([l for l, *_ in items]); axL.set_xlim(-0.6, 0.15); axL.set_xticks([-0.6, -0.4, -0.2, 0.0])
axL.set_xlabel("Correlation coefficient")
axR.set_xlim(0.4, 1.0); axR.set_xticks([0.4, 0.6, 0.8, 1.0]); axR.set_xlabel("Ratio")
axR.tick_params(axis="y", length=0); axR.set_yticklabels([]); axR.spines["left"].set_visible(False)
axL.axhline(1.5, color=LGREY, lw=0.5); axR.axhline(1.5, color=LGREY, lw=0.5)
axR.plot([], [], color=BLUE, lw=1.6, marker="o", ms=4, label="Model, estimate\nand 95% CI")
axR.plot([], [], "D", ms=4, mfc="white", mec=DARK, mew=0.9, label="Patients")
axR.legend(loc="upper right", handlelength=1.4)
tag(axL, "e", x=-0.42)
save(fig, "Figure2")

# ------------------------------------------------------------------ Figure 3: advanced disease
POOL = pd.concat([CAL, ADV], ignore_index=True)
fig = plt.figure(figsize=(180 * MM, 60 * MM))
gs = fig.add_gridspec(1, 3, wspace=0.48, left=0.07, right=0.98, top=0.86, bottom=0.2)
bands = [(40, 45), (45, 50), (50, 55), (55, 60), (60, 70), (70, 85)]
def prev(D, f):
    v = f(D); B_ = [f(D.sample(len(D), replace=True, random_state=int(rng0.integers(1e9)))) for _ in range(1000)]
    return v, *np.percentile(B_, [2.5, 97.5])
for k, (f, lab, ref, letter) in enumerate(((lambda D: 100 * (D.mPAP >= 25).mean(), "mPAP ≥ 25 mmHg (%)", 46, "a"),
                                           (lambda D: 100 * (D.mPAP > 40).mean(), "mPAP > 40 mmHg (%)", 9, "b"))):
    ax = fig.add_subplot(gs[0, k]); xs = []
    for j, (b0, b1) in enumerate(bands):
        D = POOL[(POOL.FVC >= b0) & (POOL.FVC < b1)]
        if len(D) < 8: continue
        v, l_, h_ = prev(D, f); c = 0.5 * (b0 + b1)
        ax.errorbar(c, v, yerr=[[v - l_], [h_ - v]], fmt="o", ms=3.5, color=BLUE, ecolor=BLUE, elinewidth=0.8, capsize=1.5)
        if k == 0: ax.text(c, -6, f"{len(D)}", ha="center", fontsize=5.5, color=GREY)
    ax.errorbar(49.6, ref, xerr=15.7, fmt="D", ms=4, mfc="white", mec=DARK, ecolor=DARK, elinewidth=0.6, capsize=1.5)
    ax.annotate("Transplant candidates\n(Shorr 2007, n = 2,525)", xy=(49.6, ref), xytext=(60, 72 if k == 0 else 22), fontsize=6,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=DARK))
    ax.set(xlabel="FVC (% predicted)", ylabel=lab, xlim=(32, 88), ylim=(-9, 100) if k == 0 else (-1.5, 30))
    tag(ax, letter)
ax = fig.add_subplot(gs[0, 2])
A5 = POOL[POOL.FVC < 50]
for i, (m, c, lab) in enumerate(((~A5.PH, GREY, "No PH"), (A5.PH, RED, "PH"))):
    v = A5.FVC[m]; ax.boxplot([v], positions=[i], widths=0.45, showfliers=False, medianprops=dict(color=DARK, lw=0.9),
                               boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.plot(np.full(len(v), i) + rng0.uniform(-0.12, 0.12, len(v)), v, "o", ms=2, color=c, mew=0, alpha=0.8)
ax.plot([0, 1], [51.4, 48.4], "D", ms=4, mfc="white", mec=DARK, mew=0.9, label="Transplant candidates, mean")
ax.set(xticks=[0, 1], xlim=(-0.6, 1.6), ylabel="FVC (% predicted)", ylim=(38, 56))
ax.set_xticklabels([f"No PH (n = {int((~A5.PH).sum())})", f"PH (n = {int(A5.PH.sum())})"]); ax.legend(loc="upper right", handlelength=0.8)
tag(ax, "c")
save(fig, "Figure3")

# ------------------------------------------------------------------ Figure 4: vascular phenotype
fig = plt.figure(figsize=(180 * MM, 60 * MM))
gs = fig.add_gridspec(1, 3, wspace=0.5, left=0.07, right=0.98, top=0.86, bottom=0.2)
nov = CAL[~CAL.vasc]; dl = np.polyfit(nov.FVC, np.log(nov.DLCO), 2)
res = np.log(CAL.DLCO) - np.polyval(dl, CAL.FVC)
gm = GaussianMixture(2, random_state=0).fit(res.values.reshape(-1, 1)); g1 = GaussianMixture(1, random_state=0).fit(res.values.reshape(-1, 1))
dbic = g1.bic(res.values.reshape(-1, 1)) - gm.bic(res.values.reshape(-1, 1))
ax = fig.add_subplot(gs[0, 0])
ax.hist(res, bins=28, color=LGREY, edgecolor="white", lw=0.4, density=True)
xx = np.linspace(res.min() - 0.1, res.max() + 0.1, 300); low = int(np.argmin(gm.means_.ravel()))
for k in range(2):
    w, m, s = gm.weights_[k], gm.means_.ravel()[k], np.sqrt(gm.covariances_.ravel()[k])
    ax.plot(xx, w * np.exp(-0.5 * ((xx - m) / s) ** 2) / (s * np.sqrt(2 * np.pi)), color=RED if k == low else BLUE)
ax.text(0.03, 0.95, f"Two components\nΔBIC = {dbic:.1f}", transform=ax.transAxes, va="top", fontsize=6)
ax.set(xlabel="DLCO deficit for the FVC (ln units)", ylabel="Density")
tag(ax, "a")
ax = fig.add_subplot(gs[0, 1])
sc = ax.scatter(CAL.FVC, CAL.DLCO, c=CAL.mPAP.clip(15, 35), cmap="RdBu_r", s=6, lw=0, vmin=15, vmax=35)
ff = np.linspace(CAL.FVC.min(), CAL.FVC.max(), 100); ax.plot(ff, np.exp(np.polyval(dl, ff)), color=DARK, lw=0.9)
cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03, ticks=[15, 20, 25, 30, 35]); cb.set_label("mPAP (mmHg)", fontsize=6.5); cb.ax.tick_params(labelsize=6, width=0.5, length=2)
cb.outline.set_linewidth(0.5)
ax.set(xlabel="FVC (% predicted)", ylabel="DLCO (% predicted)")
tag(ax, "b")
ax = fig.add_subplot(gs[0, 2])
def mean_auc(score_fn, n=1000, seed=0):
    rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        mp_, fv_, dl_ = noisy(CAL, rng); out.append(roc_auc_score(mp_ >= 25, score_fn(fv_, dl_)))
    return float(np.mean(out))
AUC = {"FVC/DLCO": mean_auc(lambda f, d: f / d), "DLCO": mean_auc(lambda f, d: -d), "FVC": mean_auc(lambda f, d: -f)}
mp, fv, dlc = noisy(CAL, np.random.default_rng(5)); y = mp >= 25
for score, c, lab in ((fv / dlc, BLUE, "FVC/DLCO"), (-dlc, GREEN, "DLCO"), (-fv, GREY, "FVC")):
    fpr, tpr, _ = roc_curve(y, score); ax.plot(fpr, tpr, color=c, label=f"{lab}, AUC {AUC[lab]:.2f}")
ax.plot([0, 1], [0, 1], color=LGREY, lw=0.6)
ax.set(xlabel="1 − specificity", ylabel="Sensitivity", xlim=(0, 1), ylim=(0, 1), aspect="equal")
ax.legend(loc="lower right", handlelength=1.2)
ax.text(0.03, 0.98, "Patients:\nFVC/DLCO 0.69–0.74\nDLCO 0.80", transform=ax.transAxes, va="top", ha="left", fontsize=5.5, color=DARK)
tag(ax, "c")
save(fig, "Figure4")

# ------------------------------------------------------------------ Table 1
T = {k: CI_MILD[k] for k in CI_MILD}
P0 = {"lt20": PT_MILD["lt20"], "b20_25": PT_MILD["b20_25"], "ge25": PT_MILD["ge25"], "ge35": PT_MILD["ge35"],
      "rD": PT_MILD["rDLCO"], "rF": PT_MILD["rFVC"], "ratio": PT_MILD["ratio"]}
T["rD"], T["rF"] = CI_MILD["rDLCO"], CI_MILD["rFVC"]
rows = [("mPAP < 20 mmHg, %", "lt20", "51", "ARTEMIS-IPF"), ("mPAP 20–25 mmHg, %", "b20_25", "30", "ARTEMIS-IPF"), ("mPAP ≥ 25 mmHg, %", "ge25", "19", "ARTEMIS-IPF"),
        ("mPAP ≥ 35 mmHg, %", "ge35", "4", "ARTEMIS-IPF"), ("r(mPAP, DLCO)", "rD", "≈ −0.30", "Zisman 2007"), ("r(mPAP, FVC)", "rF", "≈ 0", "Zisman 2007; ARTEMIS-IPF"),
        ("DLCO, PH / no PH", "ratio", "0.73", "Japanese cohort")]
tab = []
for lab, k, tgt, src in rows:
    l_, h_ = T[k]; fmt = "{:.0f}" if k in ("lt20", "b20_25", "ge25", "ge35") else "{:.2f}"
    tab.append(dict(Quantity=lab, Model=fmt.format(P0[k]), CI95=f"{fmt.format(l_)} to {fmt.format(h_)}", Patients=tgt, Source=src))
pd.DataFrame(tab).to_csv(os.path.join(a.out, "Table1.tsv"), sep="\t", index=False)
print(pd.DataFrame(tab).to_string(index=False))
print("AUC (mean of 1,000 noise realisations):", {k: round(v, 3) for k, v in AUC.items()})
print(f"b_best {b_best}, arterial factor {e_art:.3f}, venous factor {e_ven:.3f}, healthy mPAP {vs['PPA']:.2f}; dBIC {dbic:.1f}; n mild {len(MILD)}, adv<50 {len(A5)}")
