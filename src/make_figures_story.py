#!/usr/bin/env python3
"""Figures of the storyline "pulmonary hypertension in fibrosing lung disease arises from vascular disease in
preserved lung, not from the scar" (Story1-Story6 in results/figures/). Figure 1 is Figure1 of make_figures_erj.py.

  Story2  fibrosis alone does not produce PH, and it multiplies the effect of a lesion
  Story3  where the lesion sits: scar, field, half field, diffuse
  Story4  the diffuse lesion reproduces the disease with the UIP and the NSIP pattern
  Story5  catheterised patients                                     (needs --clinical)
  Story6  clinical consequence: flow unmasks the lesion; estimate without catheterisation  (needs --clinical and
          the tables of `noninvasive_tool.py --clinical`)

    python3 src/make_figures_story.py [--clinical patients.xlsx]
Three example lungs are simulated once (about one minute) and cached in results/figures/story_maps.npz.
"""
import argparse, logging, os, sys, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR); warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr, linregress
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import final_visits, cases
from ph_stats import bootstrap, point, prevalence_ci
import noninvasive_tool as nt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); F = os.path.join(ROOT, "results", "ph_cohorts"); OUT = os.path.join(ROOT, "results", "figures")
MM = 1 / 25.4
plt.rcParams.update({"font.family": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False, "axes.linewidth": 0.6, "savefig.dpi": 600, "pdf.fonttype": 42,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.labelpad": 2.5})

# ---- one palette for every figure (Okabe-Ito, colour-blind safe); a colour always means the same thing
DARK, NEUTRAL, LIGHT = "#222222", "#8C8C8C", "#D9D9D9"      # text and patients' measurements; no lesion; background
LESION = "#D55E00"                                           # a vascular lesion, and the diffuse lesion in preserved lung
MODEL = "#0072B2"                                            # an estimate of the model for a patient
LOC = {"scar": "#8C8C8C", "field": "#E69F00", "mix50": "#CC79A7", "diffuse": "#D55E00"}      # where the lesion sits
PATTERN = {"calibrated": "#8C8C8C", "uip": "#0072B2", "nsip": "#009E73"}                     # histological pattern
TISSUE = ["#F2F2F2", "#B0B0B0", "#4D4D4D"]                   # preserved, fibrotic, honeycomb: tissue is always grey
SEV = matplotlib.colors.LinearSegmentedColormap.from_list("sev", ["#FDE3CF", "#D55E00", "#6B2A00"])   # lesion severity
MECH = {"scar": ("Scar only", LOC["scar"], ("big2_scar",)), "field": ("Spreading", LOC["field"], ("big2_field", "adv2_field")),
        "mix50": ("Half and half", LOC["mix50"], ("big2_mix50", "adv2_mix50")), "diffuse": ("Diffuse", LOC["diffuse"], ("big2", "adv2"))}
PAT = {"calibrated": ("Calibrated", PATTERN["calibrated"], ("big2", "adv2")), "uip": ("UIP", PATTERN["uip"], ("big2_uip", "adv2_uip")),
       "nsip": ("NSIP", PATTERN["nsip"], ("big2_nsip", "adv2_nsip"))}


class Grid:
    """Fixed grid of cells on a canvas of fixed size (mm). Every panel lives in one cell; its letter sits in the
    top-left corner of the cell, so letters of a column are aligned and those of the first column are flush with
    the left edge of the figure. Axes are inset by a padding that holds tick labels, axis labels and legends, so
    panels cannot overlap. Figures are saved at the canvas size, without cropping."""

    def __init__(self, width, height, heights, widths, gap=(5.0, 4.0)):
        self.W, self.H = width, height; self.fig = plt.figure(figsize=(width * MM, height * MM))
        gx, gy = gap; w = np.array(widths, float); h = np.array(heights, float)
        w = w / w.sum() * (width - gx * (len(w) - 1)); h = h / h.sum() * (height - gy * (len(h) - 1))
        self.x0 = np.concatenate([[0], np.cumsum(w[:-1] + gx)]); self.x1 = self.x0 + w
        top = np.concatenate([[0], np.cumsum(h[:-1] + gy)]); self.y1 = height - top; self.y0 = self.y1 - h

    def box(self, r, c, rs=1, cs=1):
        return self.x0[c], self.y0[r + rs - 1], self.x1[c + cs - 1], self.y1[r]

    def cell(self, r, c, letter=None, rs=1, cs=1, pad=(12.5, 10.0, 2.5, 6.0), axis=True):
        """Axes of the cell (r, c); pad = left, bottom, right, top in mm."""
        x0, y0, x1, y1 = self.box(r, c, rs, cs); l, b, rr, t = pad
        ax = self.fig.add_axes([(x0 + l) / self.W, (y0 + b) / self.H, (x1 - x0 - l - rr) / self.W, (y1 - y0 - b - t) / self.H])
        if not axis:
            ax.set_axis_off()
        if letter:
            self.fig.text(x0 / self.W, y1 / self.H, letter, fontsize=9, fontweight="bold", va="top", ha="left")
        return ax

    def save(self, name):
        for ext in ("pdf", "png"):
            self.fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
        plt.close(self.fig); print(name, "written")


def ptxt(p):
    return f"P = {p:.3f}" if p >= 0.001 else "P < 0.001"


def trend(ax, x, y, color=DARK, lw=0.9):
    lr = linregress(x, y); xx = np.linspace(np.min(x), np.max(x), 30); ax.plot(xx, lr.intercept + lr.slope * xx, color=color, lw=lw)
    return pearsonr(x, y)


def stat(ax, r, p, n, loc="right"):
    ax.text(0.97 if loc == "right" else 0.03, 0.97, f"r = {r:.2f}, {ptxt(p)}\nn = {n}".replace("-", "−"), transform=ax.transAxes,
            ha=loc, va="top", fontsize=5.8)


def legend_below(ax, dy=-0.30, ncol=1, **kw):
    return ax.legend(loc="upper center", bbox_to_anchor=(0.5, dy), ncol=ncol, fontsize=5.6, handlelength=1.3, columnspacing=1.0, **kw)


# ---------------------------------------------------------------- diagrams
def lung_diagram(ax, x, y, w, h, scar, lesion, flow, title=None):
    """A lung as parallel channels between the pulmonary artery (left) and the veins (right). `scar`: channel is
    scar; `lesion`: 0-1, lesion on its arteriole; `flow`: relative blood flow, drawn as line width."""
    n = len(scar); ys = y + h * (np.arange(n) + 0.5) / n; yc = y + h / 2; xa, xb, xu0, xu1 = x + 0.16 * w, x + 0.84 * w, x + 0.42 * w, x + 0.70 * w
    tot = 2.6
    ax.plot([x, xa], [yc, yc], color=DARK, lw=tot, solid_capstyle="butt"); ax.plot([xb, x + w], [yc, yc], color=DARK, lw=tot, solid_capstyle="butt")
    for k in range(n):
        lw = 0.35 + 1.5 * flow[k]; c = DARK
        ax.plot([xa, xu0 - 0.12 * w, xu0], [yc, ys[k], ys[k]], color=c, lw=lw, solid_joinstyle="round", solid_capstyle="butt")
        ax.plot([xu1, xu1 + 0.03 * w, xb], [ys[k], ys[k], yc], color=c, lw=lw, solid_joinstyle="round", solid_capstyle="butt")
        ax.add_patch(matplotlib.patches.FancyBboxPatch((xu0, ys[k] - 0.32 * h / n), xu1 - xu0, 0.64 * h / n, boxstyle="round,pad=0,rounding_size=0.8",
                                                       fc=TISSUE[2] if scar[k] else TISSUE[0], ec=DARK, lw=0.5, zorder=3))
        if lesion[k] > 0:
            ax.add_patch(matplotlib.patches.Rectangle((xu0 - 0.085 * w, ys[k] - 0.26 * h / n), 0.05 * w, 0.52 * h / n, fc=SEV(0.5 * lesion[k]), ec="none", zorder=4))
    if title:
        ax.text(x + w / 2, y - 0.06 * h, title, ha="center", va="top", fontsize=5.8)


def box_diagram(ax, x, y, w, h, text, fc="white", ec=DARK, fs=5.6, bold=False):
    ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc, ec=ec, lw=0.6))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", linespacing=1.25)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=0.7, mutation_scale=6, shrinkA=0, shrinkB=0))


def visits(tags, with_lost=True):
    """Final viable visit of every case of the cohorts, with the lung lost."""
    out = []
    for t in tags:
        D = final_visits(F, t); D["tag"] = t
        if with_lost:
            lost = {i: 1 - d.healthy.iloc[-1] for i, d in cases(F, t)}; D["lost"] = D.case.map(lost)
        out.append(D)
    return pd.concat(out, ignore_index=True)


def maps():
    """Three example lungs at about 40% of lung lost: calibrated (with the lesion weights), UIP and NSIP."""
    f = os.path.join(OUT, "story_maps.npz")
    if os.path.exists(f):
        return dict(np.load(f))
    from ph_cohort import make_params
    from lupnet.sim import run, lesion_weights
    M = {}
    for pat in ("calibrated", "uip", "nsip"):
        p = make_params(20007, 8, 0.60, 0.0, 0.0, 1.0, 2.0, 0.05, pattern=pat); _, _, L = run(p, do_exercise=False)
        M[f"{pat}_ix"], M[f"{pat}_iy"], M[f"{pat}_state"] = L.ix, L.iy, L.state
        if pat == "calibrated":
            for w in ("scar", "field", "mixed", "diffuse"):
                M[f"w_{w}"] = lesion_weights(L, w, 0.25)
    np.savez(f, **M); return M


def grid(M, pat, values):
    ix, iy = M[f"{pat}_ix"].astype(int), M[f"{pat}_iy"].astype(int); G = np.full((ix.max() + 1, iy.max() + 1), np.nan); G[ix, iy] = values
    return G.T


# ---------------------------------------------------------------- Story2: fibrosis alone, and the partial dependence
def story2(clinical):
    T = nt.trajectory_points(tags=("big2", "adv2", "hom")); Fs, r2 = nt.surface(T); H = visits(("hom",)); L = np.linspace(0.05, 0.80, 40)
    g = Grid(180, 118, heights=[1, 1], widths=[1, 1, 1])
    # A: two ways to raise pressure
    ax = g.cell(0, 0, "A", pad=(2, 3, 2, 7), axis=False); ax.set(xlim=(0, 100), ylim=(0, 100))
    for y0, sc, le, fl, title in ((72, [0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], "Healthy lung"),
                                  (37, [0, 1, 1, 0], [0, 0, 0, 0], [1.7, 0.15, 0.15, 1.7], "Fibrosis: units lost, flow moves to the rest"),
                                  (2, [0, 0, 0, 0], [1, 1, 1, 1], [1, 1, 1, 1], "Vascular lesion: arterioles narrowed")):
        lung_diagram(ax, 6, y0, 88, 20, sc, le, fl); ax.text(50, y0 + 24.5, title, ha="center", va="bottom", fontsize=5.8)
    # B: no lesion, pressure against lung lost
    ax = g.cell(0, 1, "B")
    for i, d in cases(F, "hom"):
        ax.plot(100 * (1 - d.healthy), d.mPAP, color=NEUTRAL, lw=0.3, alpha=0.5)
    ax.plot(100 * L, Fs.predict(np.c_[L, np.zeros(len(L))]), color=DARK, lw=1.2); ax.axhline(25, color=DARK, lw=0.5, ls=":")
    ax.text(0.03, 0.97, f"No lesion, n = {len(H)}", transform=ax.transAxes, va="top", fontsize=5.8)
    ax.set(xlabel="Lung lost (%)", ylabel="mPAP (mmHg)", ylim=(12, 40), xlim=(0, 85))
    # C: no lesion, pressure against FVC
    ax = g.cell(0, 2, "C"); ax.plot(H.FVC, H.mPAP, "o", ms=2.6, color=NEUTRAL, mew=0); r, p = trend(ax, H.FVC, H.mPAP); stat(ax, r, p, len(H))
    ax.axhline(25, color=DARK, lw=0.5, ls=":"); ax.set(xlabel="FVC (% of baseline)", ylabel="mPAP (mmHg)", ylim=(12, 40))
    # D: patients
    ax = g.cell(1, 0, "D")
    if clinical is not None:
        V = clinical.dropna(subset=["FVC", "mPAP"]); ax.plot(V.FVC, V.mPAP, "o", ms=2.8, mfc="white", mec=DARK, mew=0.7); r, p = trend(ax, V.FVC, V.mPAP); stat(ax, r, p, len(V))
        ax.axhline(25, color=DARK, lw=0.5, ls=":"); ax.set(xlabel="FVC (% predicted)", ylabel="mPAP, catheterisation (mmHg)", ylim=(5, 82))
    else:
        ax.set_axis_off()
    # E: the surface
    ax = g.cell(1, 1, "E", pad=(12.5, 10, 9, 6))
    for a in (0, 1, 2, 4, 8):
        y = Fs.predict(np.c_[L, np.full(len(L), a)]); c = DARK if a == 0 else SEV(0.15 + 0.1 * a); ax.plot(100 * L, y, color=c, lw=1.1)
        ax.text(82, y[-1] + (-1.0 if a == 0 else 0.8 if a == 1 else 0), f"{a}", fontsize=5.8, va="center", color=c)
    ax.text(84, 51.5, "Lesion", fontsize=5.8, ha="center"); ax.axhline(25, color=DARK, lw=0.5, ls=":")
    ax.set(xlabel="Lung lost (%)", ylabel="mPAP (mmHg)", ylim=(12, 54), xlim=(0, 80)); ax.set_xticks([0, 20, 40, 60, 80])
    # F: what a lesion adds
    ax = g.cell(1, 2, "F"); lost = (0.05, 0.30, 0.60); w = 0.3
    for j, a in enumerate((2, 4)):
        d = [float(Fs.predict(np.array([[l, a]]))[0] - Fs.predict(np.array([[l, 0.0]]))[0]) for l in lost]
        ax.bar(np.arange(3) + (j - 0.5) * w * 1.08, d, w, color=SEV(0.15 + 0.1 * a), label=f"Lesion {a}")
    ax.set_xticks(range(3)); ax.set_xticklabels(["5", "30", "60"]); ax.set(xlabel="Lung lost (%)", ylabel="mPAP added by the lesion (mmHg)", ylim=(0, 14))
    ax.legend(loc="upper left", fontsize=5.8, handlelength=1)
    g.save("Story2"); print(f"  surface R2 {r2:.3f}; no lesion: {Fs.predict(np.array([[0.05, 0]]))[0]:.1f} mmHg at 5% lost, {Fs.predict(np.array([[0.60, 0]]))[0]:.1f} at 60%")


# ---------------------------------------------------------------- Story3: where the lesion sits
def story3(clinical):
    M = maps(); g = Grid(180, 172, heights=[0.60, 0.66, 1.0], widths=[1, 1, 1, 1]); scar = M["calibrated_state"] > 0
    # A: the four locations as diagrams. Channels: preserved far, preserved next to scar, scar, scar, next to scar, far
    S = [0, 0, 1, 1, 0, 0]; lo = 0.15
    D = {"scar": ([0, 0, 1, 1, 0, 0], [1, 1, lo, lo, 1, 1]), "field": ([0, 1, 1, 1, 1, 0], [1.7, 0.45, lo, lo, 0.45, 1.7]),
         "mix50": ([0.5, 1, 1, 1, 1, 0.5], [1.1, 0.75, lo, lo, 0.75, 1.1]), "diffuse": ([1, 1, 1, 1, 1, 1], [1, 1, lo, lo, 1, 1])}
    names = {"scar": "Scar only", "field": "Spreading from the scar", "mix50": "Half spreading, half diffuse", "diffuse": "Diffuse, in preserved lung"}
    for j, k in enumerate(MECH):
        ax = g.cell(0, j, "A" if j == 0 else None, pad=(3, 3, 3, 9), axis=False); ax.set(xlim=(0, 100), ylim=(0, 100))
        lung_diagram(ax, 2, 4, 96, 92, S, D[k][0], D[k][1]); g.fig.text((g.x0[j] + g.x1[j]) / 2 / g.W, (g.y1[0] - 2.2) / g.H, names[k], ha="center", va="top", fontsize=6.4, color=MECH[k][1], fontweight="bold")
    # B: the same four locations in one simulated lung
    cmap = SEV.copy(); cmap.set_bad(TISSUE[2]); cmap.set_under(TISSUE[0])
    for j, (k, w) in enumerate(zip(MECH, ("scar", "field", "mixed", "diffuse"))):
        ax = g.cell(1, j, "B" if j == 0 else None, pad=(5, 9, 9, 5)); W = np.where(scar, np.nan, np.where(M[f"w_{w}"] > 0, M[f"w_{w}"], -1.0))
        im = ax.imshow(grid(M, "calibrated", W), origin="lower", cmap=cmap, vmin=0.0, vmax=1, interpolation="nearest")
        ax.set(xticks=[], yticks=[])
        for sp in ax.spines.values(): sp.set_visible(True)
        if j == 0:
            ax.set_ylabel("Base → apex", fontsize=6)
            for c, l in ((TISSUE[2], "Scar (carries the lesion in the four cases)"), (TISSUE[0], "Preserved lung without lesion")): ax.plot([], [], "s", ms=4.5, color=c, mec=DARK, mew=0.3, label=l)
            ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.04), ncol=2, fontsize=5.6, handletextpad=0.3, columnspacing=1.2)
    x0, y0, x1, y1 = g.box(1, 3); cax = g.fig.add_axes([(x1 - 7.5) / g.W, (y0 + 15) / g.H, 1.4 / g.W, (y1 - y0 - 27) / g.H]); cb = g.fig.colorbar(im, cax=cax)
    cb.set_label("Lesion in preserved lung", fontsize=5.4, labelpad=1.5); cb.ax.tick_params(labelsize=5, width=0.4, length=1.5, pad=1); cb.outline.set_linewidth(0.4)
    # C: pressure against lung lost
    PADB = (12.5, 22, 2.5, 6); ax = g.cell(2, 0, "C", pad=PADB); bins = np.array([0.15, 0.3, 0.45, 0.6, 0.75, 0.9])
    for k, (lab, c, tags) in MECH.items():
        Dv = visits(tags); Dv = Dv[Dv.vasc_A > 0]
        mid = [Dv.lost[(Dv.lost >= a) & (Dv.lost < b)].mean() for a, b in zip(bins[:-1], bins[1:])]; med = [Dv.mPAP[(Dv.lost >= a) & (Dv.lost < b)].median() for a, b in zip(bins[:-1], bins[1:])]
        ax.plot(100 * np.array(mid), med, "-o", ms=2.4, lw=1.1, color=c)
    ax.axhline(25, color=DARK, lw=0.5, ls=":"); ax.set(xlabel="Lung lost (%)", ylabel="mPAP, cases with a lesion (mmHg)", ylim=(14, 52))
    keys = list(MECH); x = np.arange(len(keys)); cols = [MECH[k][1] for k in keys]; labs = [MECH[k][0] for k in keys]; rf, sev, pres = [], [], []
    for k in keys:
        tags = MECH[k][2]; Dk = final_visits(F, tags[0]); Mi = Dk[(Dk.FVC >= 50) & (Dk.FVC <= 90)].reset_index(drop=True)
        ci = bootstrap(Mi); rf.append((point(Mi)["rFVC"], *ci["rFVC"]))
        P = pd.concat([final_visits(F, t) for t in tags], ignore_index=True); A5 = P[P.FVC < 50]
        sev.append(prevalence_ci(A5, lambda d: d.mPAP > 40) if len(tags) > 1 else (0.0, 0.0, 0.0))
        hi = P[P.FVC >= 70]; pres.append(prevalence_ci(hi, lambda d: d.mPAP >= 25) + (len(hi),))

    def bars(ax, vals, ylabel, ref=None, reflab=None, ylim=None, extra=None):
        v = np.array([a[0] for a in vals]); lo_ = np.array([a[1] for a in vals]); hi_ = np.array([a[2] for a in vals])
        ax.bar(x, v, 0.62, color=cols); ax.errorbar(x, v, yerr=[v - lo_, hi_ - v], fmt="none", ecolor=DARK, elinewidth=0.7, capsize=1.5)
        xt, xl = list(x), list(labs)
        if extra is not None:
            ax.bar(len(keys), extra[0], 0.62, color="white", edgecolor=DARK, lw=0.9); ax.text(len(keys), extra[0] + 2, extra[1], ha="center", fontsize=5.4); xt.append(len(keys)); xl.append("Patients")
        if ref is not None:
            ax.axhline(ref, color=DARK, lw=0.9, ls="--"); ax.plot([], [], color=DARK, lw=0.9, ls="--", label=reflab); legend_below(ax, dy=-0.34)
        ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=5.4, rotation=35, ha="right", rotation_mode="anchor"); ax.set(ylabel=ylabel)
        if ylim: ax.set_ylim(*ylim)
    ax = g.cell(2, 1, "D", pad=PADB); bars(ax, rf, "r(mPAP, FVC), FVC 50–90%", 0.0, "Patients with IPF", (-0.55, 0.2))
    ax = g.cell(2, 2, "E", pad=PADB); bars(ax, sev, "mPAP > 40 mmHg, FVC < 50% (%)", 9, "Transplant candidates with IPF", (0, 45))
    extra = None
    if clinical is not None:
        V = clinical.dropna(subset=["FVC", "mPAP"]); V = V[V.FVC >= 70]; extra = (100 * (V.mPAP >= 25).mean(), f"{int((V.mPAP >= 25).sum())} of {len(V)}")
    ax = g.cell(2, 3, "F", pad=PADB); bars(ax, pres, "mPAP ≥ 25 mmHg, FVC ≥ 70% (%)", None, None, (0, 82), extra)
    g.save("Story3")
    for k, a, b, c in zip(keys, rf, sev, pres):
        print(f"  {k:8s} rFVC {a[0]:.2f} ({a[1]:.2f} to {a[2]:.2f}) | >40 mmHg advanced {b[0]:.0f}% ({b[1]:.0f}-{b[2]:.0f}) | PH with FVC>=70: {c[0]:.0f}% (n = {c[3]})")


# ---------------------------------------------------------------- Story4: UIP and NSIP
def story4(clinical):
    M = maps(); g = Grid(180, 122, heights=[1, 1], widths=[1, 1, 1]); cm3 = ListedColormap(TISSUE)
    x0, y0, x1, y1 = g.box(0, 0); g.fig.text(x0 / g.W, y1 / g.H, "A", fontsize=9, fontweight="bold", va="top", ha="left"); side = min((x1 - x0 - 12) / 2, y1 - y0 - 20)
    for j, pat in enumerate(("uip", "nsip")):
        ax = g.fig.add_axes([(x0 + 6 + j * (side + 4)) / g.W, (y0 + 13) / g.H, side / g.W, side / g.H])
        ax.imshow(grid(M, pat, M[f"{pat}_state"].astype(float)), origin="lower", cmap=cm3, vmin=0, vmax=2, interpolation="nearest"); ax.set(xticks=[], yticks=[])
        ax.set_title(PAT[pat][0], fontsize=6.6, color=PAT[pat][1], fontweight="bold", pad=2)
        for sp in ax.spines.values(): sp.set_visible(True)
        if j == 0:
            ax.set_ylabel("Base → apex", fontsize=6)
            for c, l in zip(TISSUE, ["Preserved", "Fibrotic", "Honeycomb"]): ax.plot([], [], "s", ms=4.5, color=c, mec=DARK, mew=0.3, label=l)
            ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.05), ncol=3, fontsize=5.6, handletextpad=0.2, columnspacing=0.9)
    ax = g.cell(0, 1, "B"); keys = ("lt20", "b20_25", "ge25", "ge35"); tgt = (51, 30, 19, 4); w = 0.24
    for j, (pat, (lab, c, tags)) in enumerate(PAT.items()):
        D = final_visits(F, tags[0]); Mi = D[(D.FVC >= 50) & (D.FVC <= 90)].reset_index(drop=True); pt, ci = point(Mi), bootstrap(Mi)
        v = np.array([pt[k] for k in keys]); lo = np.array([ci[k][0] for k in keys]); hi = np.array([ci[k][1] for k in keys])
        ax.bar(np.arange(4) + (j - 1) * w, v, w * 0.92, color=c, label=f"{lab} (n = {len(Mi)})"); ax.errorbar(np.arange(4) + (j - 1) * w, v, yerr=[v - lo, hi - v], fmt="none", ecolor=DARK, elinewidth=0.6, capsize=1.2)
    ax.plot(np.arange(4), tgt, "D", ms=4, mfc="white", mec=DARK, mew=0.9, label="Patients with IPF", zorder=5)
    ax.set_xticks(range(4)); ax.set_xticklabels(["< 20", "20–25", "≥ 25", "≥ 35"]); ax.set(xlabel="mPAP (mmHg), FVC 50–90%", ylabel="Share of cases (%)", ylim=(0, 100)); ax.legend(loc="upper right", fontsize=5.6, handlelength=1)
    ax = g.cell(0, 2, "C", pad=(12.5, 14, 2.5, 6)); nadv = []
    for j, (pat, (lab, c, tags)) in enumerate(PAT.items()):
        P = pd.concat([final_visits(F, t) for t in tags], ignore_index=True); A5 = P[P.FVC < 50]; v, lo, hi = prevalence_ci(A5, lambda d: d.mPAP >= 25)
        ax.bar(j, v, 0.62, color=c); ax.errorbar(j, v, yerr=[[v - lo], [hi - v]], fmt="none", ecolor=DARK, elinewidth=0.7, capsize=1.5); nadv.append(len(A5))
    ax.axhline(46, color=DARK, lw=0.9, ls="--"); ax.plot([], [], color=DARK, lw=0.9, ls="--", label="Transplant candidates with IPF")
    ax.set_xticks(range(3)); ax.set_xticklabels([f"{v[0]}\n(n = {n})" for v, n in zip(PAT.values(), nadv)], fontsize=5.6); ax.set(ylabel="mPAP ≥ 25 mmHg, FVC < 50% (%)", ylim=(0, 75)); ax.legend(loc="upper right", fontsize=5.6, handlelength=1.6)
    T = nt.trajectory_points(tags=("big2_uip", "adv2_uip", "big2_nsip", "adv2_nsip")); L = np.linspace(0.08, 0.75, 30)
    ax = g.cell(1, 0, "D", pad=(12.5, 10, 9, 6))
    for pat, ls in (("uip", "-"), ("nsip", "--")):
        Fs, _ = nt.surface(T[T.tag.str.endswith(pat)])
        for a in (0, 2, 4):
            y = Fs.predict(np.c_[L, np.full(len(L), a)]); ax.plot(100 * L, y, ls, color=PAT[pat][1], lw=1.1)
            if pat == "uip": ax.text(77, y[-1], f"{a}", fontsize=5.8, va="center")
    ax.text(79, 36.5, "Lesion", fontsize=5.8, ha="center")
    for pat, ls in (("uip", "-"), ("nsip", "--")): ax.plot([], [], ls, color=PAT[pat][1], lw=1.1, label=PAT[pat][0])
    ax.axhline(25, color=DARK, lw=0.5, ls=":"); ax.set(xlabel="Lung lost (%)", ylabel="mPAP (mmHg)", ylim=(12, 39), xlim=(5, 75)); ax.legend(loc="upper left", fontsize=5.8, handlelength=1.8)
    ax = g.cell(1, 1, "E")
    for pat in ("uip", "nsip"):
        lab, c, tags = PAT[pat]; D = pd.concat([final_visits(F, t) for t in tags], ignore_index=True)
        ax.plot(D.FVC, D.DLCO, "o", ms=1.8, color=c, mew=0, alpha=0.55); trend(ax, D.FVC, D.DLCO, color=c); ax.plot([], [], "o", ms=3, color=c, label=f"{lab} (n = {len(D)})")
    ax.set(xlabel="FVC (% of baseline)", ylabel="DLCO (% of baseline)", xlim=(20, 90), ylim=(0, 75)); ax.legend(loc="upper left", fontsize=5.8, handlelength=1)
    ax = g.cell(1, 2, "F"); rng = np.random.default_rng(0)
    for j, pat in enumerate(("uip", "nsip")):
        y = np.array([d.years.iloc[-1] for _, d in cases(F, PAT[pat][2][0])]); ax.plot(j + rng.uniform(-0.18, 0.18, len(y)), y, "o", ms=1.6, color=PAT[pat][1], mew=0, alpha=0.5)
        ax.plot([j - 0.27, j + 0.27], [np.median(y)] * 2, color=DARK, lw=1.2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["UIP", "NSIP"]); ax.set(ylabel="Years from healthy lung to the visit", xlim=(-0.5, 1.5), yscale="log"); ax.set_yticks([1, 2, 5, 10, 20, 50]); ax.set_yticklabels(["1", "2", "5", "10", "20", "50"])
    g.save("Story4")


# ---------------------------------------------------------------- Story5: catheterised patients
def story5(clinical):
    if clinical is None:
        print("Story5 skipped: needs --clinical"); return
    V = clinical.dropna(subset=["FVC", "mPAP"]).reset_index(drop=True); GC = {"UIP": (PATTERN["uip"], "o"), "NSIP": (PATTERN["nsip"], "s"), "other": (NEUTRAL, "^")}
    g = Grid(180, 130, heights=[1, 1], widths=[1, 1])
    for c_, col, lab, letter in ((0, "FVC", "FVC (% predicted)", "A"), (1, "DLCO", "DLCO (% predicted)", "B")):
        ax = g.cell(0, c_, letter, pad=(12.5, 15, 2.5, 6)); W = V.dropna(subset=[col])
        for gname, (c, mk) in GC.items():
            x = W[W.group == gname]; ax.plot(x[col], x.mPAP, mk, ms=3.2, color=c, mew=0, label=f"{gname if gname != 'other' else 'Other patterns'} (n = {len(x)})")
        r, p = trend(ax, W[col], W.mPAP); ax.axhline(25, color=DARK, lw=0.5, ls=":"); stat(ax, r, p, len(W))
        ax.set(xlabel=lab, ylabel="mPAP, catheterisation (mmHg)", ylim=(5, 82)); legend_below(ax, dy=-0.2, ncol=3)
    ax = g.cell(1, 0, "C", pad=(12.5, 19, 2.5, 6)); bands = [(20, 50), (50, 70), (70, 101)]; xb = np.arange(3); w = 0.25
    for j, k in enumerate(("field", "diffuse")):
        P = pd.concat([final_visits(F, t) for t in MECH[k][2]], ignore_index=True)
        ax.bar(xb + (j - 1) * w, [100 * (P[(P.FVC >= a) & (P.FVC < b)].mPAP >= 25).mean() for a, b in bands], w * 0.92, color=MECH[k][1], label="Simulated, lesion spreading from the scar" if k == "field" else "Simulated, diffuse lesion")
    o = [(int((V[(V.FVC >= a) & (V.FVC < b)].mPAP >= 25).sum()), int(((V.FVC >= a) & (V.FVC < b)).sum())) for a, b in bands]
    ax.bar(xb + w, [100 * a / b for a, b in o], w * 0.92, color="white", edgecolor=DARK, lw=0.9, label="Patients")
    for i, (a, b) in enumerate(o): ax.text(xb[i] + w, 100 * a / b + 2, f"{a} of {b}", ha="center", fontsize=5.6)
    ax.set_xticks(xb); ax.set_xticklabels(["< 50", "50–70", "≥ 70"]); ax.set(xlabel="FVC (%)", ylabel="mPAP ≥ 25 mmHg (%)", ylim=(0, 112)); ax.set_yticks([0, 25, 50, 75, 100]); legend_below(ax, dy=-0.22, ncol=1)
    ax = g.cell(1, 1, "D", pad=(12.5, 19, 2.5, 6)); f = os.path.join(OUT, "TableS24_individual_report.tsv")
    if os.path.exists(f):
        J = pd.read_csv(f, sep="\t").sort_values("mPAP_RHC").reset_index(drop=True); x = np.arange(len(J))
        ax.bar(x, J.baseline, color=LIGHT, width=0.85, label="Aged lung"); ax.bar(x, J.from_lung_lost, bottom=J.baseline, color=NEUTRAL, width=0.85, label="Lung lost")
        ax.bar(x, J.from_lesion, bottom=J.baseline + J.from_lung_lost, color=LESION, width=0.85, label="Vascular lesion"); ax.plot(x, J.mPAP_RHC, "o", ms=2.3, mfc="white", mec=DARK, mew=0.6, label="Catheterisation")
        ax.axhline(25, color=DARK, lw=0.5, ls=":"); ax.set(xlabel=f"Patients, ranked by measured mPAP (n = {len(J)})", ylabel="mPAP (mmHg)", xticks=[], ylim=(0, 72)); legend_below(ax, dy=-0.14, ncol=2)
    g.save("Story5")


# ---------------------------------------------------------------- Story6: clinical consequence
def story6(clinical):
    g = Grid(180, 126, heights=[1, 1], widths=[1, 1, 1]); rng = np.random.default_rng(1)
    D = pd.read_csv(os.path.join(OUT, "pressure_flow_cases.tsv"), sep="\t"); D["vasc"] = D.vasc_A > 0; N = D[D.mPAP_rest < 25]
    ax = g.cell(0, 0, "A")
    for v, c in ((False, NEUTRAL), (True, LESION)):
        for _, r in N[N.vasc == v].iterrows(): ax.plot([r.Q_rest, r["Q_7.5"], r["Q_10"]], [r.mPAP_rest, r["mPAP_7.5"], r["mPAP_10"]], "-o", ms=1.5, lw=0.45, color=c, alpha=0.6, mew=0)
        ax.plot([], [], "-o", ms=2.5, lw=0.8, color=c, label=("Lesion" if v else "No lesion") + f" (n = {(N.vasc == v).sum()})")
    ax.axhline(25, color=DARK, lw=0.5, ls=":"); ax.set(xlabel="Cardiac output (L/min)", ylabel="mPAP (mmHg)", xticks=[5, 7.5, 10], ylim=(15, 40)); ax.legend(loc="upper left", fontsize=5.8, handlelength=1.4)
    ax = g.cell(0, 1, "B"); x = 0
    for Q in (7.5, 10.0):
        for v, c in ((False, NEUTRAL), (True, LESION)):
            y = (N[N.vasc == v][f"mPAP_{Q:g}"] >= 25).values; b = [100 * rng.choice(y, len(y)).mean() for _ in range(2000)]; m = 100 * y.mean()
            ax.bar(x, m, 0.8, color=c); ax.errorbar(x, m, yerr=[[m - np.percentile(b, 2.5)], [np.percentile(b, 97.5) - m]], fmt="none", ecolor=DARK, elinewidth=0.7, capsize=1.5)
            ax.text(x, np.percentile(b, 97.5) + 2.5, f"{int(y.sum())} of {len(y)}", ha="center", fontsize=5.4); x += 1
        x += 0.6
    ax.set_xticks([0.5, 3.1]); ax.set_xticklabels(["7.5 L/min", "10 L/min"]); ax.set(ylabel="Cases reaching 25 mmHg (%)", ylim=(0, 118)); ax.set_yticks([0, 25, 50, 75, 100])
    # C: how the estimate is made
    ax = g.cell(0, 2, "C", pad=(2, 3, 2, 7), axis=False); ax.set(xlim=(0, 100), ylim=(0, 100))
    box_diagram(ax, 2, 70, 40, 26, "Required\nFVC, DLCO", fc=LIGHT); box_diagram(ax, 2, 6, 40, 56, "Optional\necho sPAP\nPA diameter (CT)\nBNP\nRV area change\nRA area")
    box_diagram(ax, 52, 60, 46, 36, "Simulated visits\nof the model\n(lung lost, lesion,\nFVC, DLCO, mPAP)", fc="#D6E9F5", ec=MODEL)
    box_diagram(ax, 52, 6, 46, 40, "For the patient\nprobability of PH\nmPAP, 90% interval\nlung lost and lesion", ec=MODEL, bold=False)
    arrow(ax, 42, 83, 52, 80); arrow(ax, 42, 50, 52, 66); arrow(ax, 75, 60, 75, 46)
    f = os.path.join(OUT, "TableS22_noninvasive_tool_cv.tsv")
    if clinical is None or not os.path.exists(f):
        print("Story6, panels D-F skipped: need --clinical and the tables of noninvasive_tool.py"); g.save("Story6"); return
    S = pd.read_csv(f, sep="\t"); S = S[S.inputs.isin(["FVC + DLCO", "+ echo sPAP", "+ BNP", "+ PA diameter", "+ PA diameter + BNP", "+ all optional"])].reset_index(drop=True)
    ax = g.cell(1, 0, "D", pad=(27, 18, 2.5, 6)); y = np.arange(len(S))[::-1]
    ax.errorbar(S.AUC_LuPNet, y, xerr=[S.AUC_LuPNet - S.AUC_min, S.AUC_max - S.AUC_LuPNet], fmt="o", ms=3.5, color=MODEL, elinewidth=0.6, capsize=0, label=f"Model-based estimate (n = {int(S.n.iloc[0])})")
    ax.plot(S.AUC_logistic_complete, y, "D", ms=3.2, mfc="white", mec=DARK, mew=0.8, label="Logistic regression, complete cases (n)")
    for yy, n in zip(y, S.n_complete): ax.text(0.99, yy, str(int(n)), fontsize=5.4, va="center", ha="right")
    lab = {"FVC + DLCO": "FVC + DLCO", "+ echo sPAP": "+ echo sPAP", "+ BNP": "+ BNP", "+ PA diameter": "+ PA diameter (CT)", "+ PA diameter + BNP": "+ PA diameter + BNP", "+ all optional": "+ all optional inputs"}
    ax.set_yticks(y); ax.set_yticklabels([lab[i] for i in S.inputs], fontsize=5.8); ax.set(xlim=(0.55, 1.0), ylim=(-0.6, len(S) - 0.4), xlabel="AUC for mPAP ≥ 25 mmHg")
    ax.legend(loc="upper right", bbox_to_anchor=(1.03, -0.27), fontsize=5.6, handlelength=1)
    cv0 = pd.read_csv(os.path.join(OUT, "cv_FVC_plus_DLCO.tsv"), sep="\t"); cv1 = pd.read_csv(os.path.join(OUT, "cv_plus_all_optional.tsv"), sep="\t")
    ax = g.cell(1, 1, "E", pad=(12.5, 18, 2.5, 6))
    for k, df in enumerate((cv0, cv1)):
        for j, ph in enumerate((False, True)):
            v = df.pPH[(df.mPAP >= 25) == ph]; ax.plot(k * 2.4 + j + rng.uniform(-0.15, 0.15, len(v)), v, "o", ms=2.6, color=LESION if ph else NEUTRAL, mew=0, alpha=0.85)
    ax.axhline(0.7, color=DARK, lw=0.4, ls=":"); ax.axhline(0.3, color=DARK, lw=0.4, ls=":")
    ax.set_xticks([0, 1, 2.4, 3.4]); ax.set_xticklabels(["No PH", "PH", "No PH", "PH"], fontsize=6); ax.set(ylabel="Estimated probability of PH", ylim=(-0.03, 1.03), xlim=(-0.6, 4), yticks=[0, 0.3, 0.5, 0.7, 1])
    ax.text(0.5, -0.24, "FVC + DLCO", ha="center", fontsize=6, transform=ax.get_xaxis_transform()); ax.text(2.9, -0.24, "+ all optional inputs", ha="center", fontsize=6, transform=ax.get_xaxis_transform())
    ax = g.cell(1, 2, "F", pad=(12.5, 18, 2.5, 6)); e = cv1.est - cv1.mPAP
    ax.errorbar(cv1.mPAP, cv1.est, yerr=[cv1.est - cv1.lo, cv1.hi - cv1.est], fmt="o", ms=2.8, color=MODEL, ecolor=MODEL, elinewidth=0.4, alpha=0.85)
    ax.plot([8, 76], [8, 76], color=DARK, lw=0.6, ls=":"); trend(ax, cv1.mPAP, cv1.est); r, p = pearsonr(cv1.mPAP, cv1.est)
    ax.text(0.03, 0.97, f"r = {r:.2f}, {ptxt(p)}, n = {len(cv1)}\nmean error {e.abs().mean():.1f} mmHg", transform=ax.transAxes, va="top", fontsize=5.8)
    ax.set(xlabel="mPAP, catheterisation (mmHg)", ylabel="mPAP, estimated (mmHg)", xlim=(8, 76), ylim=(8, 76))
    g.save("Story6")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--clinical", default=""); ap.add_argument("--only", nargs="*", default=[])
    a = ap.parse_args(); os.makedirs(OUT, exist_ok=True)
    C = nt.load_clinical(a.clinical) if a.clinical else None
    for name, fn in (("2", story2), ("3", story3), ("4", story4), ("5", story5), ("6", story6)):
        if not a.only or name in a.only:
            fn(C)
