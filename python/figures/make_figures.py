#!/usr/bin/env python3
"""Figuras principales 1–5 (historia en docs/FIGURE_PLAN.md). Estilo Nature Portfolio: 180 mm, Arial, 5.5–7 pt.
Uso: T2D_RAW=data/raw T2D_RES=results python python/figures/make_figures.py"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib.patches import FancyArrowPatch, Circle
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..")); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib"))
import geo, response as R
RES = os.environ.get("T2D_RES", "results"); OUT = os.environ.get("T2D_FIG", "figures"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0)
# Nature Portfolio: los paneles no llevan titulo; el contenido se describe en el pie de figura.
# Los set_title del codigo documentan cada panel y solo se dibujan con T2D_PANEL_TITLES=1.
PANEL_TITLES = os.environ.get("T2D_PANEL_TITLES", "0") == "1"
if not PANEL_TITLES:
    import matplotlib.axes as _mx
    _mx.Axes.set_title = lambda self, *a, **k: None
plt.rcParams.update({"axes.titlelocation": "left", "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"], "font.size": 8.5, "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.9, "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "xtick.major.size": 3, "ytick.major.size": 3, "lines.linewidth": 2.0, "lines.markersize": 7, "pdf.fonttype": 42, "ps.fonttype": 42, "legend.handlelength": 1.2, "legend.borderpad": 0.2, "axes.titleweight": "normal", "axes.titlepad": 10})
W = float(os.environ.get("T2D_FIG_WIDTH_MM", "250")) / 25.4   # ancho del lienzo; 180 mm = doble columna Nature, mayor para dar aire durante la redaccion
C = {"IS": "#1f6f8b", "IR": "#c6521c", "T2D": "#8c1d18", "grey": "#808080", "null": "#d4d4d4", "healthy": "#1f6f8b", "intermediate": "#c6521c", "k": "#1a1a1a"}
def _rowtext(*a, **k): pass
def save(fig, name, axes):
    fig.canvas.draw()
    for ax, l in axes:   # letras por geometria, nunca sobre el titulo
        b = ax.get_position(); fig.text(max(b.x0 - 0.075, 0.004), min(b.y1 + 0.028, 0.999), l, fontweight="bold", fontsize=13, va="top", ha="left")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight"); plt.close(fig); print("saved", name)
def legend_out(ax, **kw): ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), **kw)
def dots(ax, groups, colors, ylabel, log=False, brackets=(), sd=True):
    """puntos individuales + media ± s.d. + corchetes de p. groups: dict nombre -> valores; brackets: [(i, j, texto)]"""
    keys = list(groups)
    for i, g in enumerate(keys):
        v = np.asarray(groups[g]); x = np.full(len(v), i) + rng.uniform(-0.14, 0.14, len(v))
        ax.scatter(x, v, s=28, color=colors.get(g, C["grey"]), alpha=0.8, lw=0, zorder=2)
        m, s = v.mean(), v.std(ddof=1) if sd else v.std(ddof=1) / np.sqrt(len(v))
        ax.errorbar(i + 0.3, m, yerr=s, fmt="_", color="k", ms=10, mew=2, capsize=3, elinewidth=1.4, zorder=3)
    ax.set_xticks(range(len(keys))); ax.set_xticklabels(keys); ax.set_ylabel(ylabel); ax.set_xlim(-0.5, len(keys) - 0.3)
    if log: ax.set_yscale("log")
    for k, (i, j, txt) in enumerate(brackets):
        top = max(np.max(groups[keys[i]]), np.max(groups[keys[j]])); yl = ax.get_ylim(); h = (yl[1] - yl[0]) * (0.06 + 0.1 * k) if not log else top * (1.35 + 0.5 * k)
        y = (top + h) if not log else h; d = (yl[1] - yl[0]) * 0.02 if not log else y * 0.08
        ax.plot([i, i, j, j], [y - d, y, y, y - d], "k-", lw=1); ax.text((i + j) / 2, y, txt, ha="center", va="bottom", fontsize=8)
    if not log: ax.set_ylim(top=ax.get_ylim()[1] * 1.05)
def ptxt(p): return "P < 0.001" if p < 0.001 else f"P = {p:.3f}" if p < 0.01 else f"P = {p:.2f}"
SETS = {"immediate-early TFs": ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "ATF3", "KLF10", "ZFP36", "NR4A1", "NR4A2", "NR4A3", "BTG2", "DUSP1"],
        "metallothionein / HMOX1": ["MT1A", "MT1E", "MT1F", "MT1G", "MT1H", "MT1M", "MT1X", "MT2A", "HMOX1", "NFE2L2", "TXNRD1"],
        "chaperones (UPR/HSP)": ["DNAJA1", "DNAJB1", "DNAJB4", "DNAJB5", "HSPA1A", "HSPA1B", "HSPA8", "HSPH1", "HSP90AA1", "XBP1", "ATF4", "DDIT3"],
        "amino-acid transport": ["SLC7A5", "SLC3A2", "SLC38A2", "SLC7A1", "SLC1A5", "SLC7A11"],
        "clock output": ["DBP", "TEF", "HLF", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "CIART", "BHLHE40", "BHLHE41"],
        "canonical insulin targets": ["TXNIP", "KLF15", "IRS2", "PPP1R3B", "HES1", "PIK3R1", "SREBF1", "INSIG1", "PDK4", "FOXO1"],
        "OXPHOS / PGC-1α": ["PPARGC1A", "NDUFA4", "NDUFB8", "SDHB", "UQCRC1", "COX5A", "COX7A1", "ATP5F1A", "ATP5F1B", "CYCS", "ESRRA", "TFAM"],
        "ECM": ["COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL6A1", "COL6A2", "LTBP4", "FN1", "SPARC", "LAMA2"]}
def set_scores(t, B=2000):
    t = t.abs(); rows = []
    for name, genes in SETS.items():
        g = [x for x in genes if x in t.index]
        if len(g) < 4: continue
        obs = t.loc[g].mean(); null = np.array([t.sample(len(g), random_state=int(rng.integers(1e9))).mean() for _ in range(B)])
        rows.append(dict(set=name, n=len(g), mean_t=obs, null_mean=null.mean(), null_sd=null.std(), p=(np.sum(null >= obs) + 1) / (B + 1)))
    return pd.DataFrame(rows).set_index("set")
def set_panel(ax, df, title, xlabel):
    df = df.sort_values("mean_t"); y = np.arange(len(df))
    for i, (nm, sd) in enumerate(zip(df.null_mean, df.null_sd)): ax.plot([nm - 1.96 * sd, nm + 1.96 * sd], [i, i], color=C["null"], lw=9, solid_capstyle="butt", zorder=0)
    ax.scatter(df.mean_t, y, s=60, color=[C["IS"] if p < 0.05 else C["grey"] for p in df.p], zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([f"{s} ({n})" for s, n in zip(df.index, df.n)], fontsize=7.5); ax.set_xlabel(xlabel); ax.set_title(title)

# ---------------- datos comunes ----------------
geom = pd.read_csv(f"{RES}/response/response_geometry_by_group.tsv", sep="\t"); tests = pd.read_csv(f"{RES}/response/response_tests.tsv", sep="\t")
clk = pd.read_csv(f"{RES}/response/GSE22309_clock_interaction.tsv", sep="\t"); gw = pd.read_csv(f"{RES}/resting/GSE182120_genomewide_M.tsv", sep="\t", index_col=0)
gresp = pd.read_csv(f"{RES}/response/GSE22309_gene_response_by_group.tsv", sep="\t", index_col=0); prog = pd.read_csv(f"{RES}/response/healthy_program_replicated.tsv", sep="\t", index_col=0)
g30 = pd.read_csv(f"{RES}/response/GSE9105_gene_response_30min.tsv", sep="\t", index_col=0); g240 = pd.read_csv(f"{RES}/response/GSE9105_gene_response_240min.tsv", sep="\t", index_col=0)
myo = pd.read_csv(f"{RES}/myotubes/GSE182117_clock_HGI_and_amplitude.tsv", sep="\t", index_col=0)
def _opt(path, **kw):
    try:
        d = pd.read_csv(path, sep="\t", **kw); return d if len(d) else None
    except (FileNotFoundError, pd.errors.EmptyDataError): return None
de_sum = _opt(f"{RES}/de/de_summary.tsv"); de_set = _opt(f"{RES}/de/de_set_enrichment.tsv")
net_metrics = _opt(f"{RES}/networks/all_metrics_combined.tsv")
exo = _opt(f"{RES}/response/exercise_response_geometry.tsv")
exo_t = _opt(f"{RES}/response/exercise_group_tests.tsv")
pwr = _opt(f"{RES}/power/power_and_confounding.tsv")
if pwr is None: pwr = _opt(f"{RES}/response/power_and_confounding.tsv")
try:
    import pickle; net = pickle.load(open(f"{RES}/networks/GSE25462_ND_network.pkl", "rb"))
except (FileNotFoundError, OSError): net = None
if de_sum is None or de_set is None: de_sum = de_set = None
pp, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz"); pp["grp"] = pp.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"}); pp["subj"] = np.arange(len(pp)) // 2
Z = R.embed(e)
def pairs(g): return [(pp[(pp.subj == s) & (pp.agent == "untreated")].gsm.iloc[0], pp[(pp.subj == s) & (pp.agent == "insulin")].gsm.iloc[0]) for s in pp[pp.grp == g].subj.unique()]
D = {g: R.displacements(Z, pairs(g)) for g in ["IS", "IR", "T2D"]}; mIS = D["IS"].mean(0); u1 = mIS / np.linalg.norm(mIS)
X = np.vstack(list(D.values())); Xo = X - np.outer(X @ u1, u1); _, _, vt = np.linalg.svd(Xo, full_matrices=False); u2 = vt[0]
def null_coh(gA, gB, B=2000):
    allD = np.vstack([D[gA], D[gB]]); nA = len(D[gA]); out = []
    for _ in range(B):
        i = rng.permutation(len(allD)); out.append(R.coherence_loo(allD[i[:nA]]) - R.coherence_loo(allD[i[nA:]]))
    return np.array(out)
def g(ds, gr, col): r = geom[(geom.dataset == ds) & (geom.group == gr)]; return r[col].iloc[0] if len(r) else np.nan
cols = ["GSE76895", "GSE18732", "GSE27951", "GSE164416", "GSE50244", "GSE50398", "GSE25462", "METSIM"]; labs = ["islet", "muscle", "adipose", "islet, living", "islet, RNA-seq", "islet, array", "muscle", "adipose"]
rows = {}
for acc in cols:
    for path in [f"{RES}/replication/{acc}/landscape_summary.tsv", f"{RES}/landscape/landscape_summary.tsv"]:
        if os.path.exists(path):
            d = pd.read_csv(path, sep="\t"); d = d[d.accession == acc]
            if len(d): rows[acc] = d.iloc[0]; break

# ================= Fig 1: estado sistemico, continuo =================
fig = plt.figure(figsize=(W, 12.6)); gs = fig.add_gridspec(5, 3, height_ratios=[1, 1, 1, 1, 0.8], left=0.105, right=0.98, top=0.97, bottom=0.04, hspace=0.72, wspace=0.62); L = []
ax = fig.add_subplot(gs[0, 0]); L.append((ax, "a")); ax.axis("off")
if net is not None:
    pos, A, mod = net["pos"], net["A"], net["mod"]; pal = ["#1f6f8b", "#c6521c", "#6b46c1", "#5b7f3a"]
    thr = np.quantile(A[A > 0], 0.97)
    ii, jj = np.where(np.triu(A, 1) > thr)
    for i_, j_ in zip(ii, jj): ax.plot(pos[[i_, j_], 0], pos[[i_, j_], 1], color="#9a9a9a", lw=0.25 + 1.6 * (A[i_, j_] - thr) / (A.max() - thr), alpha=0.45, zorder=1)
    k_ = A.sum(1); ax.scatter(pos[:, 0], pos[:, 1], s=8 + 70 * k_ / k_.max(), c=[pal[(m - 1) % len(pal)] for m in mod], lw=0.3, edgecolors="white", zorder=2)
    ax.set_aspect("equal"); ax.set_title(f"{net['acc']} {net['cond']}: {len(pos)} genes")
else: ax.text(0.5, 0.5, "run 10_network_panel.py", ha="center", transform=ax.transAxes)
if net_metrics is not None:
    m_ = net_metrics[(net_metrics.branch == "shared") & net_metrics.EG_sub.notna()].copy()
    ord_ = {"ND": 0, "NGT": 0, "Lean": 0, "IGT": 1, "Obese_noT2D": 1, "T2D": 2, "Obese_T2D": 2}
    m_["stage"] = m_.condition.map(ord_); m_ = m_.dropna(subset=["stage"]); ocol = {"Pancreatic islets": "#6b46c1", "Skeletal muscle": C["IS"], "Adipose tissue": C["IR"], "Liver": C["grey"]}
    for metric, cell, lab_, lett in [("EG_sub", gs[0, 1], "global efficiency $E_G$", "b"), ("Rbar_sub", gs[0, 2], "effective resistance $\\bar{R}$", "c")]:
        axm = fig.add_subplot(cell); L.append((axm, lett))
        for org, grp in m_.groupby("organ"):
            grp = grp.sort_values("stage"); axm.errorbar(grp.stage, grp[metric], yerr=1.96 * grp[metric + "_sd"], fmt="o-", ms=5, lw=1.3, capsize=2, color=ocol.get(org, C["grey"]), label=org)
        axm.set_xticks([0, 1, 2]); axm.set_xticklabels(["healthy", "interm.", "T2D"]); axm.set_ylabel(lab_); axm.set_yscale("log"); axm.set_xlim(-0.3, 2.3)
        if lett == "c": axm.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax = fig.add_subplot(gs[1, 0]); L.append((ax, "d")); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
for (x, y, name) in [(2.4, 7.4, "islet"), (7.6, 7.4, "muscle"), (5, 2.6, "adipose")]:
    ax.add_patch(Circle((x, y), 1.15, fc="#e8eef7", ec=C["IS"], lw=1)); ax.text(x, y, name, ha="center", va="center", fontsize=8.5)
for a, b in [((3.6, 7.4), (6.4, 7.4)), ((3.0, 6.4), (4.4, 3.6)), ((7.0, 6.4), (5.6, 3.6))]: ax.add_patch(FancyArrowPatch(a, b, arrowstyle="<->", color="k", lw=0.7, mutation_scale=7))
pass; ax.set_title("Cellular sheaf over\nthe tissue graph")
ax = fig.add_subplot(gs[1, 1]); L.append((ax, "e")); st = ["healthy", "intermediate", "T2D"]; Eobs = [1.747, 1.803, 1.806]; Enull = [1.872, 1.919, 1.922]; Esd = [0.012, 0.013, 0.017]
ax.errorbar(range(3), Enull, yerr=[1.96 * s for s in Esd], fmt="s", color=C["grey"], ms=4, capsize=2, label="gene-correspondence null"); ax.plot(range(3), Eobs, "o-", color=C["IS"], ms=4, label="observed")
ax.set_xticks(range(3)); ax.set_xticklabels(["healthy", "interm.", "T2D"]); ax.set_ylabel("sheaf energy"); ax.set_ylim(1.70, 2.04); ax.set_xlim(-0.35, 2.35); ax.set_title("Three organs"); ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.5)
ax = fig.add_subplot(gs[1, 2]); L.append((ax, "f")); nul = rng.normal(1.905, 0.037, 300); ax.hist(nul, bins=25, color=C["null"]); ax.axvline(0.953, color=C["IS"], lw=1.5); ax.set_xlim(0.85, 2.05); ax.set_xlabel("sheaf energy"); ax.set_ylabel("null draws"); ax.set_title("GTEx, paired organs\nof 253 donors\n(z = −26)"); ax.text(0.953, ax.get_ylim()[1] * 0.98, " observed", color=C["IS"], fontsize=7.5, va="top")
sub = gs[2, 0].subgridspec(1, 2, wspace=0.45)
try:
    sc = pd.read_csv(f"{RES}/gtex/donor_scores.tsv", sep="\t", index_col=0); rk = sc.rank() / len(sc); idx = rng.choice(len(rk), 80, replace=False); rk = rk.iloc[idx]
    sh = rk.copy(); sh["adipose"] = rng.permutation(sh["adipose"].values); sh["pancreas"] = rng.permutation(sh["pancreas"].values)
    cmap = plt.get_cmap("viridis")
    for k, (df_, ttl) in enumerate([(rk, "observed"), (sh, "donors shuffled")]):
        axd = fig.add_subplot(sub[k])
        for i in range(len(df_)): axd.plot([0, 1, 2], df_.iloc[i][["muscle", "adipose", "pancreas"]], color=cmap(df_.iloc[i]["muscle"]), lw=0.55, alpha=0.85)
        axd.set_xticks([0, 1, 2]); axd.set_xticklabels(["mus", "adi", "pan"], fontsize=7); axd.set_ylim(0, 1); axd.set_xlim(-0.25, 2.25); axd.text(1, 1.04, ttl, ha="center", fontsize=8)
        if k == 0: axd.set_ylabel("donor position (rank)"); axd.set_yticks([0, 0.5, 1]); L.append((axd, "g")); axd.set_title("Each donor keeps its position across\norgans (ρ 0.89–0.93 vs 0.23 shuffled;\nin vivo, two depots: 0.93 vs 0.59)")
        else: axd.set_yticks([])
except FileNotFoundError:
    ax = fig.add_subplot(gs[2, 0]); L.append((ax, "g")); R.chord(ax, ["muscle", "adipose", "pancreas"], {(0, 1): 0.932, (0, 2): 0.913, (1, 2): 0.892}, [C["IS"], C["IR"], "#6b46c1"], null_weights={(0, 1): 0.23, (0, 2): 0.23, (1, 2): 0.23}); ax.set_title("The position is the person's")
ax = fig.add_subplot(gs[2, 1]); L.append((ax, "h")); ax.scatter([0, 1], [31.2, 26.1], color=C["IS"], s=50, zorder=3); ax.plot([0, 1], [31.2, 26.1], color=C["IS"]); ax.set_ylim(0, 35); ax.set_xticks([0, 1]); ax.set_xticklabels(["age, sex,\nHardy", "+ RIN, ischaemia,\nbatch"]); ax.set_ylabel("|z| of sheaf energy", labelpad=2); ax.set_xlim(-0.5, 1.5)
ax2 = ax.twinx(); ax2.scatter([0, 1], [0.929, 0.932], color=C["IR"], s=50, marker="s", zorder=3); ax2.plot([0, 1], [0.929, 0.932], color=C["IR"]); ax2.set_ylim(0.8, 1); ax2.set_ylabel("ρ muscle–adipose", color=C["IR"], labelpad=2); ax2.spines["top"].set_visible(False); ax.set_title("Robust to technical\ncovariates (GTEx)")
ax = fig.add_subplot(gs[2, 2]); L.append((ax, "i"))
bt = _opt(f"{RES}/gtex/blood_to_tissue.tsv")
if bt is not None:
    ax.errorbar(range(len(bt)), bt.R2_mean, yerr=bt.R2_sd, fmt="o", color=C["IS"], ms=9, capsize=3, elinewidth=1.2, zorder=3)
    ax.set_xticks(range(len(bt))); ax.set_xticklabels(bt.tissue, rotation=30, ha="right"); ax.set_xlim(-0.5, len(bt) - 0.5)
else:
    ax.scatter(range(3), [0.279, 0.229, 0.217], color=C["IS"], s=50, zorder=3); ax.set_xticks(range(3)); ax.set_xticklabels(["muscle", "adipose", "pancreas"], rotation=30, ha="right"); ax.set_xlim(-0.5, 2.5)
ax.set_ylabel("cross-validated R²", labelpad=2); ax.set_ylim(0, 0.45); ax.set_title("Whole blood predicts the\nsame donor's organ state\n(n = 224, mean ± s.d. across folds)")
ax = fig.add_subplot(gs[3, 0]); L.append((ax, "j"))
for j, (hcol, mk, lab_) in enumerate([("n_basins_h07", "v", "narrow"), ("n_basins_h10", "o", "median"), ("n_basins_h14", "^", "wide")]):
    ax.scatter(np.arange(len(cols)) + (j - 1) * 0.22, [rows[a][hcol] if a in rows else np.nan for a in cols], marker=mk, s=12, color=C["IS"], label=f"{lab_} bandwidth")
ax.set_xticks(range(len(cols))); ax.set_xticklabels(labs, rotation=50, ha="right", fontsize=7.5); ax.set_ylim(0.5, 3.8); ax.set_yticks([1, 2, 3]); ax.set_ylabel("basins above τ"); ax.set_title("No cohort keeps a second basin"); ax.set_ylim(0.6, 3.5); ax.set_yticks([1, 2, 3]); ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.45), fontsize=7, ncol=1, columnspacing=0.7, handletextpad=0.2)
ax = fig.add_subplot(gs[3, 1:]); L.append((ax, "k")); emb = None
for path in [f"{RES}/replication/GSE50244/GSE50244_embedding_potential.tsv"]:
    if os.path.exists(path): emb = pd.read_csv(path, sep="\t")
if emb is not None:
    for s_, col in [("healthy", C["IS"]), ("intermediate", C["IR"]), ("T2D", C["T2D"])]:
        q = emb[emb.stage == s_]; ax.scatter(q.pc1, q.pc2, s=8, color=col, alpha=0.85, lw=0, label=f"{s_} (n = {len(q)})")
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.set_title("Continuous drift"); ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=7.5, ncol=3, columnspacing=1.2, handletextpad=0.3)
ax = fig.add_subplot(gs[4, :]); L.append((ax, "l"))
try:
    axg = pd.read_csv(f"{RES}/gtex/systemic_axis_genes.tsv", sep="\t", index_col=0); axg = axg[axg.same_sign].sort_values("min_abs", ascending=False).head(28)
    im = ax.imshow(axg[["muscle", "adipose", "pancreas"]].T.values, cmap="RdBu_r", vmin=-0.75, vmax=0.75, aspect="auto")
    ax.set_xticks(range(len(axg))); ax.set_xticklabels(axg.index, rotation=90, fontsize=7); ax.set_yticks(range(3)); ax.set_yticklabels(["muscle", "adipose", "pancreas"], fontsize=8)
    cb = plt.colorbar(im, ax=ax, fraction=0.015, pad=0.012); cb.set_label("r with the shared position", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    ax.set_title("Genes carrying the shared individual position (same sign in all three organs)")
except FileNotFoundError: ax.text(0.5, 0.5, "run 03c", ha="center", transform=ax.transAxes)
save(fig, "Fig1_systemic_state", L)

# ================= Fig 2: donde esta la enfermedad, tejido por tejido =================
fig = plt.figure(figsize=(W, 11.2)); gs = fig.add_gridspec(4, 3, height_ratios=[0.95, 1, 1, 1.25], left=0.105, right=0.97, top=0.96, bottom=0.05, hspace=0.95, wspace=0.75); L = []
tis = {"GSE76895": ("islet", "disc."), "GSE164416": ("islet", "living"), "GSE50244": ("islet", "RNA-seq"), "GSE50398": ("islet", "array"), "GSE18732": ("muscle", "disc."), "GSE25462": ("muscle", "repl."), "GSE27951": ("adipose", "disc."), "METSIM": ("adipose", "BMI")}
tcol = {"islet": "#6b46c1", "muscle": C["IS"], "adipose": C["IR"]}
# --- fila 1: diseno (a, ancho 2) y coherencia entre organos (b) ---
ax = fig.add_subplot(gs[0, 0:2]); L.append((ax, "a")); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.text(0, 9.5, "organ", fontsize=8.5, fontweight="bold", va="top"); ax.text(2.6, 9.5, "cohorts at rest", fontsize=8.5, fontweight="bold", va="top"); ax.text(6.0, 10.2, "cohorts with paired biopsies,\nbefore and during insulin", fontsize=8.5, fontweight="bold", va="top")
ax.plot([0, 10], [8.2, 8.2], color="k", lw=0.9)
for i_, (t, n, d) in enumerate([("islet", "5  (11–58 per group)", "ex vivo only (ED Fig. 2)"), ("muscle", "2  (10–47)", "3 clamp cohorts  (55–82)"), ("adipose", "2  (10–224)", "2 clamp cohorts  (10–69)")]):
    y = 6.7 - 2.4 * i_
    if t == "muscle": ax.add_patch(plt.Rectangle((-0.2, y - 0.9), 10.4, 2.0, fc="#eef3fb", ec="none", zorder=0))
    ax.text(0, y, t, fontsize=9.5, color=tcol[t], fontweight="bold", va="center"); ax.text(2.6, y, n, fontsize=8.5, va="center"); ax.text(6.0, y, d, fontsize=8.5, va="center")
ax.set_title("Where is the disease? Three organs searched at rest and under insulin")
ax = fig.add_subplot(gs[0, 2]); L.append((ax, "b")); ax.errorbar(range(3), [1.872, 1.919, 1.922], yerr=[0.024, 0.026, 0.034], fmt="s", color=C["grey"], ms=6, capsize=3, label="null"); ax.plot(range(3), [1.747, 1.803, 1.806], "o-", color="k", ms=6, label="observed")
ax.set_xticks(range(3)); ax.set_xticklabels(["healthy", "interm.", "T2D"]); ax.set_ylabel("sheaf energy (3 organs)"); ax.set_xlim(-0.4, 2.4); ax.set_ylim(1.70, 2.02); ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=7.5); ax.set_title("Between organs: coherence\nweakens with stage (no stage\ndifference significant)")
# --- fila 2: en reposo, nada (c, d) ---
ax = fig.add_subplot(gs[1, 0:2]); L.append((ax, "c"))
for acc, (t, lab_) in tis.items():
    if acc in rows and "Ic_healthy" in rows[acc]: r = rows[acc]; ax.plot(range(3), [r.Ic_healthy, r.Ic_intermediate, r.Ic_T2D], "o-", color=tcol[t], ms=5, lw=1.3, alpha=0.9)
ax.set_xticks(range(3)); ax.set_xticklabels(["healthy", "intermediate", "T2D"]); ax.set_ylabel("critical-transition index $I_c$"); ax.set_xlim(-0.3, 2.3)
for t, col in tcol.items(): ax.plot([], [], "o-", color=col, ms=5, label=t)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3); ax.set_title("At rest: no stage is critical in any organ\n(one line per cohort, 8 cohorts)")
ax = fig.add_subplot(gs[1, 2]); L.append((ax, "d")); k = 0; ticks = []
for acc, (t, lab_) in tis.items():
    if acc in rows and "contraction_healthy_minus_T2D" in rows[acc]:
        r = rows[acc]; ratio = r.dispersion_T2D / r.dispersion_healthy; lo = (r.dispersion_healthy - r.contraction_CI_hi) / r.dispersion_healthy; hi = (r.dispersion_healthy - r.contraction_CI_lo) / r.dispersion_healthy
        ax.errorbar(ratio, k, xerr=[[max(ratio - max(lo, 0.05), 0)], [max(hi - ratio, 0)]], fmt="o", color=tcol[t], ms=6, capsize=3, elinewidth=1.2); ticks.append(f"{t}, {lab_}"); k += 1
ax.set_yticks(range(k)); ax.set_yticklabels(ticks); ax.axvline(1, color="k", lw=0.9, ls="--"); ax.set_xlabel("dispersion, T2D / healthy (95% CI)"); ax.set_xscale("log"); ax.set_ylim(-0.6, k - 0.4); ax.set_title("At rest: the occupied region neither\ncontracts nor expands consistently")
# --- fila 3: bajo insulina si hay senal (e) y mapa de evidencia (f) ---
ax = fig.add_subplot(gs[2, 0:2]); L.append((ax, "e"))
ent = [("muscle", "GSE22309_muscle_4h", "IS", "IR", "sensitive → resistant"), ("muscle", "GSE22309_muscle_4h", "IS", "T2D", "sensitive → T2D"), ("adipose", "Ryden2016_adipose_clamp", "non_obese", "obese_before", "non-obese → obese"), ("adipose", "GSE26637_adipose_clamp", "sensitive", "resistant", "sensitive → resistant")]
for i_, (t, ds, a_, b_, lab_) in enumerate(ent):
    mr = g(ds, b_, "magnitude") / g(ds, a_, "magnitude"); ca, cb = g(ds, a_, "coherence_loo"), g(ds, b_, "coherence_loo")
    ax.scatter(mr, i_ - 0.17, color=tcol[t], s=70, marker="s", zorder=3); ax.scatter(cb - ca + 1, i_ + 0.17, color=tcol[t], s=70, marker="o", facecolors="none", lw=1.8, zorder=3)
    ax.text(-0.02, i_, f"{t}: {lab_}\n(n = {int(g(ds, a_, 'n'))}/{int(g(ds, b_, 'n'))})", va="center", ha="right", fontsize=8, transform=ax.get_yaxis_transform())
ax.axvline(1, color="k", lw=0.9, ls="--"); ax.set_yticks([]); ax.set_xlim(0.3, 1.75); ax.set_ylim(-0.7, len(ent) - 0.3); ax.invert_yaxis(); ax.spines["left"].set_visible(False)
ax.scatter([], [], color="k", marker="s", s=70, label="magnitude ratio (resistant / sensitive)"); ax.scatter([], [], color="k", marker="o", facecolors="none", lw=1.8, s=70, label="1 + Δ coherence (resistant − sensitive)")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2); ax.set_xlabel("relative to the sensitive group")
ax.set_title("During insulin there is signal, and it differs by organ: adipose loses magnitude and keeps\ndirection; muscle keeps magnitude and loses direction")
ax = fig.add_subplot(gs[2, 2]); L.append((ax, "f")); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
for i_, (org, rest, resp, fails) in enumerate([("islet", "11 genes (n = 85)", "ex vivo, coherence 0.4–0.9", "not comparable"), ("adipose", "none", "magnitude ↓", "amplitude"), ("muscle", "none", "coherence ↓, clock uncoupled", "coordination")]):
    y = 0.88 - 0.33 * i_
    if org == "muscle": ax.add_patch(plt.Rectangle((0, y - 0.26), 1, 0.35, fc="#eef3fb", ec="none", zorder=0))
    ax.text(0.02, y, org, fontsize=9, color=tcol[org], fontweight="bold", va="center")
    ax.text(0.02, y - 0.085, f"at rest: {rest}", fontsize=7, va="center"); ax.text(0.02, y - 0.155, f"under perturbation: {resp}", fontsize=7, va="center")
    if fails != "—": ax.text(0.02, y - 0.225, f"what differs: {fails}", fontsize=7, va="center", style="italic")
ax.set_title("Muscle is where the response\nfragments and where paired data\nexist in all three states")
# --- fila 4: la lectura canonica (expresion diferencial en reposo) ---
if de_sum is not None:
    ax = fig.add_subplot(gs[3, 0]); L.append((ax, "g")); lab2 = [f"{r.tissue}\n{r.acc}" for r in de_sum.itertuples()]
    ax.scatter(range(len(de_sum)), de_sum.n_fdr10, s=80, color=[tcol[t] for t in de_sum.tissue], zorder=3); ax.vlines(range(len(de_sum)), 1, de_sum.n_fdr10.clip(lower=1), color=[tcol[t] for t in de_sum.tissue], lw=2)
    ax.set_yscale("symlog", linthresh=1); ax.set_ylim(0, 6e4); ax.set_xticks(range(len(de_sum))); ax.set_xticklabels(lab2, fontsize=7, rotation=30, ha="right"); ax.set_ylabel("genes at FDR < 0.1"); ax.set_xlim(-0.65, len(de_sum) - 0.35)
    for i_, r in enumerate(de_sum.itertuples()): ax.text(i_, max(r.n_fdr10, 1) * 2.0, f"{r.n_fdr10:,}", ha="center", fontsize=7.5)
    ax.set_title("Classical differential expression at rest")
    ax = fig.add_subplot(gs[3, 1:]); L.append((ax, "h"))
    piv = de_set.pivot_table(index="set", columns="acc", values="z"); pp_ = de_set.pivot_table(index="set", columns="acc", values="p")
    order = [c for c in ["GSE164416", "GSE50244", "GSE159984", "GSE25462", "METSIM"] if c in piv.columns]; piv = piv[order]; pp_ = pp_[order]
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{de_sum.set_index('acc').tissue[c]}\n{c}" for c in order], fontsize=7.5); ax.tick_params(axis="x", pad=2); ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=7.5)
    for i_ in range(piv.shape[0]):
        for j_ in range(piv.shape[1]):
            if pp_.values[i_, j_] < 0.05: ax.text(j_, i_, "*", ha="center", va="center", fontsize=10, color="w" if abs(piv.values[i_, j_]) > 2.5 else "k")
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("z of mean |t| vs random sets")
    ax.set_title("Set enrichment of the same contrasts")
save(fig, "Fig2_tissue_search", L)


# ================= ED Fig 2: la respuesta del islote ex vivo =================
isl = _opt(f"{RES}/de/islet_response_geometry.tsv")
if isl is not None:
    fig = plt.figure(figsize=(W * 0.62, 3.2)); gs = fig.add_gridspec(1, 2, left=0.14, right=0.97, top=0.9, bottom=0.26, wspace=0.55); L = []
    ax = fig.add_subplot(gs[0, 0]); L.append((ax, "a")); colr = {"palmitate": C["IR"], "high glucose": "#6b46c1", "palmitate + glucose": C["T2D"]}
    for i_, r in enumerate(isl.itertuples()):
        ax.scatter(r.coherence_loo, i_, color=colr.get(r.stimulus, C["grey"]), s=70, zorder=3); ax.hlines(i_, 0, r.coherence_loo, color=colr.get(r.stimulus, C["grey"]), lw=2)
        ax.text(r.coherence_loo + 0.03, i_, f"n = {int(r.n)}", va="center", fontsize=8)
    ax.set_yticks(range(len(isl))); ax.set_yticklabels([f"{r.stimulus}\n{r.phase}" for r in isl.itertuples()], fontsize=7.5); ax.set_xlim(0, 1.15); ax.set_xlabel("coherence (LOO)"); ax.invert_yaxis()
    ax.axvline(0, color="k", lw=0.8); ax.set_title("Islet, ex vivo perturbation")
    ax = fig.add_subplot(gs[0, 1]); L.append((ax, "b"))
    for i_, r in enumerate(isl.itertuples()): ax.scatter(r.magnitude, r.coherence_loo, color=colr.get(r.stimulus, C["grey"]), s=70, zorder=3)
    for k_, (lab_, col) in enumerate(colr.items()): ax.scatter([], [], color=col, s=70, label=lab_)
    ax.set_xlabel("|displacement| (group mean)"); ax.set_ylabel("coherence (LOO)"); ax.set_ylim(0, 1.05); ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=1)
    save(fig, "ED_Fig2_islet_response", L)

# ================= Fig 3: el basal es ciego =================
fig = plt.figure(figsize=(W, 6.2)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], width_ratios=[1.15, 1, 1], left=0.105, right=0.97, top=0.94, bottom=0.08, hspace=0.6, wspace=0.62); L = []
ax = fig.add_subplot(gs[0, 0]); L.append((ax, "a")); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.text(-0.08, 0.9, "GSE182120", fontsize=9.5, fontweight="bold", va="top")
for k_, t in enumerate(["49 individuals, two laboratories", "24 normal glucose tolerance", "25 type 2 diabetes", "clamp M-value measured in all", "resting biopsy, no perturbation"]): ax.text(-0.08, 0.72 - 0.105 * k_, "· " + t, fontsize=8, va="center")
ax.text(-0.08, 0.1, "model: expression ~ laboratory + age + BMI,\nthen ρ(residual, M-value) genome-wide", fontsize=8, va="center", bbox=dict(boxstyle="round,pad=0.4", fc="#f4f4f4", ec="none"))
ax.set_title("Design: does the resting muscle\ntranscriptome know how insulin-\nsensitive its owner is?")
ax = fig.add_subplot(gs[0, 1]); L.append((ax, "b")); p = np.sort(gw.p_M.values); exp = -np.log10((np.arange(len(p)) + 0.5) / len(p)); ax.scatter(exp, -np.log10(p), s=4, color=C["grey"], lw=0); lim = max(exp.max(), (-np.log10(p)).max()); ax.plot([0, lim], [0, lim], "k-", lw=0.9)
for gname, dx, dy in [("PPARGC1A", 8, 6), ("PDK4", 8, -10), ("TXNIP", -8, 2)]:
    if gname in gw.index:
        xx = -np.log10((gw.p_M.rank().loc[gname] - 0.5) / len(gw)); yy = -np.log10(gw.p_M.loc[gname]); ax.scatter(xx, yy, s=40, color=C["IS"], zorder=3)
        ax.annotate(gname, (xx, yy), fontsize=8, xytext=(dx, dy), textcoords="offset points", ha="right" if dx < 0 else "left")
ax.set_xlabel("expected −log10 P"); ax.set_ylabel("observed −log10 P"); ax.set_title(f"No gene tracks insulin sensitivity\n({len(gw):,} genes, 0 at FDR < 0.1;\nknown markers in blue)")
ax = fig.add_subplot(gs[0, 2]); L.append((ax, "c")); ax.scatter(gw.t_T2D_vs_NGT, -np.log10(gw.p_T2D), s=4, color=C["grey"], lw=0); ax.axhline(-np.log10(0.05), color="k", lw=0.8, ls="--")
for gname, dx, dy in [("PDK4", 8, 4), ("PPARGC1A", -8, -6)]:
    if gname in gw.index: ax.scatter(gw.t_T2D_vs_NGT.loc[gname], -np.log10(gw.p_T2D.loc[gname]), s=40, color=C["IS"], zorder=3); ax.annotate(gname, (gw.t_T2D_vs_NGT.loc[gname], -np.log10(gw.p_T2D.loc[gname])), fontsize=8, xytext=(dx, dy), textcoords="offset points", ha="right" if dx < 0 else "left")
ax.set_xlabel("t, T2D vs NGT (resting)"); ax.set_ylabel("−log10 P"); ax.set_title("Nor does the diagnosis itself\n(0 genes at FDR < 0.1)")
ax = fig.add_subplot(gs[1, 0:2]); L.append((ax, "d")); set_panel(ax, set_scores(gw.rho_M), "Not even as gene sets: only the oxidative programme moves at all\n(grey: random sets of equal size, 95%)", "mean |ρ| with M-value")
ax = fig.add_subplot(gs[1, 2]); L.append((ax, "e")); cl = [x for x in SETS["clock output"] if x in gw.index]; ax.scatter(gw.loc[cl, "rho_M"], range(len(cl)), s=50, color=C["grey"], zorder=3); ax.hlines(range(len(cl)), 0, gw.loc[cl, "rho_M"], color=C["grey"], lw=1.5)
ax.set_yticks(range(len(cl))); ax.set_yticklabels(cl); ax.axvline(0, color="k", lw=0.9); ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.6, len(cl) - 0.4); ax.set_xlabel("ρ with M-value"); ax.set_title("The clock-output genes, central to\nwhat follows, are flat at rest\n(all P > 0.2)")
save(fig, "Fig3_resting_blind", L)

# ================= Fig 4: respuesta sana =================
fig = plt.figure(figsize=(W, 8.6)); gs = fig.add_gridspec(3, 2, width_ratios=[1, 1], height_ratios=[1, 1, 1], left=0.105, right=0.97, top=0.95, bottom=0.05, hspace=0.75, wspace=0.55); L = []
axsub = gs[0, :].subgridspec(1, 3, wspace=0.55)
ax = fig.add_subplot(axsub[0]); L.append((ax, "a")); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.add_patch(Circle((2, 4.0), 0.3, fc="k")); ax.text(2, 3.0, "basal", ha="center", fontsize=8)
for ang, col in [(0.2, C["IS"]), (0.0, C["IS"]), (-0.25, C["IS"]), (1.2, C["IR"]), (-1.4, C["IR"])]: ax.add_patch(FancyArrowPatch((2, 4.0), (2 + 5 * np.cos(ang), 4.0 + 2.8 * np.sin(ang)), arrowstyle="->", color=col, lw=0.9, mutation_scale=8))
ax.text(0, 10.2, "displacement = insulin − basal, per person\ncoherence = mean cosine of each displacement\nto the leave-one-out group mean", ha="left", va="top", fontsize=7.5); ax.text(0, -0.3, "blue: shared direction (coherence ≈ 1)\norange: idiosyncratic (≈ 0)", ha="left", va="bottom", fontsize=7.5); ax.set_title("The measure")
ax = fig.add_subplot(axsub[1]); L.append((ax, "b"))
for d in D["IS"]: ax.add_patch(FancyArrowPatch((0, 0), (d @ u1, d @ u2), arrowstyle="->", color=C["IS"], lw=0.7, alpha=0.75, mutation_scale=6))
Lm = 1.1 * np.abs(D["IS"] @ np.column_stack([u1, u2])).max(); ax.set_xlim(-Lm, Lm); ax.set_ylim(-Lm, Lm); ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4); ax.set_xlabel("along group direction"); ax.set_ylabel("orthogonal"); ax.set_title(f"Insulin-sensitive, 4-h clamp\n(n = {len(D['IS'])}), coherence {R.coherence_loo(D['IS']):.2f}")
ax = fig.add_subplot(axsub[2]); L.append((ax, "c")); tc = [("GSE9105, 30 min", "GSE9105_muscle_healthy", "30min", 0.5, 12), ("GSE231509, meal 1 h", "GSE231509_muscle_meal_1h", "H", 1, 7), ("GSE7146, 2 h", "GSE7146_muscle_healthy_2h", "healthy", 2, 6), ("GSE9105, 4 h", "GSE9105_muscle_healthy", "240min", 4, 12), ("GSE22309, 4 h", "GSE22309_muscle_4h", "IS", 4.2, 20)]
for lab_, ds, gr, t, n in tc: ax.scatter(t, g(ds, gr, "coherence_loo"), s=8 + n, color=C["IS"], zorder=3, lw=0)
ax.plot([0.5, 1, 2, 4.1], [g("GSE9105_muscle_healthy", "30min", "coherence_loo"), g("GSE231509_muscle_meal_1h", "H", "coherence_loo"), g("GSE7146_muscle_healthy_2h", "healthy", "coherence_loo"), (g("GSE9105_muscle_healthy", "240min", "coherence_loo") + g("GSE22309_muscle_4h", "IS", "coherence_loo")) / 2], color=C["IS"], lw=0.6, ls=":")
offs = {"GSE9105, 30 min": (2, -14), "GSE231509, meal 1 h": (4, 6), "GSE7146, 2 h": (6, -4), "GSE9105, 4 h": (6, -8), "GSE22309, 4 h": (6, 4)}
for lab_, ds, gr, t, n in tc: ax.annotate(f"{lab_} (n = {n})", (t, g(ds, gr, "coherence_loo")), fontsize=7.5, xytext=offs[lab_], textcoords="offset points")
ax.set_xlabel("hours after stimulus"); ax.set_ylabel("coherence (LOO)"); ax.set_ylim(0, 1.05); ax.set_xlim(-0.2, 7.5); ax.set_title("Coordination is assembled\nover hours (healthy,\nfour independent datasets)")
ax = fig.add_subplot(gs[1, 0]); L.append((ax, "d")); common = g30.index.intersection(g240.index); ax.scatter(g30.loc[common, "t"], g240.loc[common, "t"], s=2, color=C["null"], lw=0); pg = [x for x in prog.index if x in common]; ax.scatter(g30.loc[pg, "t"], g240.loc[pg, "t"], s=9, color=C["IS"], lw=0, label="55-gene replicated\nprogramme")
ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4); ax.set_xlabel("paired t, 30 min"); ax.set_ylabel("paired t, 4 h"); ax.legend(frameon=False, loc="upper left"); ax.set_title("The programme is a 4-h event\n(GSE9105, 12 healthy men)")
ax = fig.add_subplot(gs[1, 1]); L.append((ax, "e")); set_panel(ax, set_scores(gresp["t"].dropna()), "Gene sets in the healthy 4-h response (GSE22309 IS)\n(grey: random sets of equal size, 95%)", "mean |paired t|")
ax = fig.add_subplot(gs[2, 0]); L.append((ax, "f")); cl = [x for x in ["DBP", "TEF", "HLF", "PER2", "NR1D2", "BHLHE40"] if x in g240.index and x in gresp.index]; y = np.arange(len(cl))
for j, (src, mk, col, lab_) in enumerate([(g30, "v", C["null"], "GSE9105, 30 min"), (g240, "o", "#7fa7d8", "GSE9105, 4 h"), (gresp, "s", C["IS"], "GSE22309 IS, 4 h")]): ax.scatter([src.loc[x, "t"] for x in cl], y + (j - 1) * 0.25, marker=mk, s=16, color=col, label=lab_, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(cl); ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("paired t (insulin − basal)"); ax.set_title("Clock output resets at 4 h\nin two cohorts"); ax.legend(frameon=False, loc="lower right", fontsize=7.5)
ax = fig.add_subplot(gs[2, 1]); L.append((ax, "g")); ieg = [x for x in SETS["immediate-early TFs"] if x in g30.index]; clo = [x for x in SETS["clock output"] if x in g30.index]
vals = np.array([[g30.loc[ieg, "t"].mean(), g240.loc[ieg, "t"].mean(), gresp.loc[[x for x in ieg if x in gresp.index], "t"].mean()], [g30.loc[clo, "t"].mean(), g240.loc[clo, "t"].mean(), gresp.loc[[x for x in clo if x in gresp.index], "t"].mean()]])
im = ax.imshow(vals, cmap="RdBu_r", vmin=-5, vmax=5, aspect="auto"); ax.set_xticks(range(3)); ax.set_xticklabels(["30 min\n(GSE9105)", "4 h\n(GSE9105)", "4 h\n(GSE22309)"], fontsize=7.5); ax.set_yticks([0, 1]); ax.set_yticklabels(["IEG", "clock"])
for i in range(2):
    for j in range(3): ax.text(j, i, f"{vals[i, j]:.1f}", ha="center", va="center", fontsize=6, color="w" if abs(vals[i, j]) > 3 else "k")
cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04); cb.set_label("mean paired t"); ax.set_title("Check: IEG (repeat-biopsy\nprone) respond at any time;\nclock output only at 4 h")


save(fig, "Fig4_healthy_coordinated", L)

# ================= Fig 5: IR fragmenta (estilo: 3 paneles por fila, grandes) =================
fig = plt.figure(figsize=(W, 12.2)); gs = fig.add_gridspec(5, 3, height_ratios=[1, 1.05, 1.05, 1, 1], left=0.08, right=0.98, top=0.98, bottom=0.035, hspace=0.62, wspace=0.5); L = []
# a magnitud
ax = fig.add_subplot(gs[0, 0]); L.append((ax, "a")); dots(ax, {k: np.linalg.norm(D[k], axis=1) for k in D}, C, "|displacement| per person", log=True); ax.set_title("Insulin-resistant muscle responds\nas strongly as healthy muscle")
# b coherencia con nulo
ax = fig.add_subplot(gs[0, 1]); L.append((ax, "b")); nIR = null_coh("IS", "IR"); nT = null_coh("IS", "T2D"); parts = ax.violinplot([nIR, nT], positions=[0, 1], showextrema=False, widths=0.75)
for pc in parts["bodies"]: pc.set_facecolor(C["null"]); pc.set_alpha(1)
oIR = R.coherence_loo(D["IS"]) - R.coherence_loo(D["IR"]); oT = R.coherence_loo(D["IS"]) - R.coherence_loo(D["T2D"]); ax.scatter([0, 1], [oIR, oT], color=[C["IR"], C["T2D"]], zorder=3, s=90)
pIR = tests[(tests.subset == "full") & (tests.test == "coherence IS - IR")].p.iloc[0]; pT = tests[(tests.subset == "full") & (tests.test == "coherence IS - T2D")].p.iloc[0]
ax.text(0, oIR + 0.05, ptxt(pIR), ha="center"); ax.text(1, oT + 0.05, ptxt(pT), ha="center"); ax.text(0.5, -0.45, "label\npermutation", ha="center", color=C["grey"], fontsize=8)
ax.set_xticks([0, 1]); ax.set_xticklabels(["IS − IR", "IS − T2D"]); ax.set_ylabel("loss of coherence (LOO)"); ax.axhline(0, color="k", lw=0.8); ax.set_xlim(-0.7, 1.7); ax.set_title("…but no longer in\na shared direction")
# c mismo lote
ax = fig.add_subplot(gs[0, 2]); L.append((ax, "c"))
for i, k in enumerate(["IS", "IR", "T2D"]):
    full = g("GSE22309_muscle_4h", k, "coherence_loo"); same = g("GSE22309_muscle_4h_sameRun", k, "coherence_loo"); nsame = int(g("GSE22309_muscle_4h_sameRun", k, "n"))
    ax.plot([i - 0.18, i + 0.18], [full, same], color=C[k], lw=1.5); ax.scatter(i - 0.18, full, s=70, facecolors="white", edgecolors=C[k], lw=1.8, zorder=3); ax.scatter(i + 0.18, same, s=70, color=C[k], zorder=3); ax.text(i + 0.18, same - 0.11, f"n = {nsame}", ha="center", fontsize=8)
ax.text(0.25, 0.93, "open: all pairs   filled: same-batch pairs", fontsize=8, color=C["grey"], transform=ax.transAxes); ax.set_xticks(range(3)); ax.set_xticklabels(["IS", "IR", "T2D"]); ax.set_ylim(0, 1.05); ax.set_xlim(-0.6, 2.6); ax.set_ylabel("coherence (LOO)"); ax.set_title("The effect survives restriction\nto same-batch biopsy pairs")
# d rosas (3 polares en una celda)
cs = {}
for k in ["IS", "IR", "T2D"]: cs[k] = [R.cosine(D[k][j], np.delete(D[k], j, 0).mean(0)) for j in range(len(D[k]))] if k == "IS" else [R.cosine(d, mIS) for d in D[k]]
cosT = tests[(tests.subset == "full") & (tests.test == "cos(dir IS, dir T2D)")]; sub = gs[1, 0:2].subgridspec(1, 3, wspace=0.3)
for k, key in enumerate(["IS", "IR", "T2D"]):
    axp = fig.add_subplot(sub[k], projection="polar"); R.rose(axp, cs[key], C[key]); axp.set_xticklabels(["healthy\ndirection", "", "opposite", ""] if k == 0 else ["", "", "", ""], fontsize=8); axp.set_title(f"{key}  (n = {len(cs[key])})", color=C[key], loc="center", pad=8, fontsize=9.5)
    if k == 0: L.append((axp, "d")); _rowtext(f"Each person's response direction, healthy direction at north: healthy people point together,\ninsulin-resistant people point everywhere (group cosine IS–T2D {cosT.observed.iloc[0]:.2f}, null {cosT.null_mean.iloc[0]:.2f}, {ptxt(cosT.p.iloc[0])})", fontsize=9.5, va="bottom")
# e replicacion
ax = fig.add_subplot(gs[1, 2]); L.append((ax, "e")); rep = [("IS, GSE22309", "GSE22309_muscle_4h", "IS", C["IS"]), ("healthy 4 h, GSE9105", "GSE9105_muscle_healthy", "240min", C["IS"]), ("IR, GSE22309", "GSE22309_muscle_4h", "IR", C["IR"]), ("prediabetes, GSE157988\nplacebo arm", "GSE157988_muscle_prediabetes_clamp", "Placebo_before", C["IR"]), ("prediabetes, GSE157988\nNMN arm", "GSE157988_muscle_prediabetes_clamp", "NMN_before", C["IR"]), ("T2D, GSE22309", "GSE22309_muscle_4h", "T2D", C["T2D"])]
for i, (lab_, ds, gr, col) in enumerate(rep): v = g(ds, gr, "coherence_loo"); ax.hlines(i, 0, v, color=col, lw=2.2); ax.scatter(v, i, color=col, s=70, zorder=3); ax.text(max(v, 0) + 0.025, i, f"n = {int(g(ds, gr, 'n'))}", va="center", ha="left", fontsize=8)
ax.set_yticks(range(len(rep))); ax.set_yticklabels([r[0] for r in rep], fontsize=8); ax.yaxis.tick_right(); ax.spines["left"].set_visible(False); ax.spines["right"].set_visible(True); ax.tick_params(axis="y", length=0); ax.set_xlim(-0.05, 1.25); ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("coherence (LOO)"); ax.invert_yaxis(); ax.set_title("Low coherence replicates in\nprediabetes (independent lab)")
genes = [x for x in ["DBP", "TEF", "HLF", "PER2", "NR1D2", "BHLHE40"] if x in clk.gene.values]; th = np.linspace(0, 2 * np.pi, len(genes), endpoint=False); vmax = 0.55
for k, key in enumerate(["IS", "IR", "T2D"]):
    axp = fig.add_subplot(gs[2, k], projection="polar"); vals = np.array([clk[clk.gene == gn].resp_IS.iloc[0] if key == "IS" else clk[(clk.gene == gn) & (clk.vs == key)].resp_other.iloc[0] for gn in genes])
    light = {"IS": "#9ec9d8", "IR": "#f0b48e", "T2D": "#d9a39f"}[key]
    axp.bar(th, np.minimum(np.abs(vals), vmax), width=2 * np.pi / len(genes) * 0.75, bottom=0, color=[C[key] if v < 0 else light for v in vals], edgecolor="white", lw=0.5)
    axp.set_theta_zero_location("N"); axp.set_theta_direction(-1); axp.set_xticks(th); axp.set_xticklabels(genes, fontsize=8); axp.set_ylim(0, vmax); axp.set_yticks([0.25, 0.5]); axp.set_yticklabels(["0.25", "0.5"], fontsize=7); axp.spines["polar"].set_linewidth(0.6); axp.tick_params(axis="x", pad=6); axp.set_rlabel_position(22)
    axp.text(-0.12, 1.16, key, transform=axp.transAxes, ha="left", color=C[key], fontsize=9.5, fontweight="bold")
    for i, gn in enumerate(genes):
        r = clk[(clk.gene == gn) & (clk.vs == "IR")]; r2 = clk[(clk.gene == gn) & (clk.vs == "T2D")]
        if key == "IR" and len(r) and r.fdr.iloc[0] < 0.05: axp.text(th[i], vmax * 0.85, "*", ha="center", va="center", fontsize=12)
        if key == "T2D" and len(r2) and r2.fdr.iloc[0] < 0.1: axp.text(th[i], vmax * 0.85, "†", ha="center", va="center", fontsize=9)
    if k == 0: L.append((axp, "f")); _rowtext("Insulin resets the clock output in healthy muscle and no longer does in insulin resistance: |response| per gene\n(dark: repressed, light: induced). * FDR < 0.05 IR vs IS; † FDR < 0.1 T2D vs IS", fontsize=9.5, va="bottom")
ax = fig.add_subplot(gs[3, 0]); L.append((ax, "g")); mg = [x for x in ["DBP", "TEF", "PER2", "BHLHE40", "TXNIP"] if x in myo.index]; y = np.arange(len(mg)); ax.scatter(myo.loc[mg, "t_NGT"], y - 0.17, color=C["IS"], s=60, marker="o", zorder=3); ax.scatter(myo.loc[mg, "t_T2D"], y + 0.17, color=C["T2D"], s=60, marker="^", zorder=3)
ax.scatter([], [], color=C["IS"], marker="o", s=60, label="NGT (n = 7)"); ax.scatter([], [], color=C["T2D"], marker="^", s=60, label="T2D (n = 5)"); ax.set_ylim(-0.7, len(mg) - 0.3); ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2)
ax.set_yticks(y); ax.set_yticklabels(mg); ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("t, chronic insulin vs control"); ax.set_title("In cells (GSE182117): DBP responds\nin NGT, BHLHE40 in both")
ax = fig.add_subplot(gs[3, 1]); L.append((ax, "h")); ag = [x for x in ["DBP", "TEF", "HLF", "PER2", "PER3", "NR1D2"] if x in myo.index]; y = np.arange(len(ag))
ax.errorbar(myo.loc[ag, "amp_NGT"], y - 0.17, xerr=myo.loc[ag, "amp_NGT_sd"], fmt="o", color=C["IS"], ms=6, capsize=3, elinewidth=1.2, label="NGT"); ax.errorbar(myo.loc[ag, "amp_T2D"], y + 0.17, xerr=myo.loc[ag, "amp_T2D_sd"], fmt="^", color=C["T2D"], ms=6, capsize=3, elinewidth=1.2, label="T2D")
ax.set_yticks(y); ax.set_yticklabels(ag); ax.set_xlabel("24-h amplitude (cosinor), mean ± s.d."); ax.set_xlim(0, None); ax.set_title("Clock-output amplitude is\n~30% lower in T2D myotubes"); ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2)
sub = gs[4, 0:2].subgridspec(1, 2, wspace=0.55)
def fate(sfx):
    same = np.sign(prog["mean"]) == np.sign(prog[f"mean_{sfx}"]); strong = prog[f"t_{sfx}"].abs() > 2; return {"kept": (same & strong).mean(), "lost": (~strong).mean(), "inverted": (~same & strong).mean()}
for k, sfx in enumerate(["IR", "T2D"]):
    axa = fig.add_subplot(sub[k]); R.alluvial(axa, fate(sfx), "healthy\nprogramme", ["kept", "lost", "inverted"], {"kept": C["IS"], "lost": C["null"], "inverted": C["T2D"]}); axa.text(0.4, 1.06, sfx, transform=axa.transAxes, ha="center", color=C[sfx], fontsize=9.5, fontweight="bold")
    if k == 0: L.append((axa, "i"))

if exo is not None:
    axj = fig.add_subplot(gs[4, 2]); L.append((axj, "j"))
    ins = [("insulin, 4 h", g("GSE22309_muscle_4h", "IS", "coherence_loo"), g("GSE22309_muscle_4h", "T2D", "coherence_loo"))]
    ex = exo[(exo.tissue == "muscle") & (exo.timepoint == "recovery")]
    ins.append(("exercise, 3 h\n(recovery)", float(ex[ex.group == "NGT"].coherence.iloc[0]), float(ex[ex.group == "T2D"].coherence.iloc[0])))
    for i_, (lab_, h, d_) in enumerate(ins):
        axj.plot([i_ - 0.17, i_ + 0.17], [h, d_], color="k", lw=1, zorder=2)
        axj.scatter(i_ - 0.17, h, color=C["IS"], s=80, zorder=3); axj.scatter(i_ + 0.17, d_, color=C["T2D"], s=80, zorder=3)
    axj.set_xticks(range(len(ins))); axj.set_xticklabels([x[0] for x in ins]); axj.set_ylim(0, 1); axj.set_xlim(-0.5, len(ins) - 0.5); axj.set_ylabel("coherence (LOO)")
    axj.scatter([], [], color=C["IS"], s=60, label="healthy / NGT"); axj.scatter([], [], color=C["T2D"], s=60, label="T2D")
    axj.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2)
    if exo_t is not None:
        tt_ = exo_t[(exo_t.tissue == "muscle") & (exo_t.timepoint == "recovery")]
        if len(tt_): axj.text(1, 0.93, f"P = {tt_.p_coherence.iloc[0]:.2f}", ha="center", fontsize=8)
    pIS_T2D = tests[(tests.subset == "full") & (tests.test == "coherence IS - T2D")].p.iloc[0]
    axj.text(0, 0.93, ptxt(pIS_T2D), ha="center", fontsize=8)
save(fig, "Fig5_IR_fragments", L)


# ================= ED Fig 3: sesgo, confusion por lote y potencia =================
if pwr is not None:
    fig = plt.figure(figsize=(W * 0.95, 3.4)); gs = fig.add_gridspec(1, 3, left=0.08, right=0.98, top=0.9, bottom=0.26, wspace=0.5); L = []
    n0 = pwr[pwr.analysis == "null coherence vs n"]
    ax = fig.add_subplot(gs[0, 0]); L.append((ax, "a"))
    ax.plot(n0.n, n0.naive_mean, "o-", color=C["T2D"], ms=4, label="naive, mean"); ax.plot(n0.n, n0.naive_p95, "--", color=C["T2D"], lw=1, label="naive, 95th pct")
    ax.plot(n0.n, n0.loo_mean, "o-", color=C["IS"], ms=4, label="leave-one-out, mean"); ax.plot(n0.n, n0.loo_p95, "--", color=C["IS"], lw=1, label="LOO, 95th pct")
    ax.axhline(0, color="k", lw=0.8); ax.set_xlabel("individuals per group"); ax.set_ylabel("coherence when no direction is shared"); ax.legend(frameon=False, fontsize=7, loc="upper right")
    ps = pwr[pwr.analysis == "pseudo-group by batch"]; br = pwr[pwr.analysis == "batch restriction"]
    ax = fig.add_subplot(gs[0, 1]); L.append((ax, "b"))
    for i_, r in enumerate(ps.itertuples()): ax.scatter(i_, r.coherence, color=C["grey"], s=70, zorder=3); ax.text(i_, r.coherence + 0.05, f"n={int(r.n)}", ha="center", fontsize=7)
    k0 = len(ps)
    for j_, r in enumerate(br.itertuples()):
        ax.scatter(k0 + j_ - 0.15, r.coherence, facecolors="none", edgecolors=C[r.group], s=70, lw=1.6, zorder=3)
        ax.scatter(k0 + j_ + 0.15, r.coherence_same, color=C[r.group], s=70, zorder=3); ax.plot([k0 + j_ - 0.15, k0 + j_ + 0.15], [r.coherence, r.coherence_same], color=C[r.group], lw=1)
    ax.set_xticks(range(k0 + len(br))); ax.set_xticklabels(list(ps.group) + list(br.group), fontsize=7, rotation=30, ha="right")
    ax.set_ylim(0, 1.15); ax.set_ylabel("coherence (LOO)"); ax.axvline(k0 - 0.5, color="k", lw=0.6, ls=":")
    ax.text(0.02, 0.10, "groups defined by\nbatch alone", transform=ax.transAxes, fontsize=7, color=C["grey"], va="bottom"); ax.text(0.62, 0.10, "clinical groups\nopen: all pairs\nfilled: same batch", transform=ax.transAxes, fontsize=7, va="bottom")
    po = pwr[pwr.analysis == "power"]
    ax = fig.add_subplot(gs[0, 2]); L.append((ax, "c"))
    for d_, grp in po.groupby("delta_fraction_randomised"):
        ax.plot(grp.n, grp.power, "o-", ms=4, label=f"{int(d_*100)}% randomised")
    ax.axhline(0.8, color="k", lw=0.8, ls="--"); ax.set_xlabel("individuals per group"); ax.set_ylabel("power to detect the loss"); ax.set_ylim(0, 1); ax.legend(frameon=False, fontsize=7, loc="upper left")
    save(fig, "ED_Fig3_power_and_confounding", L)

# ================= Fig 6: adiposo e intervenciones =================
fig = plt.figure(figsize=(W, 6.0)); gs = fig.add_gridspec(2, 3, left=0.09, right=0.97, top=0.96, bottom=0.08, hspace=0.55, wspace=0.6); L = []
def mc(ax, rows_, labels, title):
    x = np.arange(len(rows_)); ax.plot(x, rows_.magnitude, color=C["grey"], lw=1.5, zorder=2); ax.scatter(x, rows_.magnitude, color=C["grey"], s=70, marker="s", zorder=3); ax.set_ylabel("|displacement| (group mean)", color=C["grey"])
    ax2 = ax.twinx(); ax2.plot(x, rows_.coherence_loo, color=C["IS"], lw=1.5, zorder=2); ax2.scatter(x, rows_.coherence_loo, color=C["IS"], s=70, zorder=3); ax2.set_ylim(0, 1); ax2.set_ylabel("coherence (LOO)", color=C["IS"]); ax2.spines["top"].set_visible(False)
    ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_xlim(-0.4, len(rows_) - 0.6); ax.set_ylim(0, max(rows_.magnitude) * 1.25); ax.set_title(title)
    for i_, n in enumerate(rows_.n): ax.text(i_, max(rows_.magnitude) * 0.06, f"n = {int(n)}", ha="center", fontsize=8, color=C["grey"])
ax = fig.add_subplot(gs[0, 0]); L.append((ax, "a")); ry = geom[geom.dataset == "Ryden2016_adipose_clamp"].set_index("group").loc[["non_obese", "obese_before", "obese_2y_after"]]; mc(ax, ry, ["non-\nobese", "obese", "obese,\n2 y after\nsurgery"], "Adipose tissue, clamp (Rydén 2016):\nsurgery restores the magnitude of the\nresponse, not its coherence")
ax = fig.add_subplot(gs[0, 1]); L.append((ax, "b")); g26 = geom[geom.dataset == "GSE26637_adipose_clamp"].set_index("group").loc[["sensitive", "resistant"]]; mc(ax, g26, ["insulin-\nsensitive", "insulin-\nresistant"], "Adipose tissue, clamp (GSE26637):\nhalf the magnitude, same direction\n(independent cohort)")
ax = fig.add_subplot(gs[0, 2]); L.append((ax, "c")); ent = [("adipose, Rydén", "Ryden2016_adipose_clamp", "non_obese", "obese_before"), ("adipose, GSE26637", "GSE26637_adipose_clamp", "sensitive", "resistant"), ("muscle, GSE22309", "GSE22309_muscle_4h", "IS", "IR")]
for i_, (lab_, ds, a_, b_) in enumerate(ent):
    mr = g(ds, b_, "magnitude") / g(ds, a_, "magnitude"); col = C["IR"] if "adipose" in lab_ else C["IS"]
    ax.scatter(mr, i_, color=col, s=80, marker="s", zorder=3); ax.hlines(i_, 1, mr, color=col, lw=2)
ax.axvline(1, color="k", lw=0.9, ls="--"); ax.set_yticks(range(len(ent))); ax.set_yticklabels([e[0] for e in ent]); ax.yaxis.tick_right(); ax.spines["left"].set_visible(False); ax.spines["right"].set_visible(True); ax.tick_params(axis="y", length=0); ax.set_xlim(0.45, 1.3); ax.set_ylim(-0.6, len(ent) - 0.4); ax.invert_yaxis(); ax.set_xlabel("magnitude ratio (resistant / sensitive)"); ax.set_title("Two failure modes side by side:\nadipose halves its response,\nmuscle does not")
ax = fig.add_subplot(gs[1, 0]); L.append((ax, "d")); nm_ = geom[geom.dataset == "GSE157988_muscle_prediabetes_clamp"].set_index("group")
for j_, arm in enumerate(["Placebo", "NMN"]):
    b_, a_ = nm_.loc[f"{arm}_before", "coherence_loo"], nm_.loc[f"{arm}_after", "coherence_loo"]; col = C["IR"] if arm == "NMN" else C["grey"]
    ax.plot([j_ - 0.18, j_ + 0.18], [b_, a_], "-", color=col, lw=1.5); ax.scatter([j_ - 0.18, j_ + 0.18], [b_, a_], color=col, s=70, zorder=3)
    pv = tests[(tests.dataset == "GSE157988") & (tests.subset == arm)].p.iloc[0]; yb = max(b_, a_) + 0.14; ax.plot([j_ - 0.18, j_ - 0.18, j_ + 0.18, j_ + 0.18], [yb - 0.03, yb, yb, yb - 0.03], "k-", lw=1); ax.text(j_, yb + 0.02, ptxt(pv), ha="center", fontsize=8)
ax.set_xticks([0, 1]); ax.set_xticklabels(["placebo\npre → post", "NMN\npre → post"]); ax.set_ylim(-0.15, 1); ax.set_xlim(-0.6, 1.6); ax.set_ylabel("coherence (LOO)"); ax.set_title("Prediabetic muscle, 10 weeks:\nNMN improves clamp-measured\nsensitivity but not coherence")
ax = fig.add_subplot(gs[1, 1]); L.append((ax, "e")); n = np.arange(5, 41); det = 2.8 * np.sqrt(2 * 0.3**2 / n); ax.plot(n, det, color="k")
ax.axvline(11, color=C["IR"], ls="--", lw=1); ax.axvline(23, color=C["IS"], ls="--", lw=1); ax.text(12, 0.86, "NMN trial\nn = 11", fontsize=8, color=C["IR"], va="top"); ax.text(25, 0.60, "surgery\ncohort, n = 23", fontsize=8, color=C["IS"], va="top")
ax.set_xlabel("n per arm"); ax.set_ylabel("detectable Δ coherence\n(80% power, s.d. 0.3)"); ax.set_ylim(0, 0.9); ax.set_title("What the present cohorts could\ndetect: restorations below ~0.25\nwould be missed")
ax = fig.add_subplot(gs[1, 2]); L.append((ax, "f")); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.text(2.9, 9.9, "muscle", ha="center", fontsize=9.5, fontweight="bold", color=C["IS"], va="top"); ax.text(7.9, 9.9, "adipose", ha="center", fontsize=9.5, fontweight="bold", color=C["IR"], va="top")
ax.text(0.0, 7.0, "healthy", fontsize=8.5, rotation=90, va="center", ha="center"); ax.text(0.0, 2.8, "insulin-\nresistant", fontsize=8.5, rotation=90, va="center", ha="center")
for ang in [0.14, 0, -0.14]: ax.add_patch(FancyArrowPatch((1.3, 7.0), (1.3 + 3.2 * np.cos(ang), 7.0 + 2.2 * np.sin(ang)), arrowstyle="->", color=C["IS"], lw=1.4, mutation_scale=9))
for ang in [1.0, 0.2, -0.9]: ax.add_patch(FancyArrowPatch((1.3, 2.8), (1.3 + 3.2 * np.cos(ang), 2.8 + 2.2 * np.sin(ang)), arrowstyle="->", color=C["IS"], lw=1.4, mutation_scale=9))
for ang in [0.14, 0, -0.14]: ax.add_patch(FancyArrowPatch((6.3, 7.0), (6.3 + 3.2 * np.cos(ang), 7.0 + 2.2 * np.sin(ang)), arrowstyle="->", color=C["IR"], lw=1.4, mutation_scale=9))
for ang in [0.14, 0, -0.14]: ax.add_patch(FancyArrowPatch((6.3, 2.8), (6.3 + 1.5 * np.cos(ang), 2.8 + 1.0 * np.sin(ang)), arrowstyle="->", color=C["IR"], lw=1.4, mutation_scale=9))
ax.text(2.8, 0.4, "loses direction", ha="center", fontsize=8.5); ax.text(7.9, 0.4, "loses magnitude", ha="center", fontsize=8.5)
ax.set_title("Summary: each organ fails\nin its own way")
save(fig, "Fig6_adipose_interventions", L)
