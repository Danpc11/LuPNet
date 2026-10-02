"""Progressive loss of lung units and the load indices recorded along the way.

Load indices (healthy lung = 1). Each resistance/elastance load is taken from the
network itself and normalised by functional tissue, never defined as a ratio of
the other two loads, so the identities below are checks, not definitions.

  vascular   z_P = dP / dP0                        (mPAP - LAP)
             z_F = (Q / m_Q) / (Q0 / m_Q0)          flow per perfused capillary mass
             z_R = (R m_Q) / (R0 m_Q0)              R = equivalent network resistance
             z_P = z_F z_R  exactly
  airway     z_dP = dP_L / dP_L0                    transpulmonary driving pressure (V_T E_L)
             z_V = (V_T / FRC) / (V_T0 / FRC0)      tidal strain of the aerated lung
             z_E = (E_L FRC) / (E_L0 FRC0)          specific elastance
             z_dP = z_V z_E  exactly
"""
import numpy as np
from .lung import Lung
from .septal import make_lung
from .params import Params


def hazards(L: Lung, st):
    p = L.p
    s = st["air"]["strain"]
    nfrac = L.stiff_nbr_frac()
    x = np.clip(p.beta * (s / L.s_ref - 1.0), -50, 50)
    lam_f = np.where(L.state == 0, p.lam0 * np.exp(x) * (1 + p.kn * nfrac), 0.0)
    if p.parenchyma == "surfactant":
        h = L.state == 0
        # constant injury flux: per-unit rates scale as healthy_fraction^(-flux_gamma)
        fx = min(max(float(h.mean()), 1e-3) ** (-p.flux_gamma), 50.0)
        lam_f = np.where(h, fx * p.lam0_direct * np.exp(x) * (1 + p.kn * nfrac), 0.0)
        lam_f = lam_f + np.where(h & (L.rd == 2), p.lam_ci, 0.0) + np.where(h & (L.rd == 1), p.lam_ci_cyc, 0.0)
        if "events" in st["air"]:
            lam_f = lam_f + np.where(h, p.lam_ae * st["air"]["events"], 0.0)
        if p.epi_bio:
            # aberrant basaloid niche activates fibroblasts; fibroblast activation grows with local stiffness
            # (YAP/TAZ, PIEZO1/2 mechanosensing): Hill function of the stiff-neighbour fraction
            ab = L.aberrant.astype(float)
            ab_nb = np.array([(ab[i] + ab[n].sum()) / (1 + len(n)) for i, n in enumerate(L.nbrs)])
            yap = 1 + p.k_yap * nfrac ** 2 / (p.yap_K ** 2 + nfrac ** 2)
            lam_f = lam_f + np.where(h, p.k_fib_ab * ab_nb * yap, 0.0)
        lam_f = lam_f * p.antifib
    lam_h = np.where(L.state == 1, p.lam_hc, 0.0)
    return lam_f, lam_h


def surfactant_step(L, st, dt):
    """AT2 dysfunction: surfactant quality decays at a unit-specific rate that grows with tidal strain
    and with stiff or collapsed neighbours (Wu 2020 tension; Lutz 2015 sequence)."""
    p = L.p
    s = st["air"]["strain"]
    bad = ((L.state > 0) | (L.rd == 2)).astype(float)
    nb = np.array([bad[n].mean() for n in L.nbrs])
    fx = min(max(float((L.state == 0).mean()), 1e-3) ** (-p.flux_gamma), 50.0)
    if p.epi_bio:
        s = s * (1 + p.pleural_stress * (L.pleural_dist <= 1))     # stress concentration under the pleura
    rate = fx * p.lam_s * L.at2_rate * (1 + p.beta_s * np.maximum(s / L.s_ref - 1, 0) + p.kn_s * nb)
    if p.epi_bio:
        epi_step(L, s, dt, rate)
        rate = rate * (1 + p.g_sen * L.sen)                        # senescent AT2 cells renew surfactant poorly
    L.sq = np.where(L.state == 0, L.sq * np.exp(-rate * dt), L.sq)


def epi_step(L, s, dt, rate):
    """AT2 senescence (age / telomere background, stretch, SASP from senescent neighbours; Wu 2020, ERR
    2026 review) and conversion of senescent, injured, stretched epithelium into aberrant basaloid cells
    (arrested AT2 -> AT1 differentiation with sustained YAP/TAZ; Habermann 2020, Adams 2020, 2025-26 data)."""
    p = L.p
    h = L.state == 0
    sen_nb = np.array([L.sen[n].mean() for n in L.nbrs])
    dsen = p.k_sen0 * L.at2_rate + p.k_sen_strain * np.maximum(s / L.s_ref - 1, 0) + p.k_sasp * sen_nb
    L.sen = np.where(h, np.minimum(1.0, L.sen + dsen * dt * (1 - L.sen)), L.sen)
    lam_ab = p.k_ab * L.sen * (1 - L.sq) * np.maximum(s / L.s_ref, 0.0)
    rng = np.random.default_rng(int(1e6 * (L.sen.sum() + 1)) % (2 ** 32))
    L.aberrant = L.aberrant | (h & (rng.random(L.N) < 1 - np.exp(-lam_ab * dt)))


def record(L: Lung, st, t, ex=None):
    p, b = L.p, L.base_state
    vs, air, g = st["vs"], st["air"], st["gas"]
    vs0, air0 = b["vs"], b["air"]
    mQ = float(st["cf"].sum()); mQ0 = float(b["cf"].sum())
    Qn, Q0 = st["Q"], b["Q"]
    zF = (Qn / mQ) / (Q0 / mQ0) if mQ > 0 else np.inf
    zR = (vs["Req"] * mQ) / (vs0["Req"] * mQ0)
    zP = vs["dP"] / vs0["dP"]
    zV = (air["VT"] / air["FRC"]) / (air0["VT"] / air0["FRC"])
    zE = (air["EL"] * air["FRC"]) / (air0["EL"] * air0["FRC"])
    zdP = air["dPL"] / air0["dPL"]
    pv = L.pvv(vs)
    # decomposition: structural remodelling at baseline pressure vs distension of baseline structure
    P0 = b["vs"]["P"]
    r_rem_a = L.ru_a * (1 + p.alpha * P0["Pm_a"]); r_rem_v = L.ru_v * (1 + p.alpha * P0["Pm_v"])
    pv_rem = L.pvv(dict(r_a=r_rem_a, r_v=r_rem_v))
    r_dis_a = L.ru0[0] * (1 + p.alpha * vs["P"]["Pm_a"]); r_dis_v = L.ru0[1] * (1 + p.alpha * vs["P"]["Pm_v"])
    pv_dis = L.pvv(dict(r_a=r_dis_a, r_v=r_dis_v))
    Ql, VAi = vs["Ql"], air["VAi"]
    shunt = float(Ql[VAi < 1e-6 * VAi.max()].sum() / Ql.sum())
    alvD = float(VAi[Ql < 1e-9].sum() / VAi.sum())
    s = L.state
    lost = s > 0
    y = L.y
    base = y < 1 / 3; apex = y > 2 / 3
    sub = L.pleural_dist <= 2; core = L.pleural_dist >= 6
    frac = lambda m: float(lost[m].mean())
    r = dict(ae=0, death="", Dm=L._Dm / L.Dm0 if getattr(L, "Dm0", None) else np.nan,
             Vc=L._Vc / L.Vc0 if getattr(L, "Vc0", None) else np.nan, years=t * p.years_per_unit, FVC_mL=air["C"] / air0["C"] * p.FVC0_mL,
             events=float(np.mean(air["events"])) if "events" in air else 0.0, t=t, healthy=float((s == 0).mean()), fibrotic=float((s == 1).mean()),
             honeycomb=float((s == 2).mean()),
             FVC=air["C"] / air0["C"], FRC=air["FRC"] / air0["FRC"], DLCO=st["DL"] / L.DL0,
             PaO2=g["PaO2"], SaO2=g["SaO2"], PvO2=g["Pv"], gas_ok=g["ok"],
             mPAP=vs["PPA"], dP=vs["dP"], zP=zP, zF=zF, zR=zR, zdP=zdP, zV=zV, zE=zE,
             f=air["f"], VT=air["VT"], dPL=air["dPL"],
             strain_med_healthy=float(np.median(air["strain"][s == 0])) if (s == 0).any() else np.nan,
             shunt=shunt, alv_deadspace=alvD, hpv_mean=float((st["h"] * Ql).sum() / Ql.sum()),
             PVV_large=pv["large"] / L.PVV0["large"], PVV_small=pv["small"] / L.PVV0["small"],
             PVV_large_remodel=pv_rem["large"] / L.PVV0["large"], PVV_small_remodel=pv_rem["small"] / L.PVV0["small"],
             PVV_large_distend=pv_dis["large"] / L.PVV0["large"], PVV_small_distend=pv_dis["small"] / L.PVV0["small"],
             lost_base=frac(base), lost_apex=frac(apex), lost_sub=frac(sub), lost_core=frac(core),
             eq_it=st["it"], zerr=st["zerr"], converged=bool(st["zerr"] < 1e-3),
             tether_mean=float(np.mean(st["Pt"])) if np.ndim(st["Pt"]) else 0.0,
             PVR=vs["dP"] / Qn, Q=Qn)
    # field effect: healthy units next to scar vs healthy units with no stiff neighbour
    hl = s == 0
    adj = hl & (L.stiff_nbr_frac() > 0); far = hl & ~adj
    q0 = float(b["vs"]["Ql"].mean()); pc0 = float(b["vs"]["Pcap"].mean())
    mean = lambda x, m: float(x[m].mean()) if m.any() else np.nan
    r.update(Pcap_adj=mean(vs["Pcap"], adj) - pc0, Pcap_far=mean(vs["Pcap"], far) - pc0,
             Qunit_adj=mean(Ql, adj) / q0, Qunit_far=mean(Ql, far) / q0,
             PcO2_adj=mean(g["Pc"], adj), PcO2_far=mean(g["Pc"], far),
             n_adj=int(adj.sum()), n_far=int(far.sum()))
    # the same comparison within horizontal slices (removes the height confound)
    bins = np.minimum((L.y * 8).astype(int), 7); dq, dpo, w = 0.0, 0.0, 0.0
    for k in range(8):
        ma, mf = adj & (bins == k), far & (bins == k)
        if ma.any() and mf.any():
            wk = min(ma.sum(), mf.sum())
            dq += wk * (Ql[ma].mean() / Ql[mf].mean() - 1); dpo += wk * (g["Pc"][ma].mean() - g["Pc"][mf].mean()); w += wk
    r.update(Qunit_adj_vs_far_sameheight=dq / w if w else np.nan, PcO2_adj_minus_far_sameheight=dpo / w if w else np.nan)
    if p.perf_grav:
        zn = L.zones(vs)
        r.update(zone1=float((zn == 1).mean()), zone2=float((zn == 2).mean()), zone3=float((zn == 3).mean()))
    if p.epi_bio:
        r.update(sen_mean=float(L.sen[s == 0].mean()) if (s == 0).any() else np.nan,
                 aberrant=float(L.aberrant.mean()),
                 aberrant_sub=float(L.aberrant[L.pleural_dist <= 2].mean()), aberrant_core=float(L.aberrant[L.pleural_dist >= 6].mean()))
    if p.parenchyma == "surfactant":
        hl0 = L.state == 0
        r.update(closed=float(((L.rd == 2) & hl0).mean()), cyclic=float(((L.rd == 1) & hl0).mean()),
                 sq_mean=float(L.sq[hl0].mean()) if hl0.any() else np.nan,
                 closed_base=float(((L.rd == 2) & hl0)[L.y < 1 / 3].mean()),
                 closed_apex=float(((L.rd == 2) & hl0)[L.y > 2 / 3].mean()))
    # CO2 and the ventilatory limit: ventilation to perfused units clears CO2; the chemoreflex raises
    # alveolar ventilation to keep PaCO2 at 40 mmHg; the breathing power this needs is compared with rest
    perf = vs["Ql"] > 1e-6 * vs["Ql"].max()
    VAeff = float(air["VAi"][perf].sum())
    VCO2 = 0.8 * p.VO2
    r["PaCO2"] = 0.863 * VCO2 / max(VAeff, 1e-6)
    eff = VAeff / max(float(air["VAi"].sum()), 1e-9)
    VA_req = 0.863 * VCO2 / 40.0 / max(eff, 1e-3)
    r["VA_req"] = VA_req
    r["power_ratio"] = float(L._breathe(VA_req)["power"] / air0["power"]) if hasattr(L, "_breathe") else np.nan
    rv = st.get("rv")
    if rv is not None:
        r.update(SV=rv["SV"], Ees=rv["Ees"], EDV=rv["EDV"], RAP=rv["RAP"], coupling=rv["coupling"],
                 Qcap=rv["Qcap"], rv_limited=rv["limited"])
    if ex is not None:
        if ex.get("rv") is not None:
            r.update(Q_ex=ex["Q"], Qcap_ex=ex["rv"]["Qcap"], RAP_ex=ex["rv"]["RAP"],
                     coupling_ex=ex["rv"]["coupling"], rv_limited_ex=ex["rv"]["limited"])
        r.update(SaO2_ex=ex["gas"]["SaO2"], PaO2_ex=ex["gas"]["PaO2"], mPAP_ex=ex["vs"]["PPA"],
                 f_ex=ex["air"]["f"], transit_ex=float(np.median(ex["gas"]["t_transit"][ex["vs"]["Ql"] > 0])))
    return r


def exercise(L: Lung):
    p = L.p
    h, ru = L.h.copy(), (L.ru_a.copy(), L.ru_v.copy(), L.ru_ta.copy())
    x0, rd0, bud = L.x.copy(), L.rd.copy(), L._rem_budget
    L._rem_budget = None
    ex = L.equilibrate(Qtot=p.Q * p.ex_Q, VO2=p.VO2 * p.ex_VO2, VA=p.VA * p.ex_VA, remodel=False,
                       exercise=True)
    L.h = h; L.ru_a, L.ru_v, L.ru_ta = ru
    L.x, L.rd, L._rem_budget = x0, rd0, bud
    return ex


def run(p: Params, snapshots=(0.1, 0.3, 0.5), do_exercise=True, verbose=False):
    rng = np.random.default_rng(p.seed)
    L = make_lung(p)
    st = L.base_state
    rows = [record(L, st, 0.0, exercise(L) if do_exercise else None)]
    snaps = {}
    t = 0.0
    while rows[-1]["healthy"] > p.stop_healthy and t < p.t_max and rows[-1]["mPAP"] < p.stop_mPAP \
            and rows[-1]["Q"] > p.stop_Q and rows[-1]["gas_ok"]:
        lf, lh = hazards(L, st)
        tot = lf.sum() + lh.sum()
        dt = min(0.05, p.max_frac_step * L.N / max(tot, 1e-12))
        fail = (rng.random(L.N) < 1 - np.exp(-lf * dt))
        honey = (rng.random(L.N) < 1 - np.exp(-lh * dt))
        L.state[fail] = 1
        L.state[honey & ~fail] = 2
        if p.parenchyma == "surfactant":
            surfactant_step(L, st, dt)
        t += dt
        if p.vasc_A > 0:
            sv = p.vasc_A * (1.0 - np.exp(-(t * p.years_per_unit + p.t_vasc0) / p.vasc_tau))
            L.R_lesion = np.full(L.N, L._Rta0 * sv)
            if p.cap_frac > 0:
                L.cap_vasc = 1.0 / (1.0 + p.cap_frac * sv)
        if p.k_vasc > 0:
            L.R_lesion = np.full(L.N, L._Rta0 * (np.exp(p.k_vasc * (t * p.years_per_unit + p.t_vasc0)) - 1.0))
            if p.cap_frac > 0:
                L.cap_vasc = 1.0 / (1.0 + p.cap_frac * (np.exp(p.k_vasc * (t * p.years_per_unit + p.t_vasc0)) - 1.0))
        if p.epi_bio and p.senolytic_year >= 0 and not L._senolytic_done and t * p.years_per_unit >= p.senolytic_year:
            L.sen = L.sen * (1 - p.senolytic_frac); L._senolytic_done = True
        if p.tau_rem_weeks > 0:
            L._rem_budget = dt * p.years_per_unit * 52.0 / p.tau_rem_weeks
        st = L.equilibrate()
        rows.append(record(L, st, t, exercise(L) if do_exercise else None))
        if p.ae:
            fvc = rows[-1]["FVC"]
            rate = p.ae_rate0 + p.ae_rate_slope * max(0.0, 1 - fvc)          # per year
            if rng.random() < 1 - np.exp(-rate * dt * p.years_per_unit):
                sq0 = L.sq.copy()
                hit = np.clip(p.ae_hit * rng.uniform(0.5, 1.5, L.N), 0, 0.95)
                L.sq = np.where(L.state == 0, L.sq * (1 - hit), L.sq)
                # diffuse alveolar damage: a random fraction of open units floods (perfused, unventilated)
                ff = rng.uniform(p.ae_flood_lo, p.ae_flood_hi)
                L.flood = (L.state == 0) & (L.rd != 2) & (rng.random(L.N) < ff)
                bud = L._rem_budget; L._rem_budget = None
                hmax, fio2 = p.hpv_max, p.FIO2
                p.hpv_max, p.FIO2 = hmax * p.ae_hpv, p.ae_FIO2         # inflammation blunts HPV; supplemental O2
                acute = L.equilibrate(remodel=False)                        # no time to remodel
                p.hpv_max, p.FIO2 = hmax, fio2
                ra = record(L, acute, t, None); ra["ae"] = 1; ra["flood_frac"] = float(L.flood.mean())
                rows.append(ra)
                if ra["SaO2"] < p.ae_death_SaO2:
                    rows[-1]["death"] = "exacerbation"
                    break
                # survivors: part of the flooded units organise into scar, the rest recover
                scar = L.flood & (rng.random(L.N) < p.ae_fib)
                L.state[scar] = 1
                L.flood[:] = False
                L.sq = np.where(L.state == 0, sq0 * (1 - hit * (1 - p.ae_recover)), L.sq)
                L._rem_budget = bud
                st = L.equilibrate()
                rows.append(record(L, st, t, exercise(L) if do_exercise else None))
        if p.ae:
            r_ = rows[-1]
            cause = ("hypoxaemic" if r_["SaO2"] < p.death_SaO2 else
                     "ventilatory" if r_.get("power_ratio", 0) > p.W_lim else
                     "right_heart" if (p.rv and r_.get("Q", p.Q) < p.Q_death) else "")
            if cause:
                r_["death"] = cause
                break
        lost = 1 - rows[-1]["healthy"]
        for q in snapshots:
            if q not in snaps and lost >= q:
                snaps[q] = dict(state=L.state.copy(), Ql=st["vs"]["Ql"].copy(),
                                strain=st["air"]["strain"].copy(), PA=st["gas"]["PA"].copy())
        if verbose:
            r = rows[-1]
            print(f"t={t:.3f} healthy={r['healthy']:.3f} SaO2={r['SaO2']:.3f} zP={r['zP']:.2f} it={r['eq_it']}")
    return rows, snaps, L
