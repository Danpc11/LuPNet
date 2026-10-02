"""Every physical constant in one place. Pressures in mmHg (blood) and cmH2O (air),
blood flow in L/min, gas volumes in L. Vessel mechanics are computed in SI."""
from dataclasses import dataclass, asdict


@dataclass
class Params:
    G: int = 10                    # tree depth: N = 2^G terminal lung units
    geometry: str = "symmetric"    # "symmetric" (heap tree) | "anatomical" (3D volume-filling tree)
    anat_N: int = 1024             # anatomical: number of terminal units
    anat_seed: int = 0
    anat_a: float = 7.0            # anatomical: lateral half-width, cm
    anat_c: float = 9.0            # anatomical: antero-posterior half-depth, cm
    anat_H: float = 24.0           # anatomical: height, cm
    # ---- circulation (healthy, rest)
    Q: float = 5.0                 # cardiac output, L/min
    LAP: float = 8.0               # left atrial pressure, mmHg
    dP0: float = 6.0               # healthy mPAP - LAP, mmHg
    split: tuple = (0.40, 0.35, 0.25)   # arteriole / capillary / venule share of dP0
    mu_b: float = 3.5e-3           # blood viscosity, Pa s
    r_pa: float = 0.014            # main pulmonary artery radius, m
    L_pa: float = 0.05             # main pulmonary artery length, m
    L_term: float = 0.005          # equivalent length of the lumped terminal arteriole, m
    alpha: float = 0.02            # vessel distensibility, 1/mmHg
    alpha_c: float = 0.015         # capillary sheet distensibility, 1/mmHg
    Vcap: float = 0.075            # capillary blood volume, L
    # ---- gravity in perfusion, West zones, extravascular pressures (off = earlier campaigns)
    perf_grav: bool = False
    lung_height_cm: float = 30.0   # upright base-to-apex height
    heart_frac: float = 0.5        # height of the pulmonary artery / left atrium as a fraction of lung height
    rho_g: float = 0.77            # blood hydrostatic gradient, mmHg/cm
    P_alv: float = 0.0             # mean alveolar pressure, cmH2O
    P_pl: float = -5.0             # pleural pressure, cmH2O (FRC, spontaneous breathing)
    kappa_ev: float = 0.3          # perivascular pressure below pleural by kappa_ev x transpulmonary pressure
    sheet_a: float = 0.03          # capillary sheet compliance, 1/mmHg (H = 1 + a (P - P_alv))
    sheet_w: float = 0.5           # smoothing width of sheet collapse, mmHg
    sigma_sep: float = 0.0         # septal stretch: alveolar-vessel extravascular pressure + sigma (P_L - PL_FRC)
    PL_FRC: float = 5.0            # transpulmonary pressure at FRC, cmH2O
    # ---- shear set-point rule (vascular-shear-setpoint)
    b: float = 0.75                # metabolic cost exponent
    kappa: float = 0.5
    clip: float = 0.05
    remodel: bool = True
    remodel_steps: int = 0         # >0: at most this many remodelling steps per progression step (finite rate)
    veins: str = "mirror"          # "mirror" (venous tree mirrors arteries) | "septal" | "mirror_general" (check)
    vein_occl: float = 0.0         # septal venule conductance loss per stiff neighbour fraction (septal only)
    remodel_conduits: bool = True  # False: conduits keep their structure (e.g. encased in scar)
    # ---- traction of extra-alveolar vessels by stiff parenchyma
    traction: float = 0.0          # 0 off, 1 full: perivascular pressure follows the rise in recoil
    P_el_FRC: float = 5.0          # healthy elastic recoil at FRC, cmH2O
    # ---- maladaptive (non-set-point) arterial remodelling: rest point z = 1 + M
    eta_p: float = 0.0             # pressure (stretch) drive
    eta_s: float = 0.0             # scar micro-environment drive
    mal_depth: int = 7             # conduits at this depth and deeper also respond
    # ---- hypoxic pulmonary vasoconstriction on the terminal arteriole
    hpv: bool = True
    hpv_max: float = 3.0
    hpv_P50: float = 55.0          # mmHg
    hpv_n: float = 5.0
    # ---- airways and parenchyma
    mu_a: float = 1.8e-5           # air viscosity, Pa s
    r_tr: float = 0.009            # trachea radius, m
    L_tr: float = 0.12             # trachea length, m
    aw_r_exp: float = 0.38         # airway radius ~ 2^(-aw_r_exp g) (Weibel-like)
    aw_L_exp: float = 0.80         # airway length ~ 2^(-aw_L_exp g)
    k_inertial: float = 4.0        # lumped Pedley-type factor on central airway resistance
    R_periph: float = 0.4          # small-airway resistance, cmH2O s/L, whole lung
    R_tissue: float = 0.5          # tissue resistance, cmH2O s/L, whole lung
    C_lung: float = 0.20           # lung compliance, L/cmH2O
    C_cw: float = 0.20             # chest-wall compliance, L/cmH2O
    R_cw: float = 0.3
    FRC: float = 2.5               # L
    VD: float = 0.15               # anatomical dead space, L
    grav: float = 0.4              # base-to-apex spread in unit compliance (gravity)
    grav_V: float = 0.3            # apex-to-base spread in unit FRC volume
    VA: float = 4.2                # alveolar ventilation demand, L/min (held constant)
    f_grid: tuple = (6.0, 60.0, 109)    # breaths/min searched for minimum power
    # ---- gas exchange
    PB: float = 760.0              # barometric pressure, mmHg
    FIO2: float = 0.2093
    VO2: float = 250.0             # mL/min
    Hb: float = 15.0               # g/dL
    P50: float = 26.8
    hill_n: float = 2.7
    t_eq: float = 0.25             # s: healthy membrane equilibrates in ~t_eq (3 e-folds)
    # ---- exercise probe
    ex_Q: float = 2.5
    ex_VO2: float = 3.0
    ex_VA: float = 3.0
    # ---- unit states: healthy 0, fibrotic 1, honeycomb 2
    fib_C: float = 0.15            # compliance factor
    fib_V: float = 0.40            # FRC volume factor
    fib_c: float = 0.25            # surviving capillary fraction
    fib_D: float = 0.20            # membrane diffusing factor
    hc_C: float = 0.60
    hc_V: float = 1.00
    hc_c: float = 0.0
    hc_D: float = 0.0
    nbr_D: float = 1.0             # <1: healthy units lose diffusing capacity next to scar (variant)
    # ---- right ventricle (end-systolic elastance; off by default = fixed cardiac output)
    rv: bool = False
    HR: float = 70.0
    # shared target law (InFlow): one tau0 for every vessel, exact Womersley sensing on the first harmonic
    vasc_A: float = 0.0           # bounded vasculopathy: maximum lesion resistance in units of the terminal-arteriole
                                  # resistance; R_les = R_ta0 vasc_A (1 - exp(-(t_yr + t_vasc0) / vasc_tau))
    vasc_tau: float = 5.0         # time constant of the bounded vasculopathy (years)
    cap_frac: float = 0.0         # capillary loss coupled to the vasculopathy: capillary factor 1/(1 + cap_frac (exp(k t) - 1))
    t_vasc0: float = 0.0          # years of vasculopathy already elapsed at the start (independent of fibrosis)
    k_vasc: float = 0.0           # structural small-artery vasculopathy, independent of fibrosis (1/yr):
                                  # a non-remodelable series lesion R_les = R_ta0 (exp(k_vasc t) - 1) in every unit
    law: str = "local"            # "local" (per-vessel tau0 from the healthy tree) | "shared"
    w1: float = 2.865             # weight of the first-harmonic wall-shear amplitude (Feaver 2013 fit)
    phi1_a: float = 0.8           # first-harmonic flow amplitude / mean flow, pulmonary arteries (effective value)
    phi1_v: float = 0.4           # same, pulmonary veins (assumed smaller)
    b_ref_calib: float = 0.75     # exponent of the reference tree used to fix capillary and venule resistances
    Ees0: float = 0.40             # mmHg/mL, healthy RV end-systolic elastance
    Ees_max: float = 0.70          # ceiling of hypertrophic (homeometric) adaptation
    rv_ratio: float = 2.0          # Ees/Ea the RV tries to keep by hypertrophy
    V0_rv: float = 24.0            # mL
    EDV0: float = 130.0            # mL, healthy end-diastolic volume
    EDV_min: float = 90.0
    EDV_max: float = 180.0         # dilation ceiling
    RAP0: float = 4.0              # mmHg at EDV0
    k_edpvr: float = 0.011         # 1/mL
    HR_ex: float = 150.0
    ex_inotropy: float = 2.0       # healthy inotropic reserve; shrinks to 1 as Ees -> Ees_max
    ex_EDV: float = 20.0
    stop_Q: float = 2.5            # stop a run when resting output falls below this (L/min)
    # ---- parenchyma: surfactant, recruitment/derecruitment, collapse induration (off = earlier campaigns)
    parenchyma: str = "linear"     # "linear" | "surfactant"
    pl_grad: float = 0.25          # pleural pressure gradient, cmH2O per cm of height
    PEEP: float = 0.0              # extra end-expiratory transpulmonary pressure, cmH2O (intervention)
    Pc_off: float = -2.0           # closing pressure of a normal unit, cmH2O (never reached at FRC)
    Pc_slope: float = 12.0         # rise of closing pressure as surfactant is lost
    Po_off: float = 2.0            # opening minus closing pressure, normal unit
    Po_slope: float = 10.0         # extra hysteresis as surfactant is lost
    col_C: float = 0.02            # compliance factor of a closed unit
    col_V: float = 0.2             # volume factor of a closed unit
    col_c: float = 0.5             # capillary factor of a closed unit (partial compression)
    cyc_vent: float = 0.5          # ventilating fraction of a cyclic unit
    cyc_strain: float = 2.0        # strain multiplier of reopening every breath
    surf_Vcap: float = 0.7         # fraction of capillary bed lost at zero surfactant (dog-lung data: ~70%)
    dst_mid: float = 3.5           # pericapillary pressure below alveolar at mid volume, mmHg
    gamma_gain: float = 2.0        # surface tension x (1 + gain (1 - quality))
    lam_s: float = 1.0             # AT2 surfactant decline rate (time unit 1/lam0)
    beta_s: float = 3.0            # strain dependence of AT2 decline
    kn_s: float = 1.0              # extra AT2 decline per stiff-or-collapsed neighbour fraction
    at2_sigma: float = 0.7         # lognormal heterogeneity of AT2 vulnerability
    lam_ci: float = 2.0            # closed -> fibrotic (collapse induration)
    lam_ci_cyc: float = 0.5        # cyclic -> fibrotic (atelectrauma)
    lam0_direct: float = 0.3       # direct strain-driven fibrosis in surfactant mode
    # ---- two time scales and aging
    breath_model: str = "linear"   # "linear" (frequency domain) | "cycle" (time-domain breath, Bates-Irvin)
    s_open: float = 0.5            # reopening speed, 1/(cmH2O s)
    s_close: float = 0.3           # closing speed, 1/(cmH2O s)
    lam_ae: float = 0.5            # atelectrauma hazard per reopening event per breath (per model time)
    years_per_unit: float = 1.0    # model time unit in years (calibrated to FVC decline)
    tau_rem_weeks: float = 0.0     # >0: set-point remodelling time constant in weeks (finite, physical rate)
    age_M: float = 0.0             # aging: inward set-point offset (rest z = 1 + age_M)
    age_sigma: float = 0.0         # aging: per-vessel lognormal set-point error
    FVC0_mL: float = 3800.0        # healthy FVC for mL conversion (older adult)
    # ---- gas exchange physics and acute exacerbations
    rf: bool = False               # Roughton-Forster unit diffusing capacity (membrane + capillary volume)
    mem_gain: float = 0.0          # membrane conductance loss per unit of epithelial (surfactant) damage
    mem_field: float = 0.0         # diffuse membrane loss x sqrt(damaged fraction) in not-yet-scarred units
    vc_field: float = 0.0          # diffuse capillary loss x sqrt(damaged fraction) in not-yet-scarred units
    ae: bool = False               # acute exacerbations as Poisson shocks
    ae_rate0: float = 0.05         # per year at normal FVC
    ae_rate_slope: float = 0.20    # extra per year per unit of FVC lost (lower FVC = higher risk)
    ae_hit: float = 0.6            # acute fractional surfactant loss in open units (diffuse alveolar damage)
    ae_recover: float = 0.5        # fraction of the acute loss recovered by survivors
    death_SaO2: float = 0.80       # resting SaO2 on room air below this = death from chronic respiratory failure
    ae_flood_lo: float = 0.1       # fraction of open units flooded in an exacerbation, uniform [lo, hi]
    ae_flood_hi: float = 0.5
    ae_hpv: float = 0.3            # HPV strength during an exacerbation (inflammation blunts it)
    ae_FIO2: float = 0.6           # supplemental oxygen during an exacerbation
    ae_death_SaO2: float = 0.85    # SaO2 on supplemental O2 below this = death in the exacerbation (refractory)
    ae_fib: float = 0.5            # fraction of flooded units that organise into scar in survivors
    W_lim: float = 5.0             # death when the breathing power needed for PaCO2 40 exceeds W_lim x rest
    Q_death: float = 3.5           # death when resting cardiac output falls below this (with rv)
    # ---- constant injury flux and pulsatile shear sensing
    flux_gamma: float = 0.0        # per-unit injury rates x healthy_fraction^(-gamma) (1 = constant absolute flux)
    pulsatile: bool = False        # set-point senses pulsatile (Womersley) wall shear
    pulse_frac: float = 0.8        # pulsatile fraction of arterial flow
    pulse_frac_v: float = 0.3      # pulsatile fraction on the venous side
    rho_b: float = 1060.0          # blood density, kg/m^3
    wom_c: float = 2.83            # Womersley transition (2 sqrt 2)
    pulse_damp: float = 0.0        # >0: pulsatile fraction decays as exp(-generation / pulse_damp)
    # ---- epithelial mechanobiology: senescence, aberrant basaloid cells, fibroblast mechanosensing
    epi_bio: bool = False
    k_sen0: float = 0.3            # background AT2 senescence (age / telomere), per model time
    k_sen_strain: float = 1.0      # stretch-induced senescence per unit excess strain
    k_sasp: float = 0.5            # paracrine senescence from senescent neighbours (SASP)
    g_sen: float = 2.0             # surfactant decline x (1 + g_sen x senescence)
    k_ab: float = 1.0              # conversion to aberrant basaloid: k_ab x senescence x injury x stretch
    k_fib_ab: float = 1.0          # fibrosis hazard per aberrant fraction of the niche
    k_yap: float = 2.0             # fibroblast mechanosensing gain (YAP/TAZ, PIEZO)
    yap_K: float = 0.3             # stiff-neighbour fraction at half activation
    pleural_stress: float = 0.0    # extra strain on subpleural units (stress concentration at the pleura)
    senolytic_year: float = -1.0   # >= 0: clear a fraction of senescent cells once at this year
    senolytic_frac: float = 0.5
    antifib: float = 1.0           # multiplier on fibroblast-driven conversion hazards (antifibrotic)
    # ---- progression
    lam0: float = 1.0              # background hazard per healthy unit (sets the time unit)
    beta: float = 3.0              # strain feedback: hazard ~ exp(beta (s/s_ref - 1))
    kn: float = 0.5                # extra hazard per stiff grid neighbour fraction
    lam_hc: float = 0.5            # fibrotic -> honeycomb rate
    max_frac_step: float = 0.02    # at most ~2% of units change per step
    stop_healthy: float = 0.25
    stop_mPAP: float = 50.0        # stop a run in decompensated pulmonary hypertension
    t_max: float = 3.0
    seed: int = 0
    eq_maxit: int = 4000           # cap on equilibrium iterations (non-convergence is recorded as zerr)

    def dict(self):
        return asdict(self)
