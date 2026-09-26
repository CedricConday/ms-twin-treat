"""Port of the Pernice 2020 RRMS immune-response net, as a deterministic ODE.

    Pernice S, Follia L, Maglione A, Pennisi M, Pappalardo F, Novelli F,
    Clerico M, Beccuti M, Cordero F, Rolla S. "Computational modeling of the
    immune response in multiple sclerosis using epimod framework."
    BMC Bioinformatics 2020;21(Suppl 17):550. PMID 33308135, PMC7734848. CC BY.

Sources, all vendored under docs/research/pernice2020/ and nothing else:

    model_description.txt      the paper's model, data and calibration sections
    additional_file_1.txt      S1.1 general transitions, Table S1, Table S2
    fig3_net.jpg               the net diagram (arcs)
    figS2_reproduction_target.png   deterministic 30-day HD vs MS solution

Scope, licence and the decision to port are in docs/PERNICE_PORT_SCOPE.md. No
code from the authors' repositories was read or used (no licence there).
NOTHING IS TUNED: every rate is Table S1 or S2 with the table beside it, the
initial marking is the paper's Table 2, the injection schedule is the paper's.
Where the diagram does not settle an arc, the two readings and the one taken
are marked `# UNRESOLVED`, and the reproduction of Figure S2 is what decided.

TIME UNIT IS THE HOUR
---------------------
Table S2 gives the six T/NK turnover rates as "1/24 h^-1" and DACDegradation
as 0.001444057, which is a 20-day half-life only if the unit is the hour
(Keizer 2010 is the citation for it). The cytokine clearance rates (0.03 to
0.09) return the Figure S2 spikes to baseline within about three days only at
per-hour scale. So the state advances in hours; the paper's "30 days" run is
720 h and the antigen injection "at the second day" is at 48 h.

THE NET, AS TRANSCRIBED
-----------------------
26 places: 21 named plus ODC unfolded into five myelination levels (L5 = Lmax,
the paper's initial 500, down to L1 = irreversibly damaged). 48 named
transitions in the diagram; the paper's count of 55 is reached if the coloured
instances of TeffKillsODC (4) and Remyelinization (3) are counted separately,
which is how Table S1 lists them. The 15 general transitions are the ones in
S1.1 and are transcribed from its functions; every other transition is mass
action on its input places at the Table S1/S2 rate.

Peripheral compartment (suffix _out), CNS (suffix _in), coupled through BBB:
Teff and Treg cross at pPass_BBB_teff = 0.005 and pPass_BBB_treg = 0.45 per
unit of BBB permeability, which is the 90-fold asymmetry the scope document
names as the reason to port.

READOUTS
--------
`irreversibly_damaged`  the ODC L1 count, the paper's own damage readout
`lesion_load`           the time-integrated count of ODC below Lmax; a
                        level-count proxy for lesion burden, continuous in
                        time, used for the treated/untreated ratio through
                        bricks/sormani in the same seat as the other models
Both are stated; neither is fitted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

# --------------------------------------------------------------------------- #
# Places. Order is the state-vector order.
# --------------------------------------------------------------------------- #
PLACES: tuple[str, ...] = (
    "Antigen", "Resting_Teff_temp", "Resting_Teff_out", "Teff_out", "EffectorMemory",
    "Resting_Treg_temp", "Resting_Treg_out", "Treg_out", "NK_out",
    "IL10_out", "IL17_out", "IFNg_out", "DAC", "BBB",
    "Resting_Teff_in", "Teff_in", "Resting_Treg_in", "Treg_in",
    "IL10_in", "IL17_in", "IFNg_in",
    "ODC_L1", "ODC_L2", "ODC_L3", "ODC_L4", "ODC_L5",
)
IDX = {name: i for i, name in enumerate(PLACES)}
N_ODC_LEVELS = 5
ODC_SLICE = slice(IDX["ODC_L1"], IDX["ODC_L5"] + 1)

# Table 2 of the paper: the only non-zero initial markings (cells per mm^3).
# Same marking for HD and MS. ODC: "all the 500 ODC with level Lmax".
INITIAL_MARKING: dict[str, float] = {
    "Resting_Teff_out": 1689.0,
    "Resting_Treg_out": 63.0,
    "NK_out": 30.0,
    "IL17_out": 8.0,
    "IL10_out": 13.0,
    "IFNg_out": 42.0,
    "IL17_in": 1.0,
    "IL10_in": 1.0,
    "IFNg_in": 1.0,
    "ODC_L5": 500.0,
}

# --------------------------------------------------------------------------- #
# Parameters. Table S1 (calibrated; MS column where it differs) and Table S2
# (fixed). Per hour. Names are the paper's.
# --------------------------------------------------------------------------- #
HOURLY_TURNOVER = 1.0 / 24.0   # Table S2 "1/24 h^-1"

HEALTHY_PARAMS: dict[str, float] = {
    # Table S1
    "pTeff_Activation": 0.015,      # TeffActivation_out/in, MemActivation   (MS: 0.018)
    "pTreg_Activation": 4e-04,      # TregActivation_in/out                 (MS: 7e-05)
    "pTreg_Dup": 0.006,             # TregDup_in/out
    "pTeff_Dup": 0.04,              # TeffDup_Asym_out, TeffDup_Sym_in/out
    "pTeff_KillsODC": 6e-04,        # TeffKillsODC, all four level instances
    "pTrkTe": 0.02,                 # TregKillsTeff_in/out
    "pTekA": 6e-04,                 # TeffKillsA
    "pPass_BBB_treg": 0.45,         # Treg_pass_BBB
    "pPass_BBB_teff": 0.005,        # Teff_pass_BBB
    "pNKkillsTeff": 0.01,           # NKkillsTeff_out
    "pNK_prod_IFNg": 0.03,          # NK_prod_IFNg
    "pNK_prod_IL10": 0.045,         # NK_prod_IL10
    "pIL17_BBB": 0.0115,            # IL17_BBB
    "pIL10_BBB": 0.0765,            # IL10_BBB
    "pRemyelinization": 0.01,       # Remyelinization, three level instances
    "pIL10Consuption": 0.09,        # IL10Consuption_out/in
    "pIL17Consuption": 0.03,        # IL17Consuption_out/in
    "pIFNgConsuption": 0.05,        # IFNgConsuption_out/in
    "Cifn": 20.0,                   # in the activation coefficient
    "CIL10": 10.0,                  # in the TregKillsTeff coefficient
    # Table S2
    "FromTimoREG": 0.317,
    "FromTimoEFF": 0.296,
    "NKdup": HOURLY_TURNOVER,
    "NKDegradation": HOURLY_TURNOVER,
    "Teff_death": HOURLY_TURNOVER,
    "Teff_to_NLT": HOURLY_TURNOVER,
    "Treg_death": HOURLY_TURNOVER,
    "Treg_to_NLT": HOURLY_TURNOVER,
    "Treg_prod_IL10": 0.05556,
    "Teff_prod_IL17": 0.00895,
    "Teff_prod_IFNg": 0.0466,
    "DACDegradation": 0.001444057,
    # S1.1 constants
    "rho_dup": 2.0 / 3.0,           # symmetric duplication probability
    "NK_setpoint": 30.0,            # NKentry keeps NK_out around 30
    "NK_entry_rate": 0.267,
    # DAC potency: the paper's two scenarios are 0.01 (weak) and 0.03 (strong)
    "pDACkill": 0.01,
}

# "the two sets of parameter values ... differ only in the values of
#  pTreg_Activation and pTeff_Activation" (model_description.txt, calibration)
MS_PARAMS: dict[str, float] = {**HEALTHY_PARAMS,
                               "pTeff_Activation": 0.018,
                               "pTreg_Activation": 7e-05}

HOURS_PER_DAY = 24.0
CALIBRATION_HOURS = 30.0 * HOURS_PER_DAY
CALIBRATION_INJECTIONS_DAYS: tuple[float, ...] = (2.0,)
TWO_YEAR_HOURS = 730.0 * HOURS_PER_DAY
TWO_YEAR_INJECTIONS_DAYS: tuple[float, ...] = (2.0, 67.0, 127.0, 295.0, 300.0, 303.0, 307.0, 600.0)
ANTIGEN_PER_INJECTION = 100.0

# --------------------------------------------------------------------------- #
# Arc readings the diagram does not settle. Each is a named choice with both
# readings stated; `DEFAULT_READINGS` is what Figure S2 reproduced with, and
# the module docstring's UNRESOLVED list is generated from this table.
# --------------------------------------------------------------------------- #
READINGS_DOC: dict[str, tuple[str, str]] = {
    # UNRESOLVED: does TeffActivation_out consume an Antigen token per
    # activation, or read it? S1.1 gives the rate with x_Antigen as a factor
    # and says nothing about the arc; the diagram's arc from Antigen is not
    # distinguishable from a read arc at its resolution.
    "activation_consumes_antigen": ("Antigen is consumed 1:1", "Antigen is read (catalyst)"),
    # UNRESOLVED: TregActivation_out/in is mass action (not in S1.1) on
    # Resting_Treg; the diagram draws a second arc into it from the activated
    # effector place. Reading A: Teff is a catalyst factor. Reading B: Resting
    # Treg alone (then 63 x 4e-4 per hour cannot produce Figure S2's Treg_out).
    "treg_activation_needs_teff": ("rate = p * Resting_Treg * Teff", "rate = p * Resting_Treg"),
    # UNRESOLVED: NKkillsTeff_out — is the NK cell consumed? Figure S2 shows
    # NK_out falling from 30 to under 10 during the effector peak and NKentry
    # is the only refill, so consumption is the reading that can move NK_out.
    "nk_consumed_in_kill": ("NK + Teff -> nothing", "NK + Teff -> NK"),
    # UNRESOLVED: IL17_BBB / IL10_BBB — is the cytokine consumed when it acts on
    # the barrier? Cytokines have their own consumption transitions, so the
    # barrier transitions are read as catalytic on the cytokine.
    "bbb_cytokine_consumed": ("cytokine consumed", "cytokine read"),
    # UNRESOLVED: which effector population crosses the barrier and where it
    # lands. The diagram draws Teff_pass_BBB from Teff_out into Resting_Teff_in
    # and the text says CNS Resting_Teff_in are (re)activated there.
    "crossing_lands_resting": ("Teff_out -> Resting_Teff_in", "Teff_out -> Teff_in"),
    # UNRESOLVED: Theta as extracted from the PDF has identical numerator and
    # denominator (a sign lost in extraction). The stated range [0.5, 1.5]
    # requires (pro - anti) / (pro + anti); that is what is used.
    "theta_plus_denominator": ("(IL17+IFNg-IL10)/(IL17+IFNg+IL10)", "as extracted, always 1"),
    # UNRESOLVED: FromTimoEff/Reg. Table 2 gives the _temp places no marking,
    # so a mass-action transition from them is dead. Reading A keeps it dead
    # (the paper's 30-day calibration has no thymic influx). Reading B would
    # make it a constant source, which Table 2 gives no marking for.
    "thymic_influx": ("dead: temp places start at 0", "constant source at the Table S2 rate"),
    # UNRESOLVED: TeffDup_Asym_out stoichiometry. Asymmetric division yields
    # one effector and one memory cell: Teff -> Teff + Memory.
    "asym_keeps_teff": ("Teff -> Teff + EffectorMemory", "Teff -> EffectorMemory"),
    # UNRESOLVED: MemActivation — is the memory cell consumed when it produces
    # an effector? The text: memory cells "remain in this compartment and are
    # able to respond faster to the infection reactivation". A consumed memory
    # pool answers one re-challenge and is gone; a read arc answers every one.
    # Figure S2 (one injection) cannot see this; the two-year figures (S4, 7)
    # show damage stepping up at each injection cluster, which is the read arc.
    "mem_activation_consumes_memory": ("EffectorMemory -> Teff_out", "EffectorMemory read; Teff_out produced"),
}

# What the Figure S2 scan decided (2026-09-26, all 256 combinations run against
# FIG_S2_LANDMARKS by the snippet recorded in BUILD_PLAN §8.4):
#   decided by the figure   activation_consumes_antigen=False (True loses 9 of
#                           12 landmarks), treg_activation_needs_teff=True,
#                           nk_consumed_in_kill=True, bbb_cytokine_consumed=False,
#                           thymic_influx=False, asym_keeps_teff=False (True
#                           puts the irreversibly damaged count at 40 against
#                           the figure's 13; False puts it in range).
#   NOT decided (identical landmarks either way): crossing_lands_resting and
#                           theta_plus_denominator. The diagram's arc and the
#                           stated range are kept for those two.
#   never reproduced        the Teff_out peak: about 1000 in both configurations
#                           against the figure's 400 (MS) and 100 (HD). With the
#                           transcribed activation function and the Table 2
#                           marking, the resting pool activates within the first
#                           hour after injection at every reading. Downstream
#                           landmarks (BBB, Teff_in, ODC) are in range regardless.
# The two-year figures (S4, 7) are NOT reproduced by either memory reading
# (measured 2026-09-26, BUILD_PLAN §8.4): with the memory cell consumed, the
# MS patient has one attack and no damage after day 30; with it read, damage
# accumulates but reaches all 500 ODC by day 310 in MS AND in HD, while the
# paper's two-year MS run reaches about 400 by day 600 and HD stays near 0.
# `MEMORY_READ_READINGS` is that second reading, scored as a labelled variant.
DEFAULT_READINGS: dict[str, bool] = {
    "activation_consumes_antigen": False,
    "treg_activation_needs_teff": True,
    "nk_consumed_in_kill": True,
    "bbb_cytokine_consumed": False,
    "crossing_lands_resting": True,
    "theta_plus_denominator": True,
    "thymic_influx": False,
    "asym_keeps_teff": False,
    "mem_activation_consumes_memory": True,
}
MEMORY_READ_READINGS: dict[str, bool] = {**DEFAULT_READINGS, "mem_activation_consumes_memory": False}
READING_VARIANTS: dict[str, dict[str, bool]] = {"s2": DEFAULT_READINGS, "memread": MEMORY_READ_READINGS}


# --------------------------------------------------------------------------- #
# Drug profiles: multipliers on the port's own named rates.
# --------------------------------------------------------------------------- #
DIALS: tuple[str, ...] = (
    "pTeff_Activation", "pTreg_Activation", "pTeff_Dup", "pTreg_Dup",
    "pTeff_KillsODC", "pTrkTe", "pPass_BBB_teff", "pPass_BBB_treg",
    "pIL17_BBB", "pIL10_BBB", "pRemyelinization",
    "Teff_death", "Treg_death", "Teff_prod_IL17", "Teff_prod_IFNg", "Treg_prod_IL10",
    "pNKkillsTeff", "pIFNgConsuption", "pIL17Consuption",
)


@dataclass(frozen=True)
class PerniceProfile:
    """Multipliers on named rates of the port; 1.0 is untouched. `dac_dose` is
    the paper's own drug: DAC units injected at treatment start."""

    label: str = "untreated"
    multipliers: dict = field(default_factory=dict)
    dac_dose: float = 0.0
    dac_potency: float = 0.01
    source: str = ""

    def __post_init__(self) -> None:
        unknown = set(self.multipliers) - set(DIALS)
        if unknown:
            raise ValueError(f"{self.label}: not dials of this port: {sorted(unknown)}")
        for k, v in self.multipliers.items():
            if not math.isfinite(v) or v < 0.0:
                raise ValueError(f"{self.label}: {k} must be finite and >= 0, got {v}")

    @property
    def is_untreated(self) -> bool:
        return all(v == 1.0 for v in self.multipliers.values()) and self.dac_dose == 0.0

    def as_dict(self) -> dict:
        return {"label": self.label, "source": self.source, "dac_dose": self.dac_dose,
                **self.multipliers}


UNTREATED_PROFILE = PerniceProfile()


def params_for(profile: PerniceProfile, base: dict[str, float] | None = None) -> dict[str, float]:
    p = dict(MS_PARAMS if base is None else base)
    for k, m in profile.multipliers.items():
        p[k] = p[k] * m
    p["pDACkill"] = profile.dac_potency
    return p


# --------------------------------------------------------------------------- #
# The right-hand side.
# --------------------------------------------------------------------------- #
def _theta(il17: float, ifng: float, il10: float, plus: bool) -> float:
    pro = il17 + ifng
    den = (pro + il10) if plus else (pro - il10)
    if den <= 0.0:
        return 1.0
    return 1.0 + 0.5 * (pro - il10) / den


def rhs(t: float, x: np.ndarray, p: dict[str, float], rd: dict[str, bool],
        t2inj: float) -> np.ndarray:
    (Ag, RTe_tmp, RTe_o, Te_o, Mem, RTr_tmp, RTr_o, Tr_o, NK,
     IL10_o, IL17_o, IFNg_o, DAC, BBB,
     RTe_i, Te_i, RTr_i, Tr_i, IL10_i, IL17_i, IFNg_i) = (max(v, 0.0) for v in x[:21])
    odc = np.maximum(x[ODC_SLICE], 0.0)          # L1..L5
    d = np.zeros_like(x)

    def add(place: str, v: float) -> None:
        d[IDX[place]] += v

    # ---- peripheral: thymic influx (Table S2; see READINGS_DOC["thymic_influx"])
    if rd["thymic_influx"]:
        add("Resting_Teff_out", p["FromTimoEFF"])
        add("Resting_Treg_out", p["FromTimoREG"])
    else:
        f = p["FromTimoEFF"] * RTe_tmp
        add("Resting_Teff_temp", -f)
        add("Resting_Teff_out", f)
        f = p["FromTimoREG"] * RTr_tmp
        add("Resting_Treg_temp", -f)
        add("Resting_Treg_out", f)

    # ---- S1.1 Activation (general)
    coef_o = 0.5 + math.exp(-IFNg_o / p["Cifn"])
    f_act = p["pTeff_Activation"] * RTe_o * Ag * coef_o
    add("Resting_Teff_out", -f_act)
    add("Teff_out", f_act)
    if rd["activation_consumes_antigen"]:
        add("Antigen", -f_act)
    if t >= t2inj:
        f_mem = 2.0 * p["pTeff_Activation"] * Ag * coef_o * Mem
        if rd["mem_activation_consumes_memory"]:
            add("EffectorMemory", -f_mem)
        add("Teff_out", f_mem)

    # ---- S1.1 Duplication (general)
    r_dup_o = Te_o * p["pTeff_Dup"]
    add("Teff_out", p["rho_dup"] * r_dup_o)                       # TeffDup_Sym_out: Te -> 2 Te
    f_asym = (1.0 - p["rho_dup"]) * r_dup_o                       # TeffDup_Asym_out
    add("EffectorMemory", f_asym)
    if not rd["asym_keeps_teff"]:
        add("Teff_out", -f_asym)

    # ---- Teff_out sinks and products (Table S2 mass action)
    add("Teff_out", -p["Teff_death"] * Te_o)
    add("IL17_out", p["Teff_prod_IL17"] * Te_o)
    add("IFNg_out", p["Teff_prod_IFNg"] * Te_o)
    add("IL17_out", -p["pIL17Consuption"] * IL17_o)
    add("IFNg_out", -p["pIFNgConsuption"] * IFNg_o)
    add("IL10_out", -p["pIL10Consuption"] * IL10_o)

    # ---- S1.1 Killing (general): TeffKillsA
    theta_o = _theta(IL17_o, IFNg_o, IL10_o, rd["theta_plus_denominator"])
    add("Antigen", -p["pTekA"] * theta_o * Ag * Te_o)

    # ---- Treg, peripheral
    f_tra = p["pTreg_Activation"] * RTr_o * (Te_o if rd["treg_activation_needs_teff"] else 1.0)
    add("Resting_Treg_out", -f_tra)
    add("Treg_out", f_tra)
    f_kill = p["pTrkTe"] * Tr_o * Te_o * (1.0 - math.exp(-IL10_o / p["CIL10"]))   # TregKillsTeff_out
    add("Teff_out", -f_kill)
    add("Treg_out", p["pTreg_Dup"] * Tr_o)                                          # TregDup_out
    add("Treg_out", -p["Treg_death"] * Tr_o)
    add("IL10_out", p["Treg_prod_IL10"] * Tr_o)

    # ---- NK (S1.1 NKentry general; the rest Table S1/S2)
    if NK <= p["NK_setpoint"]:
        add("NK_out", p["NK_entry_rate"] * (p["NK_setpoint"] - NK))
    add("NK_out", p["NKdup"] * NK - p["NKDegradation"] * NK)
    f_nk = p["pNKkillsTeff"] * NK * Te_o                                             # NKkillsTeff_out
    add("Teff_out", -f_nk)
    if rd["nk_consumed_in_kill"]:
        add("NK_out", -f_nk)
    add("IFNg_out", p["pNK_prod_IFNg"] * NK)
    add("IL10_out", p["pNK_prod_IL10"] * NK)

    # ---- DAC (S1.1 general): only acts when a dose is present
    if DAC > 0.0 and (Tr_o + Te_o) > 0.0:
        tot = Tr_o + Te_o
        add("Teff_out", -p["pDACkill"] * (Tr_o / tot) * DAC)     # DACkillTeff, as written in S1.1
        add("Treg_out", -p["pDACkill"] * (Te_o / tot) * DAC)     # DACkillTreg, as written in S1.1
    add("DAC", -p["DACDegradation"] * DAC)

    # ---- BBB
    f17 = p["pIL17_BBB"] * IL17_o
    add("BBB", f17)
    f10 = p["pIL10_BBB"] * IL10_o * BBB
    add("BBB", -f10)
    if rd["bbb_cytokine_consumed"]:
        add("IL17_out", -f17)
        add("IL10_out", -f10)
    f_te_pass = p["pPass_BBB_teff"] * Te_o * BBB
    f_tr_pass = p["pPass_BBB_treg"] * Tr_o * BBB
    add("Teff_out", -f_te_pass)
    add("Treg_out", -f_tr_pass)
    add("Resting_Teff_in" if rd["crossing_lands_resting"] else "Teff_in", f_te_pass)
    add("Resting_Treg_in", f_tr_pass)

    # ---- CNS: activation on ODC (all five levels, Table S1 lists le1..le5)
    coef_i = 0.5 + math.exp(-IFNg_i / p["Cifn"])
    f_act_i = p["pTeff_Activation"] * RTe_i * float(odc.sum()) * coef_i
    add("Resting_Teff_in", -f_act_i)
    add("Teff_in", f_act_i)
    add("Teff_in", p["rho_dup"] * Te_i * p["pTeff_Dup"])                            # TeffDup_Sym_in
    add("Teff_in", -p["Teff_to_NLT"] * Te_i)
    add("IL17_in", p["Teff_prod_IL17"] * Te_i - p["pIL17Consuption"] * IL17_i)
    add("IFNg_in", p["Teff_prod_IFNg"] * Te_i - p["pIFNgConsuption"] * IFNg_i)

    f_tra_i = p["pTreg_Activation"] * RTr_i * (Te_i if rd["treg_activation_needs_teff"] else 1.0)
    add("Resting_Treg_in", -f_tra_i)
    add("Treg_in", f_tra_i)
    add("Teff_in", -p["pTrkTe"] * Tr_i * Te_i * (1.0 - math.exp(-IL10_i / p["CIL10"])))
    add("Treg_in", p["pTreg_Dup"] * Tr_i - p["Treg_to_NLT"] * Tr_i)
    add("IL10_in", p["Treg_prod_IL10"] * Tr_i - p["pIL10Consuption"] * IL10_i)

    # ---- ODC: TeffKillsODC on levels 5..2 (four instances), Remyelinization 2..4 (three)
    theta_i = _theta(IL17_i, IFNg_i, IL10_i, rd["theta_plus_denominator"])
    dodc = np.zeros(N_ODC_LEVELS)
    for lvl in range(N_ODC_LEVELS - 1, 0, -1):                 # index 4 (L5) down to 1 (L2)
        f = p["pTeff_KillsODC"] * theta_i * odc[lvl] * Te_i
        dodc[lvl] -= f
        dodc[lvl - 1] += f
    for lvl in range(1, N_ODC_LEVELS - 1):                     # index 1 (L2) to 3 (L4)
        f = p["pRemyelinization"] * odc[lvl]
        dodc[lvl] -= f
        dodc[lvl + 1] += f
    d[ODC_SLICE] += dodc
    return d


# --------------------------------------------------------------------------- #
# Simulation with injection events.
# --------------------------------------------------------------------------- #
def initial_state() -> np.ndarray:
    x0 = np.zeros(len(PLACES))
    for k, v in INITIAL_MARKING.items():
        x0[IDX[k]] = v
    return x0


def simulate(profile: PerniceProfile = UNTREATED_PROFILE, *,
             params: dict[str, float] | None = None,
             t_end_hours: float = TWO_YEAR_HOURS,
             injections_days: tuple[float, ...] = TWO_YEAR_INJECTIONS_DAYS,
             treatment_start_day: float = 0.0,
             readings: dict[str, bool] | None = None,
             sample_hours: float = 6.0,
             rtol: float = 1e-6, atol: float = 1e-8) -> dict:
    """Deterministic solution. Returns hourly-sampled trajectories per place.

    Antigen injections add ANTIGEN_PER_INJECTION at each listed day. A DAC
    dose, if the profile carries one, is added at `treatment_start_day`. The
    profile's multipliers apply from t = 0 (the exam's arms are chronic).
    """
    rd = dict(DEFAULT_READINGS if readings is None else readings)
    p = params_for(profile, params)
    x = initial_state()
    events = sorted({d * HOURS_PER_DAY for d in injections_days} |
                    ({treatment_start_day * HOURS_PER_DAY} if profile.dac_dose > 0 else set()))
    inj_hours = sorted(d * HOURS_PER_DAY for d in injections_days)
    t2inj = inj_hours[1] if len(inj_hours) > 1 else math.inf
    grid = np.arange(0.0, t_end_hours + 1e-9, sample_hours)
    ts, xs = [], []
    t0 = 0.0
    ok = True
    for t_ev in [*[e for e in events if e < t_end_hours], t_end_hours]:
        if t_ev > t0:
            seg = grid[(grid >= t0) & (grid <= t_ev)]
            sol = solve_ivp(rhs, (t0, t_ev), x, method="LSODA", args=(p, rd, t2inj),
                            t_eval=seg if seg.size else None, rtol=rtol, atol=atol)
            if not sol.success:
                ok = False
                break
            if sol.t.size:
                ts.append(sol.t)
                xs.append(sol.y.T)
            x = sol.y[:, -1].copy()
            t0 = t_ev
        if t_ev in inj_hours:
            x[IDX["Antigen"]] += ANTIGEN_PER_INJECTION
        if profile.dac_dose > 0 and t_ev == treatment_start_day * HOURS_PER_DAY:
            x[IDX["DAC"]] += profile.dac_dose
    t = np.concatenate(ts) if ts else np.array([0.0])
    y = np.vstack(xs) if xs else x[None, :]
    traj = {name: y[:, i] for name, i in IDX.items()}
    traj["t_hours"] = t
    traj["t_days"] = t / HOURS_PER_DAY
    # LSODA leaves places that reach zero a few 1e-6 below it; that is solver
    # noise, not a regime exit. A place at -1e-3 cells or below is.
    traj["in_regime"] = bool(ok and np.all(np.isfinite(y)) and np.all(y > -1e-3))
    traj["params"] = p
    traj["readings"] = rd
    traj["profile"] = profile.as_dict()
    below_max = y[:, ODC_SLICE].sum(axis=1) - y[:, IDX["ODC_L5"]]
    traj["odc_below_max"] = below_max
    traj["lesion_load"] = float(np.trapezoid(below_max, t)) / HOURS_PER_DAY   # cell-days
    traj["irreversibly_damaged"] = float(y[-1, IDX["ODC_L1"]])
    return traj


def calibration_run(ms: bool, **kw) -> dict:
    """The paper's 30-day deterministic run behind Figure S2."""
    return simulate(UNTREATED_PROFILE, params=MS_PARAMS if ms else HEALTHY_PARAMS,
                    t_end_hours=CALIBRATION_HOURS, injections_days=CALIBRATION_INJECTIONS_DAYS,
                    sample_hours=1.0, **kw)


# --------------------------------------------------------------------------- #
# Figure S2 landmarks, read from the vendored figure. Used by the test and
# by the reading scan; values are approximate plot readings and are compared
# as ranges, never as fitted numbers.
# --------------------------------------------------------------------------- #
FIG_S2_LANDMARKS: dict[str, dict[str, tuple[float, float]]] = {
    # place: {"MS": (lo, hi), "HD": (lo, hi)} of the peak over the 30 days
    "Teff_out": {"MS": (150.0, 600.0), "HD": (40.0, 200.0)},
    "Treg_out": {"MS": (1.0, 12.0), "HD": (8.0, 40.0)},
    "BBB": {"MS": (0.6, 2.0), "HD": (0.0, 0.3)},
    "Teff_in": {"MS": (10.0, 60.0), "HD": (0.0, 6.0)},
    "ODC_L1": {"MS": (6.0, 26.0), "HD": (0.0, 1.0)},
    "IFNg_out": {"MS": (150.0, 700.0), "HD": (40.0, 250.0)},
}


def peaks(traj: dict) -> dict[str, float]:
    return {k: float(np.max(traj[k])) for k in FIG_S2_LANDMARKS}


def landmark_report(readings: dict[str, bool] | None = None) -> dict:
    out = {}
    for label, ms in (("MS", True), ("HD", False)):
        tr = calibration_run(ms, readings=readings)
        pk = peaks(tr)
        out[label] = {k: (pk[k], FIG_S2_LANDMARKS[k][label],
                          FIG_S2_LANDMARKS[k][label][0] <= pk[k] <= FIG_S2_LANDMARKS[k][label][1])
                      for k in FIG_S2_LANDMARKS}
        out[label]["in_regime"] = tr["in_regime"]
        after = tr["t_days"] > CALIBRATION_INJECTIONS_DAYS[0]
        cleared = after & (tr["Antigen"] < 1.0)
        out[label]["antigen_cleared_day"] = float(tr["t_days"][np.argmax(cleared)]) if cleared.any() else math.nan
        out[label]["effector_memory_final"] = float(tr["EffectorMemory"][-1])
        out[label]["nk_out_min"] = float(tr["NK_out"].min())
    return out


def landmarks_hit(readings: dict[str, bool] | None = None) -> int:
    rep = landmark_report(readings)
    return sum(1 for lab in ("MS", "HD") for k in FIG_S2_LANDMARKS if rep[lab][k][2])


if __name__ == "__main__":
    rep = landmark_report()
    for label in ("MS", "HD"):
        print(label, "in_regime", rep[label]["in_regime"],
              "antigen < 1 at day", round(rep[label]["antigen_cleared_day"], 1),
              "memory", round(rep[label]["effector_memory_final"], 1),
              "NK min", round(rep[label]["nk_out_min"], 1))
        for k, v in rep[label].items():
            if k in FIG_S2_LANDMARKS:
                peak, (lo, hi), ok = v
                print(f"   {k:<10} peak {peak:9.2f}   target [{lo:g}, {hi:g}]   {'ok' if ok else 'MISS'}")
