"""Arms as multipliers on the Pernice port's own named rates, from pharmacology.

The rule is bricks/profiles.py's: WHICH rate and WHICH direction come from the
drug's independent pharmacology, never from its trial; HOW MUCH is the one
shared stub the exam sweeps. Every source below is the one already cited in
bricks/profiles.py for that arm; this file only says where it lands on a net
that has two compartments, a barrier, two T-cell types with their own death
and passage rates, NK cells, and a remyelination rate.

WHAT THE PORT SEPARATES THAT VELEZ COULD NOT (docs/PERNICE_PORT_SCOPE.md)
------------------------------------------------------------------------
    depletion          Teff_death and Treg_death up      alemtuzumab, cladribine
    effector deletion  Teff_death up only                Tovaxin
    transit block      pPass_BBB_teff and _treg down     natalizumab
    activation block   pTeff_Activation down             abatacept, anti-CD20
    proliferation      pTeff_Dup down                    IFN-beta, teriflunomide, ...
    the paper's drug   DAC dose through DACkillTeff/Treg daclizumab

WHAT IT STILL LUMPS
-------------------
Lymph node and blood vessel are one compartment, so an egress block has no
step to act on. The S1P modulators are lumped onto the barrier-passage dials
with natalizumab: sequestered cells cannot reach the barrier, which is the
only representable consequence. `LUMPED` names them.

Atacicept has no rate: the net has no B lineage. bricks/profiles.py lumped it
onto effector loss; with cell-type resolution that lump is no longer
defensible, so the arm is scored as "no handle" (see `UNASSIGNABLE`).

DIRECTION KEY IN THE PATTERN
----------------------------
As in bricks/profiles_minimal.py the pattern carries the sign, because both
signs of pTeff_Activation are in the arm set (abatacept down, APL up).
"""

from __future__ import annotations

from bricks.profiles import ENHANCE, STUB_MAGNITUDE, SUPPRESS
from bricks.qsp_pernice import DIALS, PerniceProfile

__all__ = ["PROFILES", "LUMPED", "UNASSIGNABLE", "STUB_MAGNITUDE", "DAC_MAX_DOSE",
           "touched_points", "profile_at"]

# The paper's largest daclizumab scenario (model_description.txt: doses 1000,
# 2000, 5000, 10000, 15000). Potency s maps to dose s * DAC_MAX_DOSE, weak
# potency 0.01, the paper's own function for both kills.
DAC_MAX_DOSE = 15000.0
DAC_WEAK_POTENCY = 0.01


def _p(name: str, source: str, dac_dose: float = 0.0, **mult: float) -> PerniceProfile:
    return PerniceProfile(label=name, source=source, multipliers=dict(mult),
                          dac_dose=dac_dose, dac_potency=DAC_WEAK_POTENCY)


LUMPED: dict[str, str] = {
    "fingolimod": "lymph-node egress block -> barrier passage, both cell types; no lymph-node compartment",
    "ponesimod": "as fingolimod",
    "ozanimod": "as fingolimod",
    "ofatumumab": "anti-CD20 -> activation, as ocrelizumab; no B cells",
    "ublituximab": "as ocrelizumab",
    "rituximab": "as ocrelizumab",
    "ocrelizumab": "anti-CD20 -> activation (Martinez-Pasamar 2013's reading); no B cells",
}

UNASSIGNABLE: dict[str, str] = {
    "atacicept": (
        "TACI-Ig blocks BAFF/APRIL and depletes plasma cells. The net has no B cell, "
        "plasma cell or antibody place. bricks/profiles.py lumped this onto effector "
        "loss; this port resolves cell types, and asserting that a B-lineage agent "
        "kills T effectors would be a mechanism claim with no source. No rate."
    ),
}

UNTREATED = _p("untreated", "control arm")

# ---- proliferation block (Velez alpha_E) -> pTeff_Dup
IFN_BETA = _p("IFN-beta", "reduces T-cell activation and proliferation; the proliferation "
              "rate is the one bricks/profiles.py chose (alpha_E), here pTeff_Dup", pTeff_Dup=SUPPRESS)
IFN_BETA_1B = _p("IFN-beta-1b", "as IFN-beta", pTeff_Dup=SUPPRESS)
TERIFLUNOMIDE = _p("teriflunomide", "DHODH inhibitor, blocks proliferation of activated "
                   "lymphocytes; pTeff_Dup. The Treg axis stays off (contested, see profiles.py)",
                   pTeff_Dup=SUPPRESS)
DIMETHYL_FUMARATE = _p("dimethyl fumarate", "reduces effector lymphocyte expansion; pTeff_Dup",
                       pTeff_Dup=SUPPRESS)
USTEKINUMAB = _p("ustekinumab", "anti-IL-12/23, fewer Th1/Th17 effectors; the port's Teff pool "
                 "is the IFNg/IL-17 producers and their expansion is pTeff_Dup", pTeff_Dup=SUPPRESS)

# ---- activation step
GLATIRAMER = _p("glatiramer acetate", "displaces myelin antigen from MHC II (activation down) "
                "and restores the regulatory population (Treg activation up); the two axes "
                "bricks/profiles.py gave it, on the port's own rates",
                pTeff_Activation=SUPPRESS, pTreg_Activation=ENHANCE)
ABATACEPT = _p("abatacept", "CTLA4-Ig blocks CD28 costimulation, the priming of resting T cells: "
               "pTeff_Activation, the same step bricks/profiles.py chose (delta). No Treg axis.",
               pTeff_Activation=SUPPRESS)
OCRELIZUMAB = _p("ocrelizumab", "anti-CD20; Martinez-Pasamar 2013 (PMC3651362) read B-cell "
                 "depletion as preventing uncontrolled T-effector activation without "
                 "strengthening Treg: pTeff_Activation. The magnitude 0.85 they fitted was a "
                 "threshold shift in a model this port does not share, so the stub applies.",
                 pTeff_Activation=SUPPRESS)
OFATUMUMAB = _p("ofatumumab", "as ocrelizumab, same target", pTeff_Activation=SUPPRESS)
UBLITUXIMAB = _p("ublituximab", "as ocrelizumab, same target", pTeff_Activation=SUPPRESS)
RITUXIMAB = _p("rituximab", "as ocrelizumab, same target", pTeff_Activation=SUPPRESS)
APL_CGP77116 = _p("APL CGP77116", "altered peptide ligand, encephalitogenic at the TCR "
                  "(Bielekova 2000): activation of resting effectors up", pTeff_Activation=ENHANCE)

# ---- barrier passage
NATALIZUMAB = _p("natalizumab", "anti-alpha4-integrin blocks VLA-4-dependent transit of "
                 "lymphocytes across the barrier. VLA-4 is on effectors and regulators alike, "
                 "so both passage rates fall; the asymmetry is the port's (0.005 vs 0.45).",
                 pPass_BBB_teff=SUPPRESS, pPass_BBB_treg=SUPPRESS)
FINGOLIMOD = _p("fingolimod", "S1P egress block, LUMPED onto barrier passage (see LUMPED)",
                pPass_BBB_teff=SUPPRESS, pPass_BBB_treg=SUPPRESS)
PONESIMOD = _p("ponesimod", "as fingolimod", pPass_BBB_teff=SUPPRESS, pPass_BBB_treg=SUPPRESS)
OZANIMOD = _p("ozanimod", "as fingolimod", pPass_BBB_teff=SUPPRESS, pPass_BBB_treg=SUPPRESS)

# ---- depletion
ALEMTUZUMAB = _p("alemtuzumab", "anti-CD52 depletes T cells of both types by CDC/ADCC: "
                 "Teff_death and Treg_death up. The Treg-biased repopulation (Cossburn) is a "
                 "timescale effect not represented; both death rates get the same stub.",
                 Teff_death=ENHANCE, Treg_death=ENHANCE)
CLADRIBINE = _p("cladribine", "lymphocyte-depleting purine analogue, both T-cell types",
                Teff_death=ENHANCE, Treg_death=ENHANCE)
TOVAXIN = _p("Tovaxin", "T-cell vaccination deletes myelin-reactive EFFECTOR clones: "
             "Teff_death up, regulators untouched", Teff_death=ENHANCE)

# ---- cytokine level
IFN_GAMMA = _p("IFN-gamma", "exogenous type II interferon raises IFN-gamma levels; the port "
               "carries IFNg as a state with two roles (slows activation, speeds killing), so "
               "the literal handle is its clearance rate down, and the sign of the outcome is "
               "the dynamics' to decide", pIFNgConsuption=SUPPRESS)

# ---- two-plus dials
LENERCEPT = _p("lenercept", "TNF neutralisation: effector expansion down (pTeff_Dup), regulatory "
               "expansion down (pTreg_Dup; TNF-deficient mice get severe EAE, Liu 1998), and "
               "TNFR2-dependent remyelination down (Arnett 2001), which this port has a rate for",
               pTeff_Dup=SUPPRESS, pTreg_Dup=SUPPRESS, pRemyelinization=SUPPRESS)

# ---- the paper's own drug
DACLIZUMAB = _p("daclizumab", "anti-CD25. The port models DAC itself: a dose that kills activated "
                "effectors and regulators through the paper's DACkillTeff/DACkillTreg functions, "
                "the two-sided effect bricks/profiles.py could only flag. The CD56bright NK "
                "expansion (Bielekova 2006) has a place here too: NK killing of effectors up.",
                dac_dose=DAC_MAX_DOSE * STUB_MAGNITUDE, pNKkillsTeff=ENHANCE)

# ---- no rate
ATACICEPT = _p("atacicept", UNASSIGNABLE["atacicept"])

PROFILES: dict[str, PerniceProfile] = {p.label: p for p in (
    UNTREATED, IFN_BETA, IFN_BETA_1B, TERIFLUNOMIDE, DIMETHYL_FUMARATE, USTEKINUMAB,
    GLATIRAMER, ABATACEPT, OCRELIZUMAB, OFATUMUMAB, UBLITUXIMAB, RITUXIMAB, APL_CGP77116,
    NATALIZUMAB, FINGOLIMOD, PONESIMOD, OZANIMOD,
    ALEMTUZUMAB, CLADRIBINE, TOVAXIN, IFN_GAMMA, LENERCEPT, DACLIZUMAB, ATACICEPT,
)}


def touched_points(profile: PerniceProfile) -> tuple[str, ...]:
    """Signed pattern over the port's dials, plus 'dac' when a dose is given."""
    out = []
    for k in DIALS:
        v = profile.multipliers.get(k, 1.0)
        if v < 1.0:
            out.append(f"{k}-")
        elif v > 1.0:
            out.append(f"{k}+")
    if profile.dac_dose > 0:
        out.append("dac")
    return tuple(out)


def profile_at(arm: str, s: float) -> PerniceProfile:
    """The arm at potency s: (1 - s) on suppressed rates, (1 + s) on enhanced,
    dose s * DAC_MAX_DOSE for the paper's drug. Direction never changes."""
    base = PROFILES[arm]
    mult = {k: ((1.0 - s) if v < 1.0 else (1.0 + s)) for k, v in base.multipliers.items()}
    dose = DAC_MAX_DOSE * s if base.dac_dose > 0 else 0.0
    return PerniceProfile(label=f"{arm}@{s:g}", source="exam v2 potency sweep",
                          multipliers=mult, dac_dose=dose, dac_potency=base.dac_potency)


if __name__ == "__main__":
    width = max(len(n) for n in PROFILES)
    for name, prof in PROFILES.items():
        tag = "  [LUMPED]" if name in LUMPED else ("  [NO RATE]" if name in UNASSIGNABLE else "")
        print(f"  {name:<{width}}  {touched_points(prof) or '-'}{tag}")
    print(f"\n  {len(PROFILES) - 1} arms, "
          f"{len({touched_points(p) for p in PROFILES.values()}) - 1} treated patterns")
