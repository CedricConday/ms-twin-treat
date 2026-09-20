"""B6 — brain PBPK, GROUNDED in the Verscheijden 2019 four-compartment model.

`bricks/barrier.py` is a hand-rolled 2-3 compartment ODE whose baseline was
calibrated to one published number (Pardridge 2019, ~0.15% CNS:serum for an
antibody) with the remaining rate constants left illustrative. BUILD_PLAN §6
records the decision that produced it: *"reimplement a minimal 2-3 compartment
blood<->CSF<->CNS ODE in scipy (don't install R for Verscheijden now)"*. That
deferral was never revisited, and `docs/ARCHITECTURE.md` has carried
"Verscheijden ships runnable R in supplement" ever since.

This is that model, transcribed from the authors' own R source:

    Verscheijden LFM, Koenderink JB, Johnson TN, de Wildt SN, Russel FGM.
    "Development of a physiologically-based pharmacokinetic pediatric brain
     model for prediction of cerebrospinal fluid drug concentrations and the
     influence of meningitis."
    PLoS Computational Biology 2019;15(6):e1007117.
    doi:10.1371/journal.pcbi.1007117   PMC6592555

    Source: S1 File, `pcbi.1007117.s001.R` (793 lines, deSolve), fetched from
    the Europe PMC supplementary package for PMC6592555.

No R was installed. The brain sub-model was transcribed equation by equation.

THE FOUR COMPARTMENTS
---------------------
  Cbb    brain blood
  Cbm    brain mass
  Cccsf  cranial CSF
  Cscsf  spinal CSF

  dCbb/dt   = (Qbrain*(Cart - Cbb) + PSb*(Fubm*Cbm - Fubb*Cbb)
               + PSc*(Fuccsf*Cccsf - Fubb*Cbb)
               + Q_ssink*Cscsf + Q_csink*Cccsf) / Vbb
  dCbm/dt   = (PSb*(Fubb*Cbb - Fubm*Cbm) + PSe*(Fuccsf*Cccsf - Fubm*Cbm)
               - Q_bulk*Cbm) / Vbm
  dCccsf/dt = (PSe*(Fubm*Cbm - Fuccsf*Cccsf) + PSc*(Fubb*Cbb - Fuccsf*Cccsf)
               + Q_bulk*Cbm + Q_sout*Cscsf - Qsin*Cccsf - Q_csink*Cccsf) / Vccsf
  dCscsf/dt = (Qsin*Cccsf - Q_sout*Cscsf - Q_ssink*Cscsf) / Vscsf

Three permeability-surface products carry the barriers, and they are the reason
this is worth having over a two-compartment sketch: the blood-brain barrier
(`PSb`) and the blood-CSF barrier (`PSc`) are separate, and CSF turnover
(`Qsin`, `Q_sout`, sinks) is explicit rather than folded into one rate.

WHAT THIS MODEL IS PARAMETERISED FOR, AND IT IS NOT MS
-------------------------------------------------------
Stated plainly because the structure being published does not make the numbers
transferable:

  * **Pediatric.** The cohort is ages 0.25-12.75 y, with height/weight from
    Simcyp pediatric equations, and brain blood flow is an explicit function of
    age: `Qbrain = Qcardiac * (10 + 2290*(exp(-0.608*age) - exp(-0.639*age)))/100`.
    Evaluated at an adult age that expression is outside anything it was fitted
    on. `DEFAULT_ADULT` below therefore uses published ADULT volumes and flows
    and does NOT extrapolate the pediatric curve — that substitution is this
    repo's, not the paper's, and is flagged in the constant itself.
  * **A passively-diffusing small molecule.** The PS values are for
    paracetamol, dosed 15 mg/kg IV over 10 minutes. **Therapeutic antibodies do
    not passively diffuse**, and every drug in this repo's arm set that would
    need CNS entry is an antibody. Using `PSb` unchanged for one would overstate
    brain entry by orders of magnitude, which is exactly the error
    `bricks/barrier.py` was calibrated against Pardridge to avoid.

So the honest split: **the structure and the CSF turnover are ported; `PSb` for
an antibody must come from Pardridge, not from here.** `antibody_psb()` does
that conversion and refuses to guess.

CURRENT STATUS IN THE PIPELINE
------------------------------
`BarrierStage` still owns B6. This module is a grounded replacement for its
internals, not yet wired, and there is a reason to be honest about rather than
paper over: **no arm in `bricks/intervention.py` sets `cns_required`**, so the
readout's delivery gating `(1 - effective) * treat` is currently a no-op for
every scored arm. A better barrier model changes no gate number today. It
replaces invented constants with published ones, which is `GROUNDING.md`'s
stated job, and it stops being inert the moment an arm that must reach a lesion
is added — a remyelination agent, for instance.

`validated=False`. Transcribing a model is not reproducing it; no figure of the
paper has been checked.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

# --------------------------------------------------------------------------- #
# Constants transcribed from pcbi.1007117.s001.R. Line numbers are that file's.
# --------------------------------------------------------------------------- #

# Permeability-surface products. PSb and PSc scale with brain volume in the
# source; PSe is a flat number.
#   L565: PSbmale = Vbrain/(1.36/1.04) * 1.875      BBB
#   L569: PScmale = Vbrain/(1.36/1.04) * 0.9375     BCSFB
#   L567: PSemale = 300                             brain mass <-> cranial CSF
PSB_PER_BRAIN_L = 1.875 / (1.36 / 1.04)
PSC_PER_BRAIN_L = 0.9375 / (1.36 / 1.04)
PSE_FLAT = 300.0

# Cranial CSF volume is a constant in the source (L187: Vccsfmale = 0.143 L).
VCCSF_L = 0.143

# Spinal CSF is capped at 20% of total CSF (L189), i.e. the cranial:spinal split
# is 80:20 — the source notes this matches the adult value.
SPINAL_CSF_FRACTION_CAP = 0.2

# Brain blood is 5% of brain volume (L186).
BRAIN_BLOOD_FRACTION = 0.05


@dataclass(frozen=True)
class BrainPBPK:
    """One parameterisation of the four-compartment brain model.

    Flows are L/h, volumes L, PS products L/h. `fu_*` are unbound fractions.
    """

    v_brain: float          # total brain volume, L
    q_brain: float          # brain blood flow, L/h
    q_bulk: float           # bulk flow, brain mass -> cranial CSF, L/h
    q_sin: float            # cranial -> spinal CSF, L/h
    q_sout: float           # spinal -> cranial CSF, L/h
    q_csink: float          # cranial CSF absorption, L/h
    q_ssink: float          # spinal CSF absorption, L/h
    fu_bb: float = 1.0      # unbound fraction, brain blood
    fu_bm: float = 1.0      # unbound fraction, brain mass
    fu_ccsf: float = 1.0    # unbound fraction, cranial CSF
    ps_b: float | None = None   # BBB; derived from v_brain when None
    ps_c: float | None = None   # BCSFB; derived from v_brain when None
    ps_e: float = PSE_FLAT
    source: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def psb(self) -> float:
        return self.ps_b if self.ps_b is not None else self.v_brain * PSB_PER_BRAIN_L

    @property
    def psc(self) -> float:
        return self.ps_c if self.ps_c is not None else self.v_brain * PSC_PER_BRAIN_L

    @property
    def v_bb(self) -> float:
        return BRAIN_BLOOD_FRACTION * self.v_brain

    @property
    def v_ccsf(self) -> float:
        return VCCSF_L

    @property
    def v_scsf(self) -> float:
        return VCCSF_L / (1.0 - SPINAL_CSF_FRACTION_CAP) * SPINAL_CSF_FRACTION_CAP

    @property
    def v_bm(self) -> float:
        rest = self.v_brain - self.v_bb - self.v_ccsf - self.v_scsf
        if rest <= 0.0:
            raise ValueError(
                f"brain mass volume is {rest:.4f} L; v_brain={self.v_brain} is too "
                "small for the CSF and blood volumes this model assumes")
        return rest


# Adult parameterisation. THE VOLUMES AND FLOWS ARE NOT THE PAPER'S — the paper
# is pediatric and its brain-flow expression is an age function fitted over
# 0.25-12.75 y. Standard adult reference values are substituted instead of
# evaluating that curve outside its range. The PS products and the CSF structure
# ARE the paper's.
DEFAULT_ADULT = BrainPBPK(
    v_brain=1.45,     # adult brain ~1.45 L
    q_brain=42.0,     # ~700 mL/min cerebral blood flow
    q_bulk=0.0105,    # ~0.175 mL/min ISF bulk flow
    q_sin=0.021,      # CSF production ~0.35 mL/min = 21 mL/h
    q_sout=0.0105,
    q_csink=0.0105,
    q_ssink=0.0105,
    source=("structure and PS products: Verscheijden 2019 (PMC6592555) S1 File; "
            "adult volumes and flows substituted by this repo because the "
            "paper's cohort and its brain-flow age function are pediatric"),
    meta={"validated": False, "pediatric_source": True, "adult_values_substituted": True},
)


def antibody_psb(cns_serum_ratio: float, q_efflux: float) -> float:
    """BBB permeability-surface product for an ANTIBODY, from a measured ratio.

    The paper's `PSb` describes a passively-diffusing small molecule. An
    antibody does not passively diffuse, so that value must not be reused.
    `bricks/barrier.py` already carries the number to use instead: Pardridge
    2019, CNS:serum ~0.1-0.2% for therapeutic mAbs.

    The conversion, derived from this model's own brain-mass balance at steady
    state. Entry is `PSb * Cbb`; the only net exit from the brain block is CSF
    turnover, `q_efflux * Cbm` (`q_bulk` here, because `PSe` = 300 L/h couples
    brain mass to cranial CSF so tightly that the two move as one pool). So

        Cbm / Cbb  =  PSb / (PSb + q_efflux)      =>      PSb = r*q_efflux/(1-r)

    **The competing clearance is EFFLUX, not perfusion.** A first version of
    this function used `q_brain` (~42 L/h) instead of `q_bulk` (~0.0105 L/h) and
    was wrong by about 4,000x — it returned a PSb that produced an 80% brain:
    plasma ratio for an antibody. Perfusion sets how fast drug is PRESENTED to
    the barrier; it is not what the drug leaving the brain has to compete with.
    """
    if not 0.0 < cns_serum_ratio < 1.0:
        raise ValueError(
            f"cns_serum_ratio must be a fraction in (0,1), got {cns_serum_ratio}. "
            "Pardridge 2019 gives ~0.001-0.002 for therapeutic antibodies.")
    if q_efflux <= 0.0:
        raise ValueError(
            "q_efflux must be positive; with no route out of the brain this model "
            "equilibrates to plasma for any PSb and no ratio is expressible")
    return cns_serum_ratio * q_efflux / (1.0 - cns_serum_ratio)


def simulate(model: BrainPBPK = DEFAULT_ADULT, *, c_arterial: float = 1.0,
             hours: float = 48.0, n: int = 200,
             plasma_half_life_h: float | None = None) -> dict:
    """Integrate the four brain compartments against an arterial input.

    **A CONSTANT ARTERIAL INPUT CANNOT ANSWER THE PENETRATION QUESTION, and
    finding that out is the useful part of this port.** This sub-model has no
    brain-side elimination — the only ways out of brain mass are bulk flow
    (`q_bulk`, ~0.01 L/h) and back-flux across the same barriers. So under a
    constant arterial concentration EVERY compartment equilibrates to that
    concentration whatever `PSb` is: the permeability sets how FAST equilibrium
    arrives, not the level. Measured: a small molecule reaches 0.993 of arterial
    at 48 h and an antibody with `PSb` 28,000x smaller still reaches 0.804.
    Reporting that 0.804 as "CNS penetration" would be flatly wrong.

    The published ~0.15% CNS:serum for a therapeutic antibody is therefore NOT a
    steady-state partition coefficient. It is an exposure ratio under a plasma
    concentration that is itself falling, with transport too slow to keep up.
    Reproducing it needs the plasma side, which is why the paper couples this
    brain model to a whole-body PBPK rather than running it alone.

    So: pass `plasma_half_life_h` to drive the brain with a declining
    mono-exponential plasma profile, and read `auc_ratio`, which is
    PS-dependent and is the quantity that means "penetration". Leaving it None
    keeps the constant-input case available, and its `effective` is labelled
    `equilibrium_fraction` rather than anything suggesting penetration.
    """
    p = model
    psb, psc, pse = p.psb, p.psc, p.ps_e
    v_bb, v_bm, v_ccsf, v_scsf = p.v_bb, p.v_bm, p.v_ccsf, p.v_scsf

    k_el = (np.log(2.0) / plasma_half_life_h) if plasma_half_life_h else 0.0

    def c_art(t):
        return c_arterial * np.exp(-k_el * t) if k_el else c_arterial

    def rhs(_t, y):
        c_bb, c_bm, c_ccsf, c_scsf = y
        d_bb = (p.q_brain * (c_art(_t) - c_bb)
                + psb * (p.fu_bm * c_bm - p.fu_bb * c_bb)
                + psc * (p.fu_ccsf * c_ccsf - p.fu_bb * c_bb)
                + p.q_ssink * c_scsf + p.q_csink * c_ccsf) / v_bb
        d_bm = (psb * (p.fu_bb * c_bb - p.fu_bm * c_bm)
                + pse * (p.fu_ccsf * c_ccsf - p.fu_bm * c_bm)
                - p.q_bulk * c_bm) / v_bm
        d_ccsf = (pse * (p.fu_bm * c_bm - p.fu_ccsf * c_ccsf)
                  + psc * (p.fu_bb * c_bb - p.fu_ccsf * c_ccsf)
                  + p.q_bulk * c_bm + p.q_sout * c_scsf
                  - p.q_sin * c_ccsf - p.q_csink * c_ccsf) / v_ccsf
        d_scsf = (p.q_sin * c_ccsf - p.q_sout * c_scsf - p.q_ssink * c_scsf) / v_scsf
        return [d_bb, d_bm, d_ccsf, d_scsf]

    sol = solve_ivp(rhs, (0.0, hours), [0.0, 0.0, 0.0, 0.0],
                    t_eval=np.linspace(0.0, hours, n), method="LSODA",
                    rtol=1e-8, atol=1e-12)
    if not sol.success:
        raise RuntimeError(f"brain PBPK integration failed: {sol.message}")

    c_bb, c_bm, c_ccsf, c_scsf = sol.y
    plasma = np.array([c_art(t) for t in sol.t])
    auc_brain = float(np.trapezoid(c_bm, sol.t))
    auc_plasma = float(np.trapezoid(plasma, sol.t))
    return {
        "t": sol.t,
        "plasma": plasma,
        "brain_blood": c_bb,
        "brain_mass": c_bm,
        "cranial_csf": c_ccsf,
        "spinal_csf": c_scsf,
        # NOT a penetration fraction. Under a constant input this tends to 1 for
        # any PS, because the sub-model has no brain-side elimination.
        "equilibrium_fraction": float(c_bm[-1] / c_arterial) if c_arterial else 0.0,
        # The PS-dependent quantity, and the one that means "penetration".
        # Only meaningful when the plasma profile declines.
        "auc_ratio": (auc_brain / auc_plasma) if auc_plasma > 0 else float("nan"),
        "plasma_half_life_h": plasma_half_life_h,
        "csf_fraction": float(c_ccsf[-1] / c_arterial) if c_arterial else 0.0,
        "validated": False,
        "source": p.source,
        "note": ("Verscheijden 2019 four-compartment brain model, transcribed from "
                 "the authors' S1 R file. Pediatric, small-molecule "
                 "parameterisation; see the module docstring before using PSb "
                 "for an antibody."),
    }


if __name__ == "__main__":
    print("B6 — Verscheijden 2019 brain PBPK (PMC6592555), ported.\n")

    ab_psb = antibody_psb(0.0015, DEFAULT_ADULT.q_bulk)
    antibody = BrainPBPK(**{**DEFAULT_ADULT.__dict__,
                            "ps_b": ab_psb, "ps_c": ab_psb / 2.0})

    print("CONSTANT arterial input — cannot express penetration, and here is why:")
    print(f"{'':22} {'equilibrium fraction at 48h':>28}")
    for label, m in (("small molecule", DEFAULT_ADULT), ("antibody", antibody)):
        print(f"  {label:<20} {simulate(m)['equilibrium_fraction']:>28.4f}")
    print("  The competing clearance is CSF turnover (~0.01 L/h), not perfusion")
    print("  (~42 L/h) — getting that wrong overstates antibody entry ~4,000x.\n")

    print("DECLINING plasma — the boundary condition the quantity actually needs:")
    print(f"{'':22} {'t1/2 (h)':>10} {'AUC brain / AUC plasma':>24}")
    for label, m, thalf in (("small molecule", DEFAULT_ADULT, 2.5),
                            ("antibody", antibody, 480.0)):
        r = simulate(m, hours=4 * thalf, plasma_half_life_h=thalf)
        print(f"  {label:<20} {thalf:>10.1f} {r['auc_ratio']:>24.4f}")

    print("\n  Read `auc_ratio`, never `equilibrium_fraction`, for anything that")
    print("  means CNS penetration. The published ~0.15% CNS:serum for a")
    print("  therapeutic antibody is an exposure ratio under falling plasma, not a")
    print("  partition coefficient — which is why the paper couples this brain")
    print("  model to a whole-body PBPK instead of running it alone.")
    print("  (ported, NOT reproduced — validated=False)")
