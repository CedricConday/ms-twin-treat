"""Probe A, second half — each arm as multipliers on the four rates of Jenner 2026.

The rule is bricks/profiles.py's, unchanged: WHICH rate a drug moves and in
WHICH direction comes from the drug's pharmacology and never from its trial
outcome; HOW MUCH is one shared stub (`bricks.profiles.STUB_MAGNITUDE`) that
the exam sweeps. Read that module's docstring first; this one only says what
changes when the model underneath has four rates and no cells.

THE FOUR DIALS, in the paper's words (section 2)
------------------------------------------------
    r      remyelination rate
    phi    disease strength: "a proxy for the strength of these cells
           [microglia, macrophages, T cells, B cells] in their ability to
           cause damage"
    eta    the myelin level at which inflammation ramps up or dampens; the
           patient's "resilience" to inflammatory events
    delta  decay rate of inflammation, "the average duration of an
           inflammatory episode"

THIS MODEL'S LIMIT, stated once so it need not be repeated per arm
------------------------------------------------------------------
There is no compartment, no cell type, no regulator and no antigen step. Any
drug that lowers the number or potency of damage-causing immune cells, by any
route, is "phi down": depletion (anti-CD20, anti-CD52, cladribine),
sequestration (S1P modulators), transit block (natalizumab), anti-
proliferation (teriflunomide, dimethyl fumarate), cytokine blockade
(ustekinumab), antigen displacement (glatiramer), effector-clone deletion
(Tovaxin) and type I interferon all land on the same dial with the same
sign. Nineteen arms share one pattern. That is the probe's ceiling and the
reason it can only answer the SIGN question (docs/FOUR_DAY_PLAN.md, day 1),
never a within-dial ranking. `LUMPED` names every arm this applies to.

Two arms touch no rate at all, and the reasons are structural, not
outcome-driven; see `UNASSIGNABLE`. They are kept in the arm set so the exam
scores them (a model that has no handle on a drug predicts "nothing", and is
graded on it) rather than dropped the way bricks/profiles.py drops laquinimod.

The pattern key carries the SIGN ("phi-", "phi+"), unlike bricks/profiles.py's
`touched_points`. There, opposite directions on one rate never coincided; here
both signs of phi are in the arm set, and the exam builds one response column
per pattern from its first arm, so an unsigned key would score the immune
challenges on a suppression curve.
"""

from __future__ import annotations

from bricks.profiles import ENHANCE, STUB_MAGNITUDE, SUPPRESS
from bricks.qsp_minimal import RATES, MinimalProfile

__all__ = ["PROFILES", "LUMPED", "UNASSIGNABLE", "STUB_MAGNITUDE",
           "touched_points", "profile_at"]


def _p(name: str, source: str, **rates: float) -> MinimalProfile:
    unknown = set(rates) - set(RATES)
    if unknown:
        raise ValueError(f"{name}: not rates of this model: {sorted(unknown)}")
    return MinimalProfile(label=name, source=source, **rates)


# Arms whose distinct mechanisms this model collapses onto "phi down". Keyed by
# arm so a report can name them without re-deriving the reason.
LUMPED: dict[str, str] = {
    "IFN-beta": "less T-cell activation and proliferation -> fewer damaging cells -> phi",
    "IFN-beta-1b": "as IFN-beta",
    "glatiramer acetate": "antigen displacement and Treg restoration -> phi; no antigen step, no Treg",
    "natalizumab": "BBB transit block -> phi; no CNS compartment",
    "fingolimod": "lymph-node sequestration -> phi; no lymph-node compartment",
    "ponesimod": "as fingolimod",
    "ozanimod": "as fingolimod",
    "teriflunomide": "proliferation block -> phi",
    "dimethyl fumarate": "reduced effector expansion -> phi",
    "ocrelizumab": "B-cell depletion -> phi; no B cells, and no ke to carry Martinez-Pasamar's number",
    "ofatumumab": "as ocrelizumab",
    "ublituximab": "as ocrelizumab",
    "rituximab": "as ocrelizumab",
    "alemtuzumab": "pan-lymphocyte depletion -> phi",
    "cladribine": "lymphocyte depletion -> phi",
    "atacicept": "BAFF/APRIL block, plasma-cell depletion -> phi; regulatory B-cell loss not representable",
    "ustekinumab": "IL-12/23 block on Th1/Th17 differentiation -> phi",
    "Tovaxin": "deletion of myelin-reactive effector clones -> phi",
}

# Arms this model has NO rate for. Both are left at identity, and the exam
# scores the resulting prediction of "no change". Neither reason mentions the
# trial.
UNASSIGNABLE: dict[str, str] = {
    "abatacept": (
        "CTLA4-Ig blocks CD28 costimulation, i.e. the PRIMING of naive T cells. "
        "bricks/profiles.py puts it on Velez's activation rate `delta`. This model "
        "has no naive pool and no activation step: phi is the damage capacity of "
        "cells that are already effectors. The nearest reading, 'fewer new "
        "effectors eventually means less phi', is a claim about the timescale of "
        "effector turnover that these equations do not contain. No rate."
    ),
    "daclizumab": (
        "Anti-CD25. On the only immune dial this model has, the drug's two "
        "target-level effects pull in opposite directions: fewer regulatory T cells "
        "(phi UP, less suppression) and less IL-2-driven expansion of activated "
        "effectors (phi DOWN, its transplant-rejection rationale). bricks/profiles.py "
        "could put the Treg effect on its own axis; here both land on phi with "
        "opposite signs and no representable magnitude for either. Asserting the "
        "net sign would be a coin flip; the same rule left teriflunomide's Treg "
        "axis off in bricks/profiles.py. No rate."
    ),
}


UNTREATED = _p("untreated", "control arm")

# ---- phi DOWN: every lumped arm. Sources are the ones in bricks/profiles.py;
# only the landing dial differs, for the reason in LUMPED.
def _phi_down(name: str) -> MinimalProfile:
    return _p(name, f"see bricks/profiles.py for the pharmacology; here: {LUMPED[name]}",
              phi=SUPPRESS)


PHI_DOWN = {name: _phi_down(name) for name in LUMPED}

# ---- phi UP: the immunogenic challenges.
IFN_GAMMA = _p(
    "IFN-gamma",
    "Recombinant type II interferon; activates macrophages and upregulates MHC "
    "class II (Panitch 1987, PMID 2882294). Macrophages are among the cells the "
    "paper names as what phi stands for, so more of their activity is phi UP.",
    phi=ENHANCE,
)

APL_CGP77116 = _p(
    "APL CGP77116",
    "Altered peptide ligand of MBP 83-99, encephalitogenic in T-cell assays "
    "(Bielekova 2000, doi:10.1038/80516): it expands the autoreactive effector "
    "pool, which is phi UP.",
    phi=ENHANCE,
)

# ---- Two dials: TNF blockade.
LENERCEPT = _p(
    "lenercept",
    "TNF receptor p55-IgG fusion protein; neutralizes TNF. TNF is an effector "
    "cytokine, so the immune attack weakens: phi DOWN. TNFR2 signalling is "
    "required for oligodendrocyte-progenitor proliferation and remyelination "
    "(Arnett 2001, PMID 11600888), so blocking it slows repair: r DOWN. Both are "
    "animal-model mechanism papers, as in bricks/profiles.py. The only arm whose "
    "outcome in this model is a contest between two dials.",
    phi=SUPPRESS, r=SUPPRESS,
)

# ---- No rate.
ABATACEPT = _p("abatacept", UNASSIGNABLE["abatacept"])
DACLIZUMAB = _p("daclizumab", UNASSIGNABLE["daclizumab"])


PROFILES: dict[str, MinimalProfile] = {p.label: p for p in (
    UNTREATED, *PHI_DOWN.values(), IFN_GAMMA, APL_CGP77116, LENERCEPT, ABATACEPT, DACLIZUMAB,
)}


def touched_points(profile: MinimalProfile) -> tuple[str, ...]:
    """Signed pattern: which rates the arm moves and which way, e.g. ('r-', 'phi-')."""
    out = []
    for rate in RATES:
        v = getattr(profile, rate)
        if v < 1.0:
            out.append(f"{rate}-")
        elif v > 1.0:
            out.append(f"{rate}+")
    return tuple(out)


def profile_at(arm: str, s: float) -> MinimalProfile:
    """The arm's profile at potency `s`: (1 - s) on suppressed rates, (1 + s) on enhanced.

    The same convention as backtest/lomo._profile_at; direction is the sourced
    one and only the size moves.
    """
    base = PROFILES[arm]
    rates = {}
    for rate in RATES:
        v = getattr(base, rate)
        if v != 1.0:
            rates[rate] = (1.0 - s) if v < 1.0 else (1.0 + s)
    return MinimalProfile(label=f"{arm}@{s:g}", source="exam v2 potency sweep", **rates)


if __name__ == "__main__":
    width = max(len(n) for n in PROFILES)
    for name, prof in PROFILES.items():
        tag = "  [LUMPED]" if name in LUMPED else ("  [NO RATE]" if name in UNASSIGNABLE else "")
        print(f"  {name:<{width}}  {touched_points(prof) or '-'}{tag}")
    patterns = {touched_points(p) for p in PROFILES.values()}
    print(f"\n  {len(PROFILES) - 1} arms, {len(patterns) - 1} treated patterns")
