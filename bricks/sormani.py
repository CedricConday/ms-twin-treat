"""Blocker (2) — the micro→clinical map, from a published meta-regression.

`bricks/readout.py` maps simulated damage to a relapse proxy through
`MAX_RELAPSE = 1.5`, a number its own docstring calls illustrative. BUILD_PLAN
§8 blocker (2) records the catch that made replacing it awkward: the obvious
published map (Kotelnikova 2017) goes damage → **EDSS**, while the clinical gate
scores **relapse rate**, so using it would have meant either adding an EDSS
endpoint or calibrating two things at once.

This module uses a different published map that does not have that problem.

    Sormani MP, Bruzzi P. "MRI lesions as a surrogate for relapses in multiple
    sclerosis: a meta-analysis of randomised trials."
    Lancet Neurology 2013;12(7):669-76. PMID 23743084.
    doi:10.1016/S1474-4422(13)70103-0

    31 trials, 18,901 patients with RRMS. Validating and extending
    Sormani et al., Ann Neurol 2009;65:268-75 (PMID 19334061; 23 trials,
    63 arms, 40 contrasts, 6,591 patients, adjusted R^2 = 0.81).

WHAT THE REGRESSION IS
----------------------
Both papers regress a **treatment effect** on a **treatment effect**, at the
level of a trial arm, on the log scale. From the 2013 methods:

    "We extracted data for the treatment effects on MRI lesions and on relapses
     from each trial, and the correlation of log transformed relative measures
     of these treatment effects was assessed with a weighted linear regression
     analysis" (weighted on trial size and duration)

and the result:

    "The regression equation ... showed a relation between the concurrent
     treatment effects on MRI lesions and relapses (slope=0.52; R^2=0.71), much
     the same as was previously estimated (p_interaction=0.45)"

So, with RR meaning a treated-vs-comparator rate ratio:

    log(RR_relapse) = intercept + SLOPE * log(RR_lesion)

which is, equivalently, `RR_relapse = exp(intercept) * RR_lesion ** SLOPE`.

Note what this is NOT. It does not map an absolute lesion count to an absolute
relapse rate, and it is a **trial-level** relation, not an individual-level one
— the 2009 paper is explicit that it explains between-trial variance. So it
belongs at the point where the gate compares two arms, never inside a
per-patient readout.

THE INTERCEPT IS NOT PUBLISHED WHERE THIS REPO CAN READ IT
-----------------------------------------------------------
The slope is in the abstract of a paywalled paper; the intercept is not, and
neither is the 2009 equation's. Rather than guess a number and let it look
sourced, `INTERCEPT` defaults to the value the relation's own boundary
condition forces:

    a treatment with NO effect on lesions (RR_lesion = 1, log = 0) must have
    NO effect on relapses (RR_relapse = 1, log = 0), which requires intercept = 0.

That is an assumption, not a citation, and it is the one number in this module
that someone with the PDF should replace. A fitted intercept would mean the
average trial in the meta-analysis showed a relapse effect not accounted for by
its lesion effect — plausible in real data, and it would shift every prediction
by a constant factor `exp(intercept)`. `SLOPE` is the published value and is
pinned by a test; `INTERCEPT` is flagged as assumed in every output.

WHAT IT DOES NOT FIX — read before quoting a prediction
--------------------------------------------------------
**Lenercept.** Its trial reported no significant MRI difference while the
relapse rate rose significantly (p=0.006-0.007) and relapses were more severe
and longer (Neurology 1999;53:457, PMID 10449104). Put RR_lesion ~ 1 through
this map and it returns RR_relapse ~ 1 — "no effect" — for a drug that harmed
people. The relation is a trial-level average and lenercept is a documented
outlier to it. Any use of this map on a harm candidate inherits that blind spot,
which is why `predict_relapse_ratio` reports `blind_spot=True` when the lesion
effect is small enough for it to bite.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Sormani & Bruzzi 2013, Lancet Neurol 12(7):669-76, PMID 23743084.
SLOPE = 0.52
R_SQUARED = 0.71
N_TRIALS = 31
N_PATIENTS = 18901

# NOT published where this repo can read it. Forced to 0 by the boundary
# condition RR_lesion = 1 -> RR_relapse = 1. See the module docstring.
INTERCEPT = 0.0
INTERCEPT_IS_ASSUMED = True

# Below this |log RR_lesion| the map cannot distinguish a drug that does nothing
# from one that harms without touching MRI. Chosen as a 10% lesion effect, which
# is where the lenercept case sits: that trial's MRI difference was not
# significant. Stated rather than hidden -- it is a reporting threshold for a
# known blind spot, not a claim about where the regression stops being valid.
BLIND_SPOT_LOG_RR = abs(math.log(0.90))


@dataclass(frozen=True)
class RelapsePrediction:
    """A predicted relapse-rate ratio, with everything needed to distrust it."""

    rr_relapse: float
    rr_lesion: float
    percent_change: float
    blind_spot: bool
    intercept_assumed: bool = INTERCEPT_IS_ASSUMED
    source: str = (
        "Sormani & Bruzzi 2013, Lancet Neurol 12(7):669-76 (PMID 23743084): "
        f"slope={SLOPE}, R^2={R_SQUARED}, {N_TRIALS} trials / {N_PATIENTS} patients. "
        "Intercept assumed 0 (not published where this repo can read it)."
    )

    @property
    def validated(self) -> bool:
        """False, always. The regression is published; this application is not.

        Sormani relates the MRI effect a TRIAL measured to the relapse effect
        that TRIAL measured. Feeding it a lesion ratio produced by a toy
        simulation is an extrapolation the paper says nothing about.
        """
        return False


def predict_relapse_ratio(rr_lesion: float, *, slope: float = SLOPE,
                          intercept: float = INTERCEPT) -> RelapsePrediction:
    """Map a treated-vs-comparator LESION rate ratio to a RELAPSE rate ratio.

    `rr_lesion` is arm-level: lesion measure on treatment divided by the same
    measure on that arm's comparator. 1.0 is no effect, <1 is fewer lesions.

    Returns a `RelapsePrediction` rather than a float on purpose — the blind
    spot and the assumed intercept have to travel with the number, because the
    one case this map gets confidently wrong (lenercept) looks exactly like a
    clean "no effect" result.
    """
    if not rr_lesion > 0.0:
        raise ValueError(f"rr_lesion must be positive, got {rr_lesion}")

    log_rr_lesion = math.log(rr_lesion)
    rr_relapse = math.exp(intercept + slope * log_rr_lesion)
    return RelapsePrediction(
        rr_relapse=rr_relapse,
        rr_lesion=rr_lesion,
        percent_change=(rr_relapse - 1.0) * 100.0,
        blind_spot=abs(log_rr_lesion) < BLIND_SPOT_LOG_RR,
        intercept_assumed=(intercept == INTERCEPT and INTERCEPT_IS_ASSUMED),
    )


def lesion_ratio(treated_lesions: float, comparator_lesions: float) -> float:
    """Arm-level lesion rate ratio, guarded at zero.

    A comparator arm with zero lesions makes the ratio undefined; a treated arm
    with zero lesions makes log(0) = -inf. Both are real possibilities in a
    small simulated cohort and neither should be papered over with a silent
    epsilon, so both raise.
    """
    if comparator_lesions <= 0.0:
        raise ValueError(
            "comparator arm has no lesions — the lesion rate ratio is undefined. "
            "A simulated comparator producing zero lesions means the cohort or the "
            "horizon is too small, not that the drug is infinitely good."
        )
    if treated_lesions < 0.0:
        raise ValueError(f"treated_lesions must be >= 0, got {treated_lesions}")
    if treated_lesions == 0.0:
        raise ValueError(
            "treated arm has no lesions at all — log(0) is undefined and the map "
            "cannot express 'perfect'. Widen the cohort or the horizon."
        )
    return treated_lesions / comparator_lesions


if __name__ == "__main__":
    print("Blocker (2) — Sormani & Bruzzi 2013 trial-level map, lesions -> relapses.")
    print(f"  log(RR_relapse) = {INTERCEPT} + {SLOPE} * log(RR_lesion)")
    print(f"  R^2 = {R_SQUARED} over {N_TRIALS} trials / {N_PATIENTS:,} RRMS patients")
    print("  intercept is ASSUMED 0, not published where this repo can read it\n")

    print(f"{'lesion effect':>14} {'RR_lesion':>10} {'RR_relapse':>11} "
          f"{'relapse effect':>15}  note")
    for pct in (-90, -80, -68, -50, -30, -10, 0, +10, +50):
        rr_les = 1.0 + pct / 100.0
        if rr_les <= 0:
            continue
        pred = predict_relapse_ratio(rr_les)
        note = "BLIND SPOT — cannot separate 'no effect' from harm" if pred.blind_spot else ""
        print(f"{pct:>13}% {rr_les:>10.2f} {pred.rr_relapse:>11.3f} "
              f"{pred.percent_change:>14.1f}%  {note}")

    print("\n  The blind spot is not hypothetical: the lenercept trial reported no")
    print("  significant MRI difference while relapses rose significantly (p=0.006).")
    print("  This map returns 'no effect' for it. (validated=False)")
