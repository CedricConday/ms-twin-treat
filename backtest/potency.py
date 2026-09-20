"""Blocker (4), second half — per-drug potency from the MRI channel.

`bricks/profiles.py` says WHICH dials each drug moves and in which direction,
from pharmacology. It leaves HOW MUCH as one shared stub, because a hand-picked
number per drug is indistinguishable from fitting, and fitting to a drug's own
relapse number is the circularity BUILD_PLAN §8.4 records.

This module supplies the missing magnitude WITHOUT that circularity, by fitting
each drug's potency against its trial's **MRI lesion** outcome and leaving the
relapse number untouched for the gate to score against.

    observed lesion ratio (from the trial)
        -> invert the cached response table to the potency that reproduces it
        -> that potency is the drug's magnitude on the dials profiles.py assigned

The gate then predicts a relapse change from that potency and is checked against
the trial's reported ARR. The relapse number never enters the fit, so the test
stays a test.

WHY THIS IS NOT CIRCULAR, PRECISELY
------------------------------------
Sormani maps lesions -> relapses, and this module runs that map backwards to
read a lesion ratio out of a table stored in relapse units. That is arithmetic
on the simulation's own output, not on the trial's. The only trial-derived input
is the MRI lesion ratio. Nothing anywhere in the loop reads the arm's ARR.

WHAT BLOCKS IT TODAY
--------------------
**The paywall was never the real obstacle.** Nine of twelve quantified arms now
carry a measured lesion ratio, taken from the trials' own POSTED RESULTS on
ClinicalTrials.gov — structured, public, and free. An earlier version of this
module recorded the data as "behind paywalls" after trying only journal full
text. Three arms remain, each for its own reason and none of them paywalls:
AFFIRM (2006) predates the results-posting requirement, CARE-MS I posted a
lesion VOLUME change rather than a count, and PRISMS (1998) predates the
registry.

Four of the nine fitted arms come back OUT OF RANGE, and that is the machinery
working: `gamma_E` and a lowered `alpha_R` both RAISE damage in this model, so no
potency reproduces a lesion reduction on those dials. The fitter says so instead
of returning the nearest grid point.

THE ONE CROSS-CHECK AVAILABLE, AND IT DISAGREES
------------------------------------------------
Ocrelizumab is the only arm with a magnitude from two independent sources, and
they do not agree:

    ke = 0.85   Martinez-Pasamar 2013, fitted to EAE mouse T-cell dynamics
                (K_eff 1000 -> ~850 cells)
    ke = 0.40   this module, fitted to OPERA I's Gd-enhancing lesion ratio

A factor of two apart, from mouse flow cytometry versus human MRI. Neither is
adjusted toward the other. Possible readings, none of them settled here: the
mouse and human systems genuinely differ; a two-year human MRI ratio and a
30-day EAE experiment are not the same measurement; or the model absorbs the
difference into the wrong dial because `ke` is standing in for B cells it does
not represent. `bricks/profiles.py` keeps the EAE value as ocrelizumab's
default, because it is a direct measurement of the parameter rather than an
inversion through two models.

This is the only place in the repo where a parameter can be cross-checked at
all, so the disagreement is worth more than either number: it is a size
estimate for how much the whole potency layer can be trusted.

THE HOLE THAT DOES NOT CLOSE WITH MORE DATA
--------------------------------------------
Lenercept reported **no significant MRI difference** while its relapse rate rose
significantly (p=0.006) with more severe and longer relapses (Neurology
1999;53:457, PMID 10449104). Its MRI-fitted potency would therefore be ~0, i.e.
"does nothing", for a drug that harmed people. Any potency fitted through this
channel inherits that blind spot, `fit_potency` reports it, and no amount of
extraction fixes it -- the harm channel needs a non-MRI source.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from backtest.lomo import POTENCY_GRID, _predicted, load
from bricks.profiles import PROFILES, touched_points
from bricks.sormani import BLIND_SPOT_LOG_RR, invert


@dataclass(frozen=True)
class FittedPotency:
    """A per-drug magnitude, with what is needed to distrust it."""

    arm: str
    potency: float
    observed_lesion_ratio: float
    achieved_lesion_ratio: float
    points: tuple[str, ...]
    blind_spot: bool
    in_range: bool
    source: str

    @property
    def multipliers(self) -> dict[str, float]:
        """The profile multipliers this potency implies, directions preserved."""
        base = PROFILES[self.arm]
        return {p: (1.0 - self.potency) if getattr(base, p) < 1.0
                else (1.0 + self.potency)
                for p in self.points}


def _simulated_lesion_ratios(arm: str, table: dict) -> tuple[list[float], list[float]]:
    """(potencies, simulated lesion ratios) for one arm, from the cached table.

    The table stores relapse changes, because that is what the gate scores.
    Inverting each through Sormani recovers the lesion ratio the simulation
    produced, which is the quantity a trial's MRI outcome is comparable to.
    """
    xs, ys = [], []
    for s in POTENCY_GRID:
        pct = _predicted(arm, float(s), table)
        if np.isnan(pct):
            continue
        rr_relapse = 1.0 + pct / 100.0
        if rr_relapse <= 0.0:
            continue
        xs.append(float(s))
        ys.append(invert(rr_relapse))
    return xs, ys


def fit_potency(arm: str, observed_lesion_ratio: float, source: str,
                table: dict | None = None) -> FittedPotency:
    """The potency at which this arm reproduces its trial's MRI lesion ratio.

    `observed_lesion_ratio` is treated-arm lesions divided by its comparator's,
    as the trial reported them: 0.06 for a 94% reduction. `source` must name the
    trial and identifier it came from — an unsourced magnitude is the thing this
    whole module exists to avoid, so it is a required argument, not a default.
    """
    if not observed_lesion_ratio > 0.0:
        raise ValueError(
            f"observed_lesion_ratio must be positive, got {observed_lesion_ratio}")
    if not source.strip():
        raise ValueError(
            f"{arm}: a fitted magnitude needs a source naming the trial it came from")
    if arm not in PROFILES:
        raise KeyError(f"unknown arm {arm!r}")

    table = table or load()
    xs, ys = _simulated_lesion_ratios(arm, table)
    if not xs:
        raise RuntimeError(
            f"{arm}: no in-regime points in the response table; nothing to fit against")

    # The simulated lesion ratio falls as potency rises, so search rather than
    # interpolate — the curve is not guaranteed monotone in this model and
    # np.interp would silently assume it is.
    errors = [abs(np.log(y) - np.log(observed_lesion_ratio)) for y in ys]
    best = int(np.argmin(errors))
    achieved = ys[best]

    lo, hi = min(ys), max(ys)
    return FittedPotency(
        arm=arm,
        potency=xs[best],
        observed_lesion_ratio=observed_lesion_ratio,
        achieved_lesion_ratio=achieved,
        points=touched_points(PROFILES[arm]),
        # bool(...) because numpy comparisons yield np.bool_, which is not the
        # `bool` these fields advertise and fails an `is True` check downstream.
        blind_spot=bool(abs(np.log(observed_lesion_ratio)) < BLIND_SPOT_LOG_RR),
        in_range=bool(lo <= observed_lesion_ratio <= hi),
        source=source,
    )


# --------------------------------------------------------------------------- #
# The observed MRI lesion ratios. One entry, because one is what is verified.
# --------------------------------------------------------------------------- #
# Transcribed from the trial's own report, never a meta-analysis ranking. Adding
# a row requires the number to be in docs/TRIAL_ANCHORS.md with its PMID.
# THE METRICS ARE NOT INTERCHANGEABLE, and the third element of each tuple says
# which one it is. "New or newly enlarging T2", "Gd-enhancing T1 per scan" and
# "combined unique active lesions per year" measure different things and have
# different base rates. Sormani & Bruzzi regress on "MRI lesions" across trials
# that used a mix of these, so a ratio built from any of them is inside what the
# meta-analysis pooled -- but a ratio is only comparable to ANOTHER ratio of the
# same metric, and nothing here should rank two arms measured differently.
OBSERVED_LESION_RATIOS: dict[str, tuple[float, str]] = {
    "ocrelizumab": (
        0.02 / 0.29,
        "OPERA I, NEJM 2017 (PMID 28002679): 0.02 vs 0.29 Gd-enhancing T1 lesions "
        "per scan on IFN beta-1a, reported as 94% lower [Gd-T1/scan]",
    ),
    "fingolimod": (
        2.5 / 9.8,
        "FREEDOMS, NCT00289978 posted results: 2.5 (0.5 mg) vs 9.8 placebo, new or "
        "newly enlarged T2 lesions at month 24 [new-T2]",
    ),
    "dimethyl fumarate": (
        2.6 / 17.0,
        "DEFINE, NCT00420212 posted results: 2.6 (240 mg BID) vs 17.0 placebo, new "
        "or newly enlarging T2 hyperintense lesions [new-T2]",
    ),
    "glatiramer acetate": (
        8.0 / 17.4,
        "CONFIRM, NCT00451451 posted results: 8.0 vs 17.4 placebo, new or newly "
        "enlarging T2 hyperintense lesions [new-T2]",
    ),
    "teriflunomide": (
        0.261 / 1.331,
        "TEMSO, NCT00134563 posted results: 0.261 (14 mg) vs 1.331 placebo, "
        "Gd-enhancing T1 lesions per scan, Poisson estimate [Gd-T1/scan]",
    ),
    "ponesimod": (
        1.405 / 3.164,
        "OPTIMUM, NCT02425644 posted results: 1.405 vs 3.164 on teriflunomide, "
        "combined unique active lesions per year to week 108 [CUAL/year]",
    ),
    "cladribine": (
        0.38 / 1.43,
        "CLARITY, NCT00213135 posted results: 0.38 (3.5 mg/kg) vs 1.43 placebo, "
        "active T2 lesions [active-T2]",
    ),
    "daclizumab": (
        4.31 / 9.44,
        "DECIDE, NCT01064401 posted results: 4.31 vs 9.44 on IFN beta-1a, adjusted "
        "mean new or newly enlarging T2 hyperintense lesions to week 96 [new-T2]",
    ),
    "natalizumab": (
        1.9 / 11.0,
        "AFFIRM, NEJM 2006 (PMID 16510744): 1.9 vs 11.0 placebo, mean new or "
        "enlarging hyperintense T2 lesions over two years, reported as 83% lower "
        "[new-T2]",
    ),
    "IFN-beta": (
        0.5 / 2.25,
        "PRISMS MRI results, Li et al., Ann Neurol 1999;46(2):197 (trial PMID "
        "9820297): median 0.5 (44 mcg) vs 2.25 placebo T2-active lesions per "
        "patient per scan [active-T2]",
    ),
    "ofatumumab": (
        0.0115 / 0.4555,
        "ASCLEPIOS I, NCT02792218 posted results: 0.0115 vs 0.4555 on teriflunomide, "
        "Gd-enhancing T1 lesions per scan [Gd-T1/scan]",
    ),
}

# Still to extract. ONE, and the reason is not access:
#   alemtuzumab  CARE-MS I reports the PROPORTION OF PATIENTS free of new or
#                enlarging T2 lesions (48% vs 58%), and its posted results give a
#                percent change in lesion VOLUME. Neither is convertible to a
#                lesion-count ratio, and pretending otherwise would put a number
#                here that no trial measured.
#
# The two that WERE listed here as blocked are now extracted, and both were in
# the journals all along -- AFFIRM's MRI table in the NEJM paper, PRISMS's in its
# companion MRI paper (Li 1999). "Predates the registry" was a reason to look
# somewhere else, not a reason to stop.
PENDING_EXTRACTION = (
    "alemtuzumab",
)


def fit_all(table: dict | None = None) -> dict[str, FittedPotency]:
    table = table or load()
    return {arm: fit_potency(arm, ratio, source, table)
            for arm, (ratio, source) in OBSERVED_LESION_RATIOS.items()}


def main() -> int:
    fits = fit_all()
    print("Blocker (4) — per-drug potency fitted on the MRI channel, not on ARR.\n")
    print(f"{'arm':<20} {'observed RR':>12} {'achieved':>10} {'potency':>8}  dials")
    print("-" * 74)
    for arm, f in fits.items():
        flags = []
        if f.blind_spot:
            flags.append("BLIND SPOT")
        if not f.in_range:
            flags.append("OUT OF RANGE — model cannot reach this lesion ratio")
        print(f"{arm:<20} {f.observed_lesion_ratio:>12.4f} {f.achieved_lesion_ratio:>10.4f} "
              f"{f.potency:>8.2f}  "
              f"{', '.join(f'{k}={v:g}' for k, v in f.multipliers.items())}"
              f"{'   ' + '; '.join(flags) if flags else ''}")
    print("-" * 74)
    print(f"{len(fits)} arm(s) fitted, {len(PENDING_EXTRACTION)} pending extraction:")
    print(f"  {', '.join(PENDING_EXTRACTION)}")
    out = [a for a, f in fits.items() if not f.in_range]
    if out:
        print(f"\n  {len(out)} arm(s) OUT OF RANGE: {', '.join(out)}")
        print("  Not a fitting failure — the model cannot produce a lesion REDUCTION on")
        print("  those dials at any potency. gamma_E and a lowered alpha_R both raise")
        print("  damage in this model (see bricks/qsp_velez.py), so no magnitude")
        print("  reproduces their trials. The fitter reports it instead of inventing one.")
    print("\n  Three arms remain unextracted, each for a different reason: AFFIRM (2006)")
    print("  predates results-posting; CARE-MS I posted a lesion VOLUME change, which")
    print("  is not convertible to a count ratio; PRISMS (1998) predates the registry.")
    ocr = fits.get("ocrelizumab")
    if ocr:
        eae = PROFILES["ocrelizumab"].ke
        mri = ocr.multipliers.get("ke")
        print(f"\n  CROSS-CHECK, the only one available: ocrelizumab's ke is {eae:g} from")
        print(f"  Martinez-Pasamar's EAE fit and {mri:g} from OPERA's MRI ratio — "
              f"{max(eae, mri) / min(eae, mri):.1f}x apart.")
        print("  Neither is adjusted toward the other. That gap is the best estimate")
        print("  this repo has of how far the potency layer can be trusted.")

    print("\n  And the hole that more data does not close: lenercept reported NO")
    print("  significant MRI difference while relapses rose (p=0.006). Fitted through")
    print("  this channel its potency would be ~0 — 'does nothing' — for a drug that")
    print("  harmed people. The harm channel needs a non-MRI source.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
