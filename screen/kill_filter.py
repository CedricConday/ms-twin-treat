"""Kill doomed candidates. Four filters, none of which needs a predicted effect.

A candidate here is a `MechanismProfile` — a set of multipliers on the Vélez
model's named rates. The screen asks four questions, cheapest first, and each one
can only ever say "this is doomed", never "this is good".

  1. OUT OF REGIME        the model runs away on it, so nothing it produces means
                          anything (bricks/qsp_velez.py has no carrying capacity
                          on E)
  2. UNREACHABLE          the dial cannot produce a benefit at ANY potency. Five
                          of the repo's own arms sit on such dials, and four came
                          back OUT OF RANGE in backtest/potency.py for this
                          reason.
  3. DEGENERATE           indistinguishable from a drug that already exists, so
                          even if it worked nothing would have been learned
  4. REGULATORY LIABILITY its target is Treg-biased (bricks/harm_channel.py). A
                          SOFT flag, not a kill — the channel's own statistic is
                          p = 0.067 on the six arm-set targets, or 0.048 after a
                          pre-registered enlargement to seven, and it does not
                          get to veto anything on that evidence. Nominally
                          crossing 0.05 by adding one target does not make a
                          post-hoc hypothesis a filter.

WHAT THIS CANNOT DO, AND WHY THAT IS ENFORCED IN CODE
------------------------------------------------------
It cannot tell you a surviving candidate is good, or order two survivors. That
needs the model to generalise to an unseen mechanism, and `backtest/lomo.py`
measures 45.9pp against a 12.3pp null — worse than answering with the training
arms' average. `rank_candidates()` therefore raises rather than returning a
sorted list, because a ranked list is the artifact somebody would screenshot.

SURVIVING IS NOT PASSING. A candidate that clears all four filters has only
avoided the failures this model can see. The largest known blind spot is that
the entire depleting / sequestering / trafficking class cannot come out
beneficial here at all — not a calibration gap, a structural property (see
bricks/qsp_velez.py). A screen over this model is blind to most modern
high-efficacy DMT mechanisms.

`validated=False` throughout.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from itertools import combinations

import numpy as np

from bricks.harm_channel import ARM_TARGET
from bricks.harm_channel import score as harm_score
from bricks.profiles import PROFILES, touched_points
from bricks.qsp_velez import INTERVENTION_POINTS, MechanismProfile, simulate

# Seeds for the regime and reachability checks. Small on purpose: these ask
# yes/no structural questions ("does it diverge", "can it ever help"), not
# "how much", so they do not need the 128-seed cohort backtest/lomo.py needs to
# estimate a magnitude through a ~90,000-fold spread.
SCREEN_SEEDS = (0, 1, 2, 3, 4, 5, 6, 7)
SCREEN_HORIZON_DAYS = 730.0

# Potencies probed when asking whether a dial can EVER help.
REACHABILITY_PROBES = (0.2, 0.5, 0.8)

# A candidate must beat untreated damage by more than the noise floor to count
# as reachable. 14% is the measured bootstrapped CV of the median at n=128
# (scripts/measure_qsp_variance.py); anything smaller is not a result.
NOISE_FLOOR = 0.14


class KillReason(Enum):
    OUT_OF_REGIME = "the model runs away on this candidate"
    UNREACHABLE = "this dial cannot produce a benefit at any potency"
    DEGENERATE = "indistinguishable from a drug that already exists"


@dataclass(frozen=True)
class ScreenResult:
    """One candidate's verdict, with everything needed to disagree with it."""

    profile: MechanismProfile
    killed_by: KillReason | None
    detail: str = ""
    regulatory_liability: float | None = None   # Treg:effector ratio, if known
    like_existing: tuple[str, ...] = ()
    median_damage: float | None = None        # at the candidate's own potency
    best_damage: float | None = None          # best over REACHABILITY_PROBES
    best_potency: float | None = None         # the potency that achieved it
    untreated_damage: float | None = None

    @property
    def survived(self) -> bool:
        return self.killed_by is None

    @property
    def validated(self) -> bool:
        return False


def _median_damage(profile: MechanismProfile) -> float | None:
    vals = []
    for seed in SCREEN_SEEDS:
        traj = simulate(profile, t_end=SCREEN_HORIZON_DAYS, seed=seed)
        if not traj["in_regime"]:
            return None
        vals.append(float(traj["total_damage"][-1]))
    return float(np.median(vals))


def enumerate_candidates(max_points: int = 2,
                         potency: float = 0.5) -> list[MechanismProfile]:
    """Every combination of up to `max_points` dials, each turned both ways.

    Deliberately exhaustive rather than clever. The candidate space this model
    admits is small — seven dials — and a search heuristic over a space you can
    enumerate is a way to hide which corners were never visited.
    """
    out = []
    for k in range(1, max_points + 1):
        for points in combinations(INTERVENTION_POINTS, k):
            for directions in range(2 ** k):
                mults = {}
                for i, p in enumerate(points):
                    down = bool(directions & (1 << i))
                    mults[p] = (1.0 - potency) if down else (1.0 + potency)
                label = ",".join(f"{p}{'-' if mults[p] < 1 else '+'}" for p in points)
                out.append(MechanismProfile(label=label, source="enumerated candidate",
                                            **mults))
    return out


def screen(candidates: list[MechanismProfile],
           target_gene: dict[str, str] | None = None) -> list[ScreenResult]:
    """Run the four filters. Returns a verdict per candidate, in input order."""
    target_gene = target_gene or {}
    untreated = _median_damage(PROFILES["untreated"])
    if untreated is None:
        raise RuntimeError("the untreated arm left the model's regime; nothing to "
                           "screen against")

    existing = {touched_points(p): name for name, p in PROFILES.items()
                if name != "untreated"}

    results = []
    for cand in candidates:
        points = touched_points(cand)

        dmg = _median_damage(cand)
        if dmg is None:
            results.append(ScreenResult(
                cand, KillReason.OUT_OF_REGIME,
                "at least one infection history diverged; damage is undefined",
                untreated_damage=untreated))
            continue

        # Reachability: can this dial pattern help at ANY probed potency? The
        # BEST probed potency is what decides, and it is carried into the result
        # -- reporting the candidate's default-potency damage while it passed at
        # a different one would show a number that did not earn the pass.
        best, best_at = dmg, None
        for probe in REACHABILITY_PROBES:
            probed = {p: (1.0 - probe) if getattr(cand, p) < 1.0 else (1.0 + probe)
                      for p in points}
            d = _median_damage(MechanismProfile(label=cand.label, source="probe",
                                                **probed))
            if d is not None and d < best:
                best, best_at = d, probe
        if best > untreated * (1.0 - NOISE_FLOOR):
            results.append(ScreenResult(
                cand, KillReason.UNREACHABLE,
                f"best probed damage {best:.4f} vs untreated {untreated:.4f} — no "
                f"improvement beyond the {NOISE_FLOOR:.0%} noise floor at any potency",
                median_damage=dmg, best_damage=best, best_potency=best_at,
                untreated_damage=untreated))
            continue

        if points in existing:
            results.append(ScreenResult(
                cand, KillReason.DEGENERATE,
                f"same intervention points as {existing[points]}",
                like_existing=(existing[points],),
                median_damage=dmg, best_damage=best, best_potency=best_at,
                untreated_damage=untreated))
            continue

        liability = None
        gene = target_gene.get(cand.label)
        if gene:
            for arm, g in ARM_TARGET.items():
                if g == gene:
                    s = harm_score(arm)
                    liability = s.treg_ratio
                    break

        results.append(ScreenResult(cand, None, "survived all four filters",
                                    regulatory_liability=liability,
                                    median_damage=dmg, best_damage=best,
                                    best_potency=best_at,
                                    untreated_damage=untreated))
    return results


def rank_candidates(*_args, **_kwargs):
    """Refuses. Ranking needs the model to generalise, and it does not.

    This function exists so the refusal is in the code rather than only in a
    docstring somebody can skip. `backtest/lomo.py` is the gate: when its MAE
    drops below its predict-the-mean null, delete this and implement ranking.
    """
    raise NotImplementedError(
        "Ranking surviving candidates by predicted effect is not supported, and "
        "this is deliberate. backtest/lomo.py reports out-of-sample MAE 45.9pp "
        "against a 12.3pp predict-the-mean null: on a mechanism it has not seen, "
        "this model answers with roughly the training arms' average whatever it "
        "is asked about. A ranking built on that is a sorted list of noise. "
        "Kill filters do not need a magnitude and are implemented; ranking does, "
        "and is gated on LOMO beating its null. Measured 2026-09-20, the one "
        "structural fix that flips the depletion sign does NOT lift the gate "
        "either: backtest/lomo_capacity.py fits a binding carrying capacity "
        "inside each fold and returns 45.6pp against the same 12.3pp null."
    )


def main() -> int:
    cands = enumerate_candidates(max_points=2, potency=0.5)
    results = screen(cands)

    kills = {r: 0 for r in KillReason}
    for res in results:
        if res.killed_by:
            kills[res.killed_by] += 1
    survivors = [r for r in results if r.survived]

    print("SCREEN — kill filter over the Velez model's intervention space")
    print(f"  {len(cands)} candidates: every combination of up to 2 dials, both "
          "directions\n")
    for reason, n in kills.items():
        print(f"  killed {n:>4}  {reason.value}")
    print(f"  survived {len(survivors):>3}\n")

    if survivors:
        print("Survivors (surviving is NOT passing — see the module docstring).")
        print("  'damage' is the BEST probed potency, which is what cleared the")
        print("  reachability filter; negative means LESS damage than untreated.\n")
        print(f"  {'candidate':<26} {'potency':>8} {'damage':>10} {'vs untreated':>13}")
        for r in sorted(survivors, key=lambda x: (x.best_damage or 0.0)):
            change = r.best_damage / r.untreated_damage - 1.0
            pot = f"{r.best_potency:.1f}" if r.best_potency is not None else "default"
            print(f"  {r.profile.label:<26} {pot:>8} {r.best_damage:>10.4f} "
                  f"{change:>12.0%}")

    print("\n  No effect size is predicted and no ordering is claimed.")
    print("  rank_candidates() raises: backtest/lomo.py reports 45.9pp against a")
    print("  12.3pp null, so a predicted magnitude would be noise.")
    print("\n  Largest known blind spot: the depleting / sequestering / trafficking")
    print("  class cannot come out beneficial in this model at all. A screen over")
    print("  it is blind to most modern high-efficacy DMT mechanisms.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
