"""Is the predictor entitled to be listened to? Measure it, do not assume it.

An `EvidenceCertificate` is the answer to one question, per scorer:

    Fitted on everything except the thing it was asked about, does this
    predictor beat answering with the training arms' average?

It is re-measured, never recalled. `certify()` runs the leave-one-arm-out
scorer live and reads the leave-one-mechanism-out scorer from its cached run
(LOMO costs minutes per mechanism curve; LOO reads a tabulated response curve
and costs milliseconds). A scorer with no measurement is `UNMEASURED`, which is
NOT a pass -- `gate/criterion.py` clause (5).

THE STATISTIC
-------------
Per held-out unit the two scorers already report `error` and `null_error`. The
paired difference `error - null_error` is negative when the predictor beat the
null on that unit. The certificate reports the mean of those differences and a
paired bootstrap CI over them, resampling the units themselves -- which is the
right resampling unit here because the folds are the independent replicates, not
the seeds inside them.

A predictor passes only if the CI lies entirely below zero AND the relative
margin clears `CRITERION.min_relative_margin`. The CI guards against a lucky
inversion on 12 arms; the margin guards against an inversion so small it is
inside the model's own noise floor.

WHAT A PASSING CERTIFICATE WOULD AND WOULD NOT MEAN
----------------------------------------------------
Would: on held-out units, this predictor's errors were smaller than the null's
by more than resampling noise. That is a statement about a simulation scored
against real trial numbers.

Would not: that the simulation is biologically correct, that the arms are a
representative sample of MS therapy, or that a candidate scoring well would work
in a person. The arm set is 12 quantified arms chosen because their trials
report ARR, five of the mechanism dials are a single lumped `gamma_E` (see
`bricks/profiles.py`), and the whole depleting/sequestering class is
structurally unable to come out beneficial on this model. A certificate cannot
certify past those.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from gate.criterion import CRITERION, AcceptanceCriterion

# Where a LOMO run leaves its scored rows for the certificate to read. Written
# by `scripts/cache_lomo_certificate.py`; absent by default, which is the
# UNMEASURED case and is handled rather than crashed on.
LOMO_CACHE = Path(__file__).resolve().parent.parent / "results" / "lomo_certificate.json"

MEASURED = "MEASURED"
UNMEASURED = "UNMEASURED"


def _paired_bootstrap(diffs: list[float], resamples: int, confidence: float,
                      seed: int = 0) -> tuple[float, float]:
    """CI on the mean of paired (error - null_error) differences.

    Resamples the held-out units with replacement. Returns (lo, hi); the
    predictor beat the null beyond noise when hi < 0.
    """
    if len(diffs) < 2:
        # One unit cannot bound its own sampling error. Report an unbounded
        # interval rather than a reassuring one.
        return (float("-inf"), float("inf"))
    arr = np.asarray(diffs, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(arr), size=(resamples, len(arr)))
    means = arr[idx].mean(axis=1)
    alpha = (1.0 - confidence) / 2.0
    return (float(np.quantile(means, alpha)), float(np.quantile(means, 1.0 - alpha)))


@dataclass(frozen=True)
class ScorerEvidence:
    """One out-of-sample scorer's measured performance against its null."""

    name: str
    status: str                      # MEASURED | UNMEASURED
    n_units: int = 0
    mae: float = float("nan")
    null_mae: float = float("nan")
    ci_low: float = float("nan")
    ci_high: float = float("nan")
    mae_subset: float | None = None       # placebo-controlled arms only, LOO
    null_mae_subset: float | None = None
    detail: str = ""

    def beats_null(self, criterion: AcceptanceCriterion = CRITERION) -> bool:
        if self.status != MEASURED:
            return False
        if not criterion.margin_met(self.mae, self.null_mae):
            return False
        if not (self.ci_high < 0):
            return False
        if criterion.require_placebo_subset and self.mae_subset is not None:
            if not criterion.margin_met(self.mae_subset, self.null_mae_subset):
                return False
        return True

    def line(self) -> str:
        if self.status != MEASURED:
            return f"{self.name:<24} UNMEASURED  ({self.detail})"
        rel = ((self.null_mae - self.mae) / self.null_mae * 100) if self.null_mae else float("nan")
        return (f"{self.name:<24} MAE {self.mae:>6.1f}pp   null {self.null_mae:>6.1f}pp   "
                f"margin {rel:>+6.1f}%   CI [{self.ci_low:+.1f}, {self.ci_high:+.1f}]   "
                f"n={self.n_units}")


@dataclass(frozen=True)
class EvidenceCertificate:
    """Everything the device is allowed to know about its own predictor."""

    scorers: dict[str, ScorerEvidence] = field(default_factory=dict)
    criterion: AcceptanceCriterion = CRITERION

    @property
    def unlocks_pass(self) -> bool:
        """May the device emit PASS at all?

        Requires the arm-holdout scorer to clear the bar, and -- because a
        screened candidate is a new mechanism -- the mechanism-holdout scorer
        to be measured and clear it too.
        """
        loo = self.scorers.get("leave-one-arm-out")
        if loo is None or not loo.beats_null(self.criterion):
            return False
        if self.criterion.require_mechanism_holdout:
            lomo = self.scorers.get("leave-one-mechanism-out")
            if lomo is None or not lomo.beats_null(self.criterion):
                return False
        return True

    def blocking_reasons(self) -> list[str]:
        """Why PASS is unavailable, in the order a reader should care."""
        out = []
        for key in ("leave-one-arm-out", "leave-one-mechanism-out"):
            sc = self.scorers.get(key)
            if sc is None:
                out.append(f"{key}: not run")
            elif sc.status != MEASURED:
                out.append(f"{key}: UNMEASURED ({sc.detail})")
            elif not sc.beats_null(self.criterion):
                if sc.mae >= sc.null_mae:
                    out.append(f"{key}: loses to the null ({sc.mae:.1f}pp vs {sc.null_mae:.1f}pp)")
                else:
                    out.append(f"{key}: beats the null but not by the required margin "
                               f"({sc.mae:.1f}pp vs {sc.null_mae:.1f}pp)")
        return out

    def report(self) -> str:
        lines = [self.criterion.describe(), ""]
        for sc in self.scorers.values():
            lines.append("  " + sc.line())
        lines.append("")
        lines.append(f"  PASS unlocked: {self.unlocks_pass}")
        for r in self.blocking_reasons():
            lines.append(f"    blocked by {r}")
        return "\n".join(lines)


def _loo_evidence(criterion: AcceptanceCriterion) -> ScorerEvidence:
    from backtest.loo import run_loo

    r = run_loo()
    diffs = [row["error"] - row["null_error"] for row in r["rows"]]
    lo, hi = _paired_bootstrap(diffs, criterion.bootstrap_resamples,
                               criterion.bootstrap_confidence)
    return ScorerEvidence(
        name="leave-one-arm-out",
        status=MEASURED,
        n_units=len(r["rows"]),
        mae=r["mae"],
        null_mae=r["null_mae"],
        ci_low=lo,
        ci_high=hi,
        mae_subset=r["mae_placebo_only"],
        null_mae_subset=r["null_mae_placebo_only"],
        detail="run live from the cached response curve",
    )


def _lomo_evidence(criterion: AcceptanceCriterion,
                   cache: Path = LOMO_CACHE) -> ScorerEvidence:
    """Read a cached LOMO run. Absent or malformed cache -> UNMEASURED, not a crash."""
    if not cache.exists():
        return ScorerEvidence(
            name="leave-one-mechanism-out", status=UNMEASURED,
            detail=f"no cached run at {cache.name}; "
                   "run scripts/cache_lomo_certificate.py (minutes, not seconds)")
    try:
        payload = json.loads(cache.read_text())
        folds = payload["folds"]
        diffs = [f["mae"] - f["null_mae"] for f in folds]
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        return ScorerEvidence(
            name="leave-one-mechanism-out", status=UNMEASURED,
            detail=f"cache unreadable ({type(exc).__name__}); treated as no measurement")
    lo, hi = _paired_bootstrap(diffs, criterion.bootstrap_resamples,
                               criterion.bootstrap_confidence)
    return ScorerEvidence(
        name="leave-one-mechanism-out",
        status=MEASURED,
        n_units=len(folds),
        mae=float(payload["mae"]),
        null_mae=float(payload["null_mae"]),
        ci_low=lo,
        ci_high=hi,
        detail=f"cached run from {payload.get('measured_on', 'an undated run')} "
               f"at commit {payload.get('commit', 'unknown')}",
    )


def certify(criterion: AcceptanceCriterion = CRITERION,
            lomo_cache: Path = LOMO_CACHE) -> EvidenceCertificate:
    """Measure both scorers and return what the device is entitled to claim."""
    scorers = {}
    loo = _loo_evidence(criterion)
    scorers[loo.name] = loo
    lomo = _lomo_evidence(criterion, lomo_cache)
    scorers[lomo.name] = lomo
    return EvidenceCertificate(scorers=scorers, criterion=criterion)


def main() -> int:
    cert = certify()
    print("EVIDENCE CERTIFICATE — may the device say PASS?\n")
    print(cert.report())
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Simulation")
    print("  numbers are proxies from toy models. Nothing here is evidence about")
    print("  multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
