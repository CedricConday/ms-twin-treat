"""Every figure the gate's docs quote is either LIVE or PINNED. Anything else is drift.

Tonight a cached mechanism-holdout number measured on a 12-arm exam survived two
arm-set changes into a 15-arm repo, and the reader-facing page quoted 12-arm
figures throughout after six arms were wired underneath it. Both were found by
reading. Reading does not scale and does not run in CI.

So: this script recomputes every quantity the gate lane reports, scans
`gate/*.py` and the gate docs for figures written as `N.Npp`, and classifies each
one.

    LIVE    matches a currently measured quantity, within TOLERANCE
    PINNED  listed in PINNED_FIGURES with a reason it must not track the live
            set -- a historical comparison, a worked example, a bound
    DRIFT    neither. Either the measurement moved or the figure was wrong.

A DRIFT is not automatically a bug: adding an arm moves numbers, and that is the
point of adding arms. It means a human has to look. Exit code 1 so CI can make
that non-optional.

WHAT THIS DOES NOT CHECK, STATED SO IT IS NOT MISTAKEN FOR COVERAGE
--------------------------------------------------------------------
It checks that quoted figures correspond to measurements. It cannot check that a
figure is quoted about the right thing -- "21.4pp" attached to the wrong scorer
passes here. Prose stays a human's job; this catches the arithmetic going stale
underneath it.

Run:  PYTHONPATH=. python3 scripts/verify_gate_docs.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOLERANCE = 0.051          # figures are quoted to one decimal
FIGURE = re.compile(r"(\d+\.\d)pp")

SCANNED = [
    "gate/__init__.py", "gate/ceiling.py", "gate/headroom.py", "gate/device.py",
    "gate/evidence.py", "gate/criterion.py", "gate/provenance.py",
    "docs/DECISION_GATE.md", "docs/RECOVERABILITY.md",
]

# Figures that must NOT track the live measurement, each with the reason.
# A figure here is a claim that something was true at a stated time, not a claim
# about the current exam.
PINNED_FIGURES = {
    "0.8": "the out-of-sample prize BEFORE six arms were wired, retained because the "
           "page reports that widening the exam moved it to 0.9pp.",
}


def live_measurements() -> dict[str, float]:
    """Recompute everything the gate lane reports. No cached values except the
    LOMO run, which carries its own arm-set fingerprint check."""
    from backtest.loo import run_loo
    from gate.ceiling import recoverability, run_ceiling
    from gate.evidence import certify
    from gate.headroom import dial_ceiling, run_headroom

    out: dict[str, float] = {}

    loo = run_loo()
    out["loo.mae"] = loo["mae"]
    out["loo.null"] = loo["null_mae"]

    cert = certify()
    lomo = cert.scorers.get("leave-one-mechanism-out")
    if lomo is not None and lomo.status == "MEASURED":
        out["lomo.mae"] = lomo.mae
        out["lomo.null"] = lomo.null_mae

    ceil = run_ceiling()
    for key in ("mae_reachable", "null_mae_reachable", "mae_unreachable",
                "null_mae_unreachable", "mae_all", "null_mae_all"):
        out[f"ceiling.{key}"] = ceil[key]
    for fit in ceil["fits"]:
        out[f"ceiling.{fit.arm}.error"] = fit.error
        out[f"ceiling.{fit.arm}.null"] = fit.null_error

    rec = recoverability()
    out["recoverability.mae"] = rec["mae"]
    out["recoverability.null"] = rec["null_mae"]
    # Per-arm rows too: docs/RECOVERABILITY.md tabulates them, and a checker that
    # only knows aggregates reports every correct per-arm figure as drift, which
    # is the fastest way to teach someone to ignore it.
    for row in rec["rows"]:
        out[f"recoverability.{row['arm']}.error"] = row["error"]
        out[f"recoverability.{row['arm']}.null"] = row["null_error"]

    # The capacity comparison lives in backtest/lomo_capacity.py's artifact rather
    # than in gate/, but gate/provenance.py and DECISION_GATE.md both quote it to
    # justify why a verdict must name its model. A figure this lane quotes is a
    # figure this lane has to keep current.
    cap = ROOT / "results" / "lomo_capacity.json"
    if cap.exists():
        import json as _json

        payload = _json.loads(cap.read_text())
        out["capacity.fold_fitted.mae"] = payload["mae"]
        out["capacity.null"] = payload["null_mae"]
        for k, row in (payload.get("fixed_k_diagnostic") or {}).items():
            label = "transcription" if k == "None" else f"K{float(k):g}"
            out[f"capacity.{label}.mae"] = row["mae"]
        base = out.get("capacity.transcription.mae")
        k2000 = out.get("capacity.K2000.mae")
        if base is not None and k2000 is not None:
            out["capacity.gap_transcription_vs_K2000"] = abs(base - k2000)

    head = run_headroom()
    out["headroom.mae"] = head["mae"]
    out["headroom.null"] = head["null_mae"]
    out["headroom.excl_singletons"] = head["mae_excl_singletons"]
    out["headroom.excl_singletons_null"] = head["null_mae_excl_singletons"]

    dial = dial_ceiling()
    for row in ("in_sample", "multi_only", "out_of_sample"):
        out[f"dial.{row}.mae"] = dial[row]["mae"]
        out[f"dial.{row}.null"] = dial[row]["null_mae"]
        out[f"dial.{row}.headroom"] = dial[row]["headroom"]

    return out


def main() -> int:
    live = live_measurements()

    print(f"LIVE MEASUREMENTS — {len(live)} quantities\n")
    for name, value in sorted(live.items()):
        if "." in name and name.count(".") == 1:      # aggregates only, per-arm is noise here
            print(f"  {name:<34} {value:>7.2f}pp")

    drift: list[tuple[str, int, str]] = []
    pinned_seen: set[str] = set()
    live_count = 0

    for rel in SCANNED:
        path = ROOT / rel
        if not path.exists():
            continue
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            for match in FIGURE.finditer(line):
                figure = match.group(1)
                if figure in PINNED_FIGURES:
                    pinned_seen.add(figure)
                    continue
                if any(abs(float(figure) - v) < TOLERANCE for v in live.values()):
                    live_count += 1
                    continue
                drift.append((rel, lineno, line.strip()[:96]))

    print(f"\n  {live_count} figure(s) match a live measurement.")
    print(f"  {len(pinned_seen)} pinned figure(s) seen, of {len(PINNED_FIGURES)} declared.")

    stale_pins = set(PINNED_FIGURES) - pinned_seen
    if stale_pins:
        print(f"  {len(stale_pins)} pinned figure(s) no longer appear anywhere and can be "
              f"dropped from PINNED_FIGURES: {', '.join(sorted(stale_pins))}")

    if not drift:
        print("\nNo drift. Every quoted figure is a current measurement or a declared pin.")
        return 0

    print(f"\nDRIFT — {len(drift)} figure(s) match neither:\n")
    for rel, lineno, text in drift:
        print(f"  {rel}:{lineno}")
        print(f"    {text}")
    print("\n  Adding an arm moves numbers, so drift is expected after wiring and is not")
    print("  automatically a bug. Re-measure, then either update the figure or pin it")
    print("  with the reason it must not track the live set.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
