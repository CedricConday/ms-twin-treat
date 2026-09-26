"""Run every gate and write results/STATE.md. One command, one dated answer.

WHY THIS EXISTS
---------------
On 2026-09-20 three separate stale numbers were found in this repo, and none of
them was stale prose — each was a measurement someone wrote down by hand and
never re-ran:

  * `backtest/clinical_velez.py` PRINTED "the ABM path reports 9/14 direction,
    1/9 magnitude" as a hardcoded string, beside its own freshly computed
    number, after the arm set grew to 17.
  * README.md led with 9/14 and 1/9, dated and wrong.
  * BUILD_PLAN.md §8.0 carried a 4-arm table as the live state, with the
    superseding numbers 400 lines below it.

The common cause is that regenerating the full picture meant running six
commands and reading six outputs, so nobody did it. This makes it one command.

WHAT IT DOES NOT DO
-------------------
It does not interpret. Every line it writes is a number a gate printed plus the
command that printed it, the git SHA it ran at, and the timestamp. No verdicts,
no trends, no "improved since" — those are judgements, they belong in
BUILD_PLAN's dated log where a human signs them.

It is also NOT a CI gate. It runs the slow scorers (the clinical gates are
minutes) and it must never become something that has to pass, because then the
temptation is to make the numbers look good rather than current.

Run:  PYTHONPATH=. python scripts/state_of_build.py
"""

from __future__ import annotations

import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "STATE.md"

# (label, module, regexes to pull out of its stdout). A gate whose output shape
# changes shows up as a MISSING line rather than a silently dropped row.
GATES: list[tuple[str, str, list[tuple[str, str]]]] = [
    ("clinical gate (ABM path)", "backtest.clinical", [
        ("direction", r"DIRECTION gate:\s*(\d+/\d+)"),
        ("magnitude", r"magnitude:\s*(\d+/\d+)"),
    ]),
    ("clinical gate (grounded stack)", "backtest.clinical_velez", [
        ("direction", r"DIRECTION:\s*(\d+/\d+)"),
        ("magnitude", r"magnitude:\s*(\d+/\d+)"),
    ]),
    ("leave-one-ARM-out", "backtest.loo", [
        ("MAE", r"out-of-sample MAE:\s*([\d.]+pp)"),
        ("null", r"predict-the-mean null:\s*([\d.]+pp)"),
    ]),
    ("leave-one-MECHANISM-out", "backtest.lomo", [
        ("MAE", r"out-of-sample MAE:\s*([\d.]+pp)"),
        ("null", r"predict-the-mean null:\s*([\d.]+pp)"),
    ]),
    # Name the variant. Quoting this row without saying which one produced a
    # wrong claim in BUILD_PLAN on 2026-09-21: the headline is K fitted PER FOLD,
    # an extra free parameter fitted inside each fold, while the screen runs at a
    # FIXED K and scores WORSE than the transcription there. Both are reported.
    ("LOMO + capacity (K per fold, then fixed K=2000)", "backtest.lomo_capacity", [
        ("MAE", r"out-of-sample MAE:\s*([\d.]+pp)"),
        ("null", r"predict-the-mean null:\s*([\d.]+pp)"),
        ("fixed K=2000", r"K=2000\.0\s+MAE\s+([\d.]+pp)"),
    ]),
    ("per-drug potency (MRI channel)", "backtest.potency", [
        ("fitted", r"(\d+) arm\(s\) fitted"),
        ("out of range", r"(\d+) arm\(s\) OUT OF RANGE"),
    ]),
    ("dial-level ceiling (out of sample)", "scripts.dial_ceiling", [
        ("MAE", r"out of sample, within group\s+([\d.]+pp)"),
        ("null", r"out of sample, within group\s+[\d.]+pp\s+([\d.]+pp)"),
        ("headroom", r"out of sample, within group\s+[\d.]+pp\s+[\d.]+pp\s+([\d.]+pp)"),
    ]),
    ("harm channel", "bricks.harm_channel", [
        ("enlarged ranking", r"N=(\d+), both harm cases still top-2"),
        ("p", r"top-2:\s*1/\d+ = ([0-9]+\.[0-9]+)"),
    ]),
    # Exam v2 (docs/EXAM_V2_PREREG.md): all 23 arms, interval-scored, per model.
    # The module name carries its arguments; run_gate splits on whitespace.
    ("exam v2, Velez transcription", "backtest.exam_v2 --model velez", [
        ("S1 MAE", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE ([\d.]+pp)"),
        ("null", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE [\d.]+pp vs null ([\d.]+pp)"),
        ("S2 headroom", r"headroom ([\d.]+pp)"),
    ]),
    ("exam v2, Jenner 2026 two-equation probe", "backtest.exam_v2 --model minimal", [
        ("S1 MAE", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE ([\d.]+pp)"),
        ("null", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE [\d.]+pp vs null ([\d.]+pp)"),
        ("S2 headroom", r"headroom ([\d.]+pp)"),
    ]),
    ("exam v2, Pernice 2020 port (Figure S2 readings)", "backtest.exam_v2 --model pernice", [
        ("S1 MAE", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE ([\d.]+pp)"),
        ("null", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE [\d.]+pp vs null ([\d.]+pp)"),
        ("S2 headroom", r"headroom ([\d.]+pp)"),
    ]),
    ("exam v2, Pernice 2020 port (memory-read variant)", "backtest.exam_v2 --model pernice --variant memread", [
        ("S1 MAE", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE ([\d.]+pp)"),
        ("null", r"S1 interval-LOMO, comparator-adjusted \(primary\): MAE [\d.]+pp vs null ([\d.]+pp)"),
        ("S2 headroom", r"headroom ([\d.]+pp)"),
    ]),
    # Exam C (docs/EXAM_COCHRANE_PREREG.md): the Cochrane 2024 network, 10 arms,
    # scored on the cached curves. The row quotes the port's Figure S2 reading,
    # the closest any scorer here has come to its null; the others are in the log.
    ("exam C (Cochrane 2024 network), Pernice port", "backtest.exam_cochrane", [
        ("CS1 MAE", r"\npernice: .*?\n\s*CS1 interval-LOMO\s+MAE ([\d.]+pp)"),
        ("null", r"\npernice: .*?\n\s*CS1 interval-LOMO\s+MAE [\d.]+pp vs null ([\d.]+pp)"),
        ("CS2 headroom", r"\npernice: (?:.*\n)*?\s*CS2 oracle .*?headroom ([\d.]+pp)"),
    ]),
]


def sha() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, timeout=30)
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                               capture_output=True, text=True, timeout=30)
        return out.stdout.strip() + (" (dirty)" if dirty.stdout.strip() else "")
    except Exception:
        return "unknown"


def run_gate(module: str) -> tuple[str, bool]:
    try:
        p = subprocess.run([sys.executable, "-m", *module.split()], cwd=ROOT,
                           capture_output=True, text=True, timeout=3600,
                           env={**__import__("os").environ, "PYTHONPATH": str(ROOT)})
        return p.stdout + p.stderr, p.returncode == 0
    except subprocess.TimeoutExpired:
        return "TIMED OUT after 3600s", False


def main() -> int:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    rows: list[str] = []
    for label, module, patterns in GATES:
        print(f"running {module} ...", flush=True)
        text, ok = run_gate(module)
        if not ok:
            rows.append(f"| {label} | `{module}` | **FAILED TO RUN** |")
            continue
        parts = []
        for name, rx in patterns:
            m = re.search(rx, text)
            parts.append(f"{name} {m.group(1)}" if m else f"{name} **MISSING**")
        rows.append(f"| {label} | `python -m {module}` | {', '.join(parts)} |")

    body = [
        "# STATE OF THE BUILD",
        "",
        f"Generated {stamp} at `{sha()}` by `python scripts/state_of_build.py`.",
        "**Every number here was produced by the command beside it, in this run.**",
        "Nothing on this page is transcribed, and nothing on it is interpreted —",
        "for what the numbers mean, read BUILD_PLAN.md §8.4, where a human signs it.",
        "",
        "| gate | command | result |",
        "|---|---|---|",
        *rows,
        "",
        "## How to read this",
        "",
        "- A gate with a null **loses to its null** wherever MAE > null. That has",
        "  been true of every out-of-sample scorer in this repo since they were built.",
        "- **MISSING** means the gate ran but its output no longer matches the pattern",
        "  this script greps for. That is a script bug, not a result — fix the pattern.",
        "- This page is regenerated, never edited. If a number here disagrees with a",
        "  number written in a doc, the doc is wrong.",
        "",
        "Nothing in this repo is evidence about multiple sclerosis.",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(body) + "\n")
    print(f"\nwritten: {OUT}")
    print("\n".join(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
