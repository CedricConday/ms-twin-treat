"""Docs must agree with what the code prints. Cheap claims only.

WHY
---
On 2026-09-20 three separate numbers in this repo were stale, and every one was
a measurement written down by hand and never re-run: a hardcoded "9/14, 1/9"
PRINTED by backtest/clinical_velez.py, the same pair on README's front page, and
BUILD_PLAN §8.0's four-arm table. `scripts/state_of_build.py` makes regenerating
them one command. This makes *disagreeing* with them a test failure.

SCOPE, DELIBERATELY NARROW
---------------------------
Only claims that can be recomputed in milliseconds are checked here — the dial
ceiling, and the screen counts, which are read from the artifact the screen
itself wrote. The clinical gates and the LOMO scorers take minutes and belong in
`scripts/state_of_build.py`, not in the test suite: a slow gate that must pass
is a gate people learn to make pass.

If one of these fails, the doc is wrong, not the test. Regenerate and commit.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _read(rel: str) -> str:
    p = ROOT / rel
    if not p.exists():
        pytest.skip(f"{rel} not present")
    # Emphasis is the doc's business, not the test's: a row the author chose to
    # bold is the same claim as one they did not.
    return p.read_text().replace("**", "").replace("*", "")


def test_the_port_scope_quotes_the_live_dial_ceiling():
    """PERNICE_PORT_SCOPE.md's table is the basis for a port/no-port decision."""
    from scripts.dial_ceiling import dial_ceiling

    r = dial_ceiling()
    doc = _read("docs/PERNICE_PORT_SCOPE.md")
    for key in ("in_sample", "singletons_dropped", "out_of_sample"):
        v = r[key]
        row = f"{v['mae']:.1f}pp | {v['null']:.1f}pp | {v['headroom']:.1f}pp"
        assert row in doc, (
            f"docs/PERNICE_PORT_SCOPE.md does not carry the live {key} row "
            f"({row}); re-run scripts/dial_ceiling.py and update the doc"
        )


def test_the_screen_doc_matches_the_screen_artifact():
    """SCREEN_RESULTS.md is generated, so a mismatch means it was hand-edited."""
    raw = _read("results/screen.json")
    counts = json.loads(raw)["counts"]
    doc = _read("docs/SCREEN_RESULTS.md")
    for verdict, n in counts.items():
        assert re.search(rf"\|\s*{verdict}\s*\|\s*{n}\s*\|", doc), (
            f"{verdict}={n} in results/screen.json is not in docs/SCREEN_RESULTS.md; "
            "regenerate with `python -m screen.report`"
        )


def test_the_capacity_screen_doc_names_its_parameter():
    """A capacity result unlabelled by K is unreproducible — the two disagree."""
    doc = _read("docs/SCREEN_RESULTS_K2000.md")
    assert "K = 2000" in doc or "K=2000" in doc
    assert "DIFFERENT MODEL" in doc.upper()


def test_the_harm_channel_statistic_is_quoted_at_its_real_size():
    """0.048 is a pre-registered enlargement, not the miscount that preceded it."""
    from bricks.harm_channel import ranking, ranking_enlarged, top_k_probability

    assert len(ranking()) == 6
    assert len(ranking_enlarged()) == 7
    p = top_k_probability(len(ranking_enlarged()), 2)
    doc = _read("docs/HARM_CHANNEL_PREREG.md")
    assert f"{p:.3f}" in doc


def test_readme_points_at_the_regenerated_state_file():
    """The front page must defer to results/STATE.md, not carry rival numbers."""
    readme = _read("README.md")
    assert "results/STATE.md" in readme
    assert "state_of_build" in readme
