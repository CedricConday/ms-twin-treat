"""Tests for the clinical readout.

One test, and it exists to stop a wiring change that BUILD_PLAN asked for.
"""

from __future__ import annotations


def test_the_clip_saturates_and_that_is_why_qsp_damage_stays_out():
    """A damage value above 1.0 pegs the proxy, so an unbounded source cannot feed it.

    BUILD_PLAN §8.4 listed "readout must read qsp_damage" as the remaining wiring
    step. Measured 2026-09-20, 78% of untreated Velez runs exceed this clip, so
    wiring it would hide most treatment effects behind the ceiling. This pins the
    mechanism cheaply, without simulating anything.
    """
    import numpy as np

    from bricks.readout import MAX_RELAPSE, ReadoutStage

    stage = ReadoutStage()
    over = stage.run({"abm_damage": np.array([3.0]), "cns_exposure": {"effective": 1.0}})
    at_ceiling = stage.run({"abm_damage": np.array([1.0]), "cns_exposure": {"effective": 1.0}})
    assert over["readout"]["relapse_proxy"] == at_ceiling["readout"]["relapse_proxy"]
    assert over["readout"]["relapse_proxy"] == MAX_RELAPSE
