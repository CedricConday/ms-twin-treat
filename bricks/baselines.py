"""Null models for the perturbation backtest — the bar every real brick must clear.

These are not meant to be good. They exist so that "good" is defined by
*beating them*, not by looking plausible. Three nulls, in increasing difficulty:

  IdentityNull        predicts nothing changes (perturbed == control).
                      The floor. A model that cannot beat this has learned
                      literally nothing.

  LeaveOneOutMeanShiftNull
                      predicts the average delta of the OTHER cell types.
                      **THE CANONICAL BAR.** See below.

  GlobalMeanShiftNull predicts the SAME shift for every cell type, averaged over
                      ALL of them — including the one being predicted.


WHICH BAR IS CANONICAL, AND WHY IT MATTERS
------------------------------------------
These two mean-shift nulls are not interchangeable, and grading different bricks
against different ones would quietly make their scores incomparable.

`GlobalMeanShiftNull` is normally constructed from
`PerturbationBenchmark.global_mean_delta()`, which averages the true delta over
every cell type **including the held-out one**. That null therefore sees data a
leave-one-cell-type-out model never gets. It is leaky *in the null's favour*.

`LeaveOneOutMeanShiftNull` is the apples-to-apples null: for cell type X it
averages only the training cell types' deltas, exactly the information a LOCTO
model is allowed. This is the **canonical bar** — the one a new brick must clear
to claim it learned anything about cell-type-specific response.

Report both. The leaky global null is retained as a *harder secondary* check
precisely because of its information advantage: clearing it too is a stronger
statement than clearing the fair bar alone. On Kang IFN-beta the numbers are

    identity                            0.0000
    leave-one-out mean shift (CANON)    0.8166   (0.8232 with the >=0 floor)
    global mean shift (leaky, harder)   0.8498
    bricks/cell_scgpt.py                0.8732   clears both

A brick that clears the canonical bar but not the leaky one has still learned
something real; say so plainly rather than picking whichever bar flatters it.
"""

from __future__ import annotations

import numpy as np


class IdentityNull:
    name = "null:identity(no-change)"

    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        return np.asarray(control_mean, dtype=float).copy()


class GlobalMeanShiftNull:
    """Applies one shared delta to every cell type.

    Construct with the average delta across cell types (e.g. from
    `PerturbationBenchmark.global_mean_delta()`). This encodes the null
    hypothesis "the perturbation is uniform across cell types" — the reference a
    cell-type-aware model must outperform to justify its complexity.
    """

    name = "null:global-mean-shift"

    def __init__(self, shared_delta: np.ndarray) -> None:
        self.shared_delta = np.asarray(shared_delta, dtype=float).ravel()

    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        control_mean = np.asarray(control_mean, dtype=float).ravel()
        return control_mean + self.shared_delta


class LeaveOneOutMeanShiftNull:
    """THE CANONICAL BAR: predict cell type X with the mean delta of the others.

    Same hypothesis as GlobalMeanShiftNull — "the perturbation is uniform across
    cell types" — but restricted to the information a leave-one-cell-type-out
    model actually has. Nothing about X's own perturbed data enters its own
    prediction, so this is the fair comparison, and it is the bar a new brick
    should be graded against first.

    Args:
        benchmark: the PerturbationBenchmark being scored.
        floor:     clamp predicted expression at >= 0. Off by default — a null
                   should stay assumption-free, and the floor is a modelling
                   choice (it is worth ~+0.007 here, 0.8166 -> 0.8232). Turn it
                   on when you want the strictest possible version of the bar.
    """

    def __init__(self, benchmark, floor: bool = False) -> None:
        self.floor = floor
        cts = benchmark.cell_types
        self._delta: dict[str, np.ndarray] = {}
        for held in cts:
            train = [ct for ct in cts if ct != held]
            self._delta[held] = np.mean(
                [benchmark.profiles[ct].true_delta for ct in train], axis=0)

    name = "null:loo-mean-shift(canonical)"

    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        control_mean = np.asarray(control_mean, dtype=float).ravel()
        delta = self._delta[cell_type]
        if self.floor:
            delta = np.maximum(delta, -control_mean)
        return control_mean + delta
