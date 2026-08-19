"""Null models for the perturbation backtest — the bar every real brick must clear.

These are not meant to be good. They exist so that "good" is defined by
*beating them*, not by looking plausible. Two nulls, in increasing difficulty:

  IdentityNull        predicts nothing changes (perturbed == control).
                      The floor. A model that cannot beat this has learned
                      literally nothing.

  GlobalMeanShiftNull predicts the SAME shift for every cell type — the average
                      perturbation effect. This is the honest bar for any claim
                      of cell-type specificity: IFN-beta does have a shared
                      signature, so reproducing the average is easy. A model
                      earns the word "specific" only by beating this.
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
