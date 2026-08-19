"""B2 — Cell brick (v0): cross-cell-type perturbation transfer.

The question this brick answers: given how *other* immune cell types respond to
IFN-beta, can we predict a held-out cell type's response? That is the real
perturbation-transfer problem, and it is the honest first contender against the
global-mean-shift bar (~0.85 on Kang).

Model (deliberately simple, interpretable, and *different from the null*):
per gene, an affine map   delta_g = a_g + b_g * control_g   fit across the
training cell types by least squares. The mean-shift null is the special case
b_g = 0 (predict the same average delta everywhere). The b_g term lets the
prediction bend with a cell type's own baseline — which is exactly the
cell-type-specific signal the null throws away. If b_g helps, we beat the bar;
if it doesn't, we tie it and say so. No faked win.

Evaluation is leave-one-cell-type-out (LOCTO): each cell type is predicted by a
model fit only on the others. That is a genuine held-out test, not a fit-and-report.

Interfaces:
  - PerturbationModel:  predict(control_mean, cell_type) -> perturbed_mean   (harness)
  - Stage:              __call__(state) -> state  writes state["cell_delta"]  (spine)

Later: swap this affine map for scGPT embeddings behind the same interface.
"""

from __future__ import annotations

import numpy as np

from backtest.harness import PerturbationBenchmark


def _fit_affine(controls: np.ndarray, deltas: np.ndarray, ridge: float = 1.0):
    """Per-gene affine fit delta = a + b*control across training cell types.

    controls, deltas : (n_train_celltypes, n_genes)
    Returns (a, b), each (n_genes,). Ridge-regularizes the slope so that with few
    cell types it degrades gracefully toward the mean-shift null (b->0), never
    blows up. That graceful degradation is the point: worst case we equal the
    null, we never do worse by overfitting 7 points.
    """
    x = controls
    y = deltas
    n = x.shape[0]
    xbar = x.mean(axis=0)
    ybar = y.mean(axis=0)
    xc = x - xbar
    yc = y - ybar
    # per-gene slope via regularized covariance / variance
    cov = (xc * yc).sum(axis=0)
    var = (xc * xc).sum(axis=0) + ridge * n
    b = cov / var
    a = ybar - b * xbar
    return a, b


class CellTransferModel:
    """LOCTO affine cross-cell-type transfer over a PerturbationBenchmark."""

    name = "cell:affine-transfer(LOCTO)"

    def __init__(self, benchmark: PerturbationBenchmark, ridge: float = 1.0) -> None:
        self.benchmark = benchmark
        self.ridge = ridge
        # Precompute a leave-one-out affine fit per cell type.
        self._ab: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        cts = benchmark.cell_types
        for held in cts:
            train = [ct for ct in cts if ct != held]
            controls = np.vstack([benchmark.profiles[ct].control_mean for ct in train])
            deltas = np.vstack([benchmark.profiles[ct].true_delta for ct in train])
            self._ab[held] = _fit_affine(controls, deltas, ridge)

    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        control_mean = np.asarray(control_mean, float).ravel()
        a, b = self._ab[cell_type]
        pred_delta = a + b * control_mean
        return control_mean + pred_delta

    # --- Stage interface (spine) -------------------------------------------
    def __call__(self, state: dict) -> dict:
        cell_delta = {}
        for ct, prof in self.benchmark.profiles.items():
            cell_delta[ct] = self.predict(prof.control_mean, ct) - prof.control_mean
        state["cell_delta"] = cell_delta  # real output, not a STANDIN
        return state


if __name__ == "__main__":
    from data.kang import to_benchmark
    from bricks.baselines import GlobalMeanShiftNull, IdentityNull

    print("B2 cell-transfer — scoring against the bar on real Kang IFN-beta data...")
    bench = to_benchmark()
    model = CellTransferModel(bench)
    ident = IdentityNull()
    mshift = GlobalMeanShiftNull(bench.global_mean_delta())

    r_model = bench.evaluate(model)["delta_pearson"]
    r_shift = bench.evaluate(mshift)["delta_pearson"]
    r_ident = bench.evaluate(ident)["delta_pearson"]

    print("\nper-cell-type delta_pearson:")
    for ct in bench.cell_types:
        tag = "  <-- beats bar" if r_model[ct] > r_shift[ct] else ""
        print(f"  {ct:<20} model={r_model[ct]:.3f}  meanshift={r_shift[ct]:.3f}{tag}")
    print(f"\naggregate: identity={r_ident.mean():.4f}  "
          f"meanshift(bar)={r_shift.mean():.4f}  model={r_model.mean():.4f}")
    verdict = "BEATS the bar" if r_model.mean() > r_shift.mean() else \
              "ties/does not beat the bar (honest: too few cell types to bend much)"
    print(f"verdict: cell-transfer {verdict}")
