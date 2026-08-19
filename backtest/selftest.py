"""Validate the ruler before measuring anything real — the backtest of the backtest.

We build synthetic control/perturbed data with a KNOWN structure:
  - a shared perturbation component (same in every cell type), plus
  - a cell-type-SPECIFIC component on a subset of genes, plus
  - observation noise.

Then we assert the harness orders three predictors correctly:
  oracle (sees truth)  >  global-mean-shift null  >  identity null (~0)

If that ordering does not hold, the harness cannot tell a good prediction from a
trivial one, and every downstream result would be meaningless. This file is the
gate: it must pass before any real dataset is trusted to the harness.

Run:  python -m backtest.selftest        (from the repo root, venv active)
"""

from __future__ import annotations

import numpy as np

from backtest.harness import PerturbationBenchmark
from bricks.baselines import GlobalMeanShiftNull, IdentityNull


class _Oracle:
    """Cheats: returns the true perturbed mean. Only for validating the ruler."""

    name = "oracle(sees-truth)"

    def __init__(self, perturbed_means: dict[str, np.ndarray]) -> None:
        self._truth = perturbed_means

    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        return self._truth[cell_type]


def _make_synthetic(seed: int = 0, n_genes: int = 300, n_types: int = 6, noise: float = 0.05):
    rng = np.random.default_rng(seed)
    gene_names = [f"g{i}" for i in range(n_genes)]
    cell_types = [f"ct{i}" for i in range(n_types)]

    # Shared perturbation: ~40 genes move the same way in every cell type.
    shared_delta = np.zeros(n_genes)
    shared_idx = rng.choice(n_genes, size=40, replace=False)
    shared_delta[shared_idx] = rng.normal(2.0, 0.5, size=40)

    control_means, perturbed_means = {}, {}
    for ct in cell_types:
        base = rng.normal(5.0, 1.0, size=n_genes)
        # Cell-type-specific perturbation: a different ~20 genes per type.
        specific = np.zeros(n_genes)
        spec_idx = rng.choice(n_genes, size=20, replace=False)
        specific[spec_idx] = rng.normal(1.5, 0.4, size=20)

        control_means[ct] = base + rng.normal(0, noise, size=n_genes)
        perturbed_means[ct] = base + shared_delta + specific + rng.normal(0, noise, size=n_genes)

    return gene_names, control_means, perturbed_means


def main() -> int:
    gene_names, control_means, perturbed_means = _make_synthetic()
    bench = PerturbationBenchmark(gene_names, control_means, perturbed_means)

    oracle = _Oracle(perturbed_means)
    identity = IdentityNull()
    mean_shift = GlobalMeanShiftNull(bench.global_mean_delta())

    r_oracle = bench.evaluate(oracle)["delta_pearson"].mean()
    r_shift = bench.evaluate(mean_shift)["delta_pearson"].mean()
    r_identity = bench.evaluate(identity)["delta_pearson"].mean()

    print(f"  oracle   delta_pearson = {r_oracle: .4f}  (expect ~1.0)")
    print(f"  meanshift delta_pearson = {r_shift: .4f}  (expect mid: gets shared, misses specific)")
    print(f"  identity delta_pearson = {r_identity: .4f}  (expect ~0.0)")

    ok = True
    checks = [
        ("oracle ~ 1.0", r_oracle > 0.98),
        ("oracle beats mean-shift", r_oracle > r_shift + 0.05),
        ("mean-shift beats identity", r_shift > r_identity + 0.05),
        ("identity ~ 0.0", abs(r_identity) < 0.05),
    ]
    for label, passed in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}] {label}")
        ok = ok and passed

    # And the report gate: the oracle must be flagged as beating both nulls.
    rep = bench.report(oracle, nulls=[identity, mean_shift])
    print(f"  report.beats_all_nulls (oracle) = {rep['beats_all_nulls']}  (expect True)")
    ok = ok and rep["beats_all_nulls"]

    print("\n  RESULT:", "HARNESS VALIDATED — the ruler measures." if ok
          else "HARNESS BROKEN — do not trust downstream results.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
