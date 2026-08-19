"""Invariant tests for the B2-SCGPT cell brick.

These do not test that the model is biologically right. They test the three
claims the scorecard rests on, which are the ones that would quietly rot:

  1. LOCTO is real — the prediction for a cell type does not move when that cell
     type's *perturbed* data changes. This is the leakage test. If temperature
     selection or the similarity ever started peeking at the held-out fold, this
     is what would catch it.
  2. The similarity is what drives the gain — scrambling it degrades the model.
  3. The non-negativity floor holds: predicted expression never goes below zero
     on log1p data.

Run:  python -m pytest tests/test_cell_scgpt.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

from backtest.harness import PerturbationBenchmark
from bricks.baselines import GlobalMeanShiftNull
from bricks.cell_scgpt import ScGPTCellModel


def _clustered_bench(seed: int = 0, n_genes: int = 300):
    """Toy benchmark where cell types fall into two lineages.

    Cell types within a lineage share a baseline profile *and* a lineage-specific
    perturbation component. That is the structure the model is supposed to
    exploit: baseline similarity predicts response similarity. The mean-shift
    null, which averages across lineages, should not be able to keep up.
    """
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    shared = np.zeros(n_genes)
    shared[rng.choice(n_genes, 40, replace=False)] = rng.normal(2.0, 0.3, 40)

    lineage_base, lineage_resp = {}, {}
    for lin in (0, 1):
        lineage_base[lin] = np.abs(rng.normal(5, 1, n_genes))
        r = np.zeros(n_genes)
        r[rng.choice(n_genes, 40, replace=False)] = rng.normal(2.0, 0.3, 40)
        lineage_resp[lin] = r

    ctrl, pert = {}, {}
    for i in range(8):
        lin = i % 2
        base = np.abs(lineage_base[lin] + rng.normal(0, 0.15, n_genes))
        ctrl[f"ct{i}"] = base
        pert[f"ct{i}"] = np.abs(base + shared + lineage_resp[lin]
                                + rng.normal(0, 0.15, n_genes))
    return PerturbationBenchmark(genes, ctrl, pert)


def _model(bench, **kw):
    # "control" so the test never depends on the scGPT embedding cache existing.
    return ScGPTCellModel(bench, similarity="control", **kw)


def test_predict_shape_and_contract():
    bench = _clustered_bench()
    m = _model(bench)
    for ct, prof in bench.profiles.items():
        out = m.predict(prof.control_mean, ct)
        assert out.shape == prof.control_mean.shape
        assert np.all(np.isfinite(out))


def test_locto_no_leakage_of_held_out_perturbed_data():
    """The prediction for cell type X must not depend on X's perturbed data.

    This is the honesty guarantee in executable form. We corrupt one cell type's
    perturbed profile beyond recognition; every prediction for that cell type
    must be bit-identical, while predictions for the others are free to move.
    """
    bench = _clustered_bench()
    held = bench.cell_types[0]
    baseline = _model(bench).predict(bench.profiles[held].control_mean, held)

    rng = np.random.default_rng(99)
    ctrl = {ct: p.control_mean.copy() for ct, p in bench.profiles.items()}
    pert = {ct: p.perturbed_mean.copy() for ct, p in bench.profiles.items()}
    pert[held] = np.abs(rng.normal(5, 2, len(bench.gene_names)))   # nonsense
    corrupted = PerturbationBenchmark(list(bench.gene_names), ctrl, pert)

    after = _model(corrupted).predict(corrupted.profiles[held].control_mean, held)
    np.testing.assert_allclose(baseline, after, rtol=0, atol=0)


def test_temperature_selection_ignores_the_held_out_fold():
    """Same corruption, but checking the *selected hyperparameter* specifically."""
    bench = _clustered_bench()
    held = bench.cell_types[0]
    before = _model(bench)._chosen_T[held]

    rng = np.random.default_rng(7)
    ctrl = {ct: p.control_mean.copy() for ct, p in bench.profiles.items()}
    pert = {ct: p.perturbed_mean.copy() for ct, p in bench.profiles.items()}
    pert[held] = np.abs(rng.normal(5, 2, len(bench.gene_names)))
    after = _model(PerturbationBenchmark(list(bench.gene_names), ctrl, pert))._chosen_T[held]
    assert before == after


def test_beats_mean_shift_null_when_lineage_structure_exists():
    bench = _clustered_bench()
    model = _model(bench)
    null = GlobalMeanShiftNull(bench.global_mean_delta())
    r_model = bench.evaluate(model)["delta_pearson"].mean()
    r_null = bench.evaluate(null)["delta_pearson"].mean()
    assert r_model > r_null, f"model {r_model:.4f} did not beat null {r_null:.4f}"


def test_scrambled_similarity_destroys_the_gain():
    """If the similarity carries the signal, breaking it must cost us."""
    bench = _clustered_bench()
    good = bench.evaluate(_model(bench))["delta_pearson"].mean()

    m = _model(bench)
    rng = np.random.default_rng(3)
    perm = rng.permutation(len(m.cell_types))
    m._sim = m._sim[np.ix_(perm, perm)]
    for ct in m.cell_types:                      # rebuild with the broken metric
        h = m._idx[ct]
        train = [j for j in range(len(m.cell_types)) if j != h]
        m._pred_delta[ct] = m._combine(train, h, m._chosen_T[ct])
    bad = bench.evaluate(m)["delta_pearson"].mean()
    assert bad < good, f"scrambling the similarity did not hurt ({bad:.4f} vs {good:.4f})"


def test_predicted_expression_respects_the_non_negativity_floor():
    bench = _clustered_bench()
    m = _model(bench, floor=True)
    for ct, prof in bench.profiles.items():
        assert m.predict(prof.control_mean, ct).min() >= -1e-9


def test_stage_interface_writes_cell_delta():
    bench = _clustered_bench()
    state = _model(bench)({})
    assert set(state["cell_delta"]) == set(bench.cell_types)
    for v in state["cell_delta"].values():
        assert v.shape == (len(bench.gene_names),)


def test_missing_scgpt_cache_raises_rather_than_silently_substituting():
    bench = _clustered_bench()
    with pytest.raises(FileNotFoundError):
        ScGPTCellModel(bench, similarity="scgpt", checkpoint="definitely-not-a-checkpoint")
