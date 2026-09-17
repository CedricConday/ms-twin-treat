"""Derive SUPPRESSIVE_STRENGTH from the Kang 2018 IFN-beta effect magnitude.

This is the arithmetic behind the one number in `bricks/grounding.py` that is no
longer set by reasoning. Run it to reproduce that number from the data:

    PYTHONPATH=. python scripts/derive_suppressive_strength.py

THE RULE (fixed before the number was looked at, so it cannot be a fit):

  numerator    mean over cell types of  || mean_IFNbeta - mean_control ||_2
               i.e. how far IFN-beta moves an immune cell's expression profile.

  denominator  mean over cell-type PAIRS of  || mean_control_i - mean_control_j ||_2
               i.e. how far apart two DIFFERENT immune cell identities are in the
               same space. This is the in-data yardstick for "a wholesale change
               of immune cell state", which is what treat=1.0 is supposed to mean.

  strength     numerator / denominator, rounded to 2 dp.

Both terms are Euclidean distances between per-cell-type mean expression vectors
over all 15,706 genes, on the log-normalized matrix the rest of the repo scores
(data/kang.py). Cell types whose measured delta fails the harness reliability bar
(PerturbationBenchmark.LOW_RELIABILITY = 0.5) are dropped from BOTH terms — the
same rule backtest/harness.py already uses, for the same reason: a delta that is
mostly sampling noise cannot ground anything. On Kang that drops Megakaryocytes
(63/69 cells, reliability 0.045) and keeps 7 cell types.

INDEPENDENCE — the property that makes this worth doing: nothing above touches a
relapse rate, a trial arm, or any clinical endpoint. Kang is an in-vitro PBMC
stimulation experiment. The number is fixed by the transcriptome and the cell-type
geometry alone, so the clinical gate stays an out-of-sample test of it.

WHAT IT IS NOT: this measures the SIZE of IFN-beta's effect on immune cells, and
the model then reads that size as "fraction of the autoreactive attack suppressed".
That equation is a modelling assumption, not a measurement — a large transcriptional
footprint (IFN-beta induces a big ISG program) is not the same thing as a large
clinical suppression. See GROUNDING.md for what that assumption costs.
"""

from __future__ import annotations

import itertools

import numpy as np

from data.kang import to_benchmark


def derive(bench=None) -> dict:
    """Recompute the strength from the Kang benchmark. Returns the arithmetic."""
    bench = to_benchmark() if bench is None else bench
    noisy = set(bench.unreliable_cell_types)
    keep = sorted(ct for ct in bench.profiles if ct not in noisy)

    deltas = {ct: float(np.linalg.norm(bench.profiles[ct].true_delta)) for ct in keep}
    identity = {
        (a, b): float(np.linalg.norm(bench.profiles[a].control_mean
                                     - bench.profiles[b].control_mean))
        for a, b in itertools.combinations(keep, 2)
    }
    num = float(np.mean(list(deltas.values())))
    den = float(np.mean(list(identity.values())))
    return {
        "cell_types": keep,
        "dropped_unreliable": sorted(noisy),
        "n_genes": len(bench.gene_names),
        "per_cell_type_delta_norm": deltas,
        "mean_delta_norm": num,
        "mean_identity_distance": den,
        "ratio": num / den,
        "strength": round(num / den, 2),
    }


def main() -> int:
    from bricks.grounding import KANG_IDENTITY_DISTANCE, KANG_IFNB_DELTA_NORM, SUPPRESSIVE_STRENGTH

    d = derive()
    print("SUPPRESSIVE_STRENGTH from the Kang 2018 IFN-beta magnitude\n")
    print(f"  genes: {d['n_genes']}   cell types kept: {len(d['cell_types'])}"
          f"   dropped (unreliable): {d['dropped_unreliable']}\n")
    for ct, n in sorted(d["per_cell_type_delta_norm"].items()):
        print(f"    {ct:<20} ||IFNb - ctrl|| = {n:8.3f}")
    print(f"\n  numerator    mean ||IFNb - ctrl||           = {d['mean_delta_norm']:.3f}")
    print(f"  denominator  mean ||ctrl_i - ctrl_j|| (pairs) = {d['mean_identity_distance']:.3f}")
    print(f"  ratio                                        = {d['ratio']:.4f}")
    print(f"  SUPPRESSIVE_STRENGTH = round(ratio, 2)       = {d['strength']:.2f}")

    print("\n  recorded in bricks/grounding.py: "
          f"{KANG_IFNB_DELTA_NORM:.3f} / {KANG_IDENTITY_DISTANCE:.3f} "
          f"-> {SUPPRESSIVE_STRENGTH:.2f}")
    drift = (abs(d["mean_delta_norm"] - KANG_IFNB_DELTA_NORM) > 0.01
             or abs(d["mean_identity_distance"] - KANG_IDENTITY_DISTANCE) > 0.01
             or d["strength"] != SUPPRESSIVE_STRENGTH)
    print("  DRIFT — the module no longer matches the data." if drift
          else "  match: the constant in the module is the one the data gives.")
    return 1 if drift else 0


if __name__ == "__main__":
    raise SystemExit(main())
