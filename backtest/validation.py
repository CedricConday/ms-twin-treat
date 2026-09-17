"""Input guards for the benchmark's scoring path.

These check SOFTWARE inputs only. They make malformed arrays fail loudly at the
point they enter the scorer instead of travelling on as a finite-looking number.
They say nothing about whether a model is biologically right — see
`docs/QUALITY.md` for that distinction.

The wording of the shape error is deliberately the same as the one
`score_delta` has always raised, so existing callers that match on it keep
working.
"""

from __future__ import annotations

import numpy as np


def validate_profile_pair(control: np.ndarray, perturbed: np.ndarray,
                          *, name: str = "profile") -> None:
    """Raise a useful error when a pair of expression/delta vectors is unusable.

    A NaN here is the failure mode that matters: it survives arithmetic, lands
    in a Pearson correlation as a quiet 0.0 or nan, and is then averaged into an
    aggregate that still looks like a score.
    """
    control = np.asarray(control)
    perturbed = np.asarray(perturbed)
    if control.ndim != 1 or perturbed.ndim != 1:
        raise ValueError(f"{name} arrays must be one-dimensional, "
                         f"got {control.ndim}d and {perturbed.ndim}d")
    if control.shape != perturbed.shape:
        raise ValueError(f"shape mismatch: {name} {control.shape} vs {perturbed.shape}")
    if not np.issubdtype(control.dtype, np.number) or not np.issubdtype(perturbed.dtype, np.number):
        raise TypeError(f"{name} arrays must be numeric, "
                        f"got {control.dtype} and {perturbed.dtype}")
    if not np.isfinite(control).all() or not np.isfinite(perturbed).all():
        raise ValueError(f"{name} arrays must contain only finite values "
                         f"(found NaN or inf)")


def validate_top_k(top_k: int) -> None:
    """Validate the ranking metric's cutoff.

    `top_k=0` used to yield a precision of 0.0 and `top_k=-1` a nonsensical
    slice, both without complaint. Neither is a meaningful request.
    """
    if isinstance(top_k, bool) or not isinstance(top_k, (int, np.integer)) or top_k < 1:
        raise ValueError(f"top_k must be a positive integer, got {top_k!r}")
