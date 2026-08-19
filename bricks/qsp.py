"""B4 — QSP / mechanistic ODE brick (v0 TOY model).

A minimal neuroinflammation ODE, integrated with scipy. This is a TOY model, not
a validated MS QSP model (no open MS QSP model exists -- see the research brain).
The *integration* is real; the *model* is illustrative. Every output is flagged
validated=False so nothing downstream mistakes a toy trajectory for biology.

State variables (dimensionless, 0..~1):
  A  activated autoreactive immune cells
  C  pro-inflammatory cytokines
  M  intact myelin (starts ~1.0; damage lowers it)

Dynamics (toy):
  dA/dt = s_A + r_CA * C * A * (1 - A) - mu_A * A - treat * A   immune activation with
                                                        logistic carrying capacity (A<=1),
                                                        cytokine feedback, clearance, drug
  dC/dt = k_C * A - mu_C * C                            cytokines from immune cells, decay
  dM/dt = -k_dmg * A * M + k_rep * (1 - M)              myelin damage by immune cells, slow repair

The (1 - A) term is not cosmetic: without a carrying capacity the cytokine-immune
feedback is unbounded positive feedback and the untreated trajectory diverges to
inf/nan. Real immune responses saturate; the toy model must too.

Intervention coupling: reads state["intervention"]["treat"] if B7 has written it
(0 = untreated, >0 = drug suppressing immune activation). Defaults to untreated.

Writes state["qsp_traj"] = {"t", "A", "C", "M", "validated": False}.
Swap-in later: a real SBML/tellurium MS QSP model behind this same Stage.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

_DEFAULTS = dict(s_A=0.02, r_CA=0.9, mu_A=0.25, k_C=0.8, mu_C=0.5,
                 k_dmg=0.6, k_rep=0.05)


def _rhs(t, y, p, treat):
    A, C, M = y
    dA = p["s_A"] + p["r_CA"] * C * A * (1.0 - A) - p["mu_A"] * A - treat * A
    dC = p["k_C"] * A - p["mu_C"] * C
    dM = -p["k_dmg"] * A * M + p["k_rep"] * (1.0 - M)
    return [dA, dC, dM]


def simulate(params: dict | None = None, treat: float = 0.0,
             t_end: float = 40.0, n: int = 200, y0=(0.15, 0.05, 1.0)) -> dict:
    """Integrate the toy neuroinflammation ODE. Returns trajectories."""
    p = {**_DEFAULTS, **(params or {})}
    sol = solve_ivp(_rhs, (0.0, t_end), list(y0), args=(p, treat),
                    t_eval=np.linspace(0, t_end, n), method="LSODA",
                    rtol=1e-6, atol=1e-9)
    if not sol.success:
        raise RuntimeError(f"QSP integration failed: {sol.message}")
    return {"t": sol.t, "A": sol.y[0], "C": sol.y[1], "M": sol.y[2]}


class QSPBrick:
    """Stage: integrate the toy QSP model, honoring any intervention in state."""

    name = "qsp:toy-neuroinflammation-ode(scipy)"

    def __init__(self, params: dict | None = None, t_end: float = 40.0) -> None:
        self.params = params
        self.t_end = t_end

    def run(self, state: dict) -> dict:
        treat = 0.0
        interv = state.get("intervention")
        if isinstance(interv, dict):
            treat = float(interv.get("treat", 0.0))
        # per-patient params from the VPop (B9) override the defaults
        params = {**(self.params or {}), **(state.get("qsp_params") or {})}
        traj = simulate(params=params, treat=treat, t_end=self.t_end)
        state["qsp_traj"] = {**traj, "treat": treat, "validated": False,
                             "note": "toy neuroinflammation ODE, not a validated MS QSP model"}
        return state

    __call__ = run


if __name__ == "__main__":
    untreated = simulate(treat=0.0)
    treated = simulate(treat=0.6)
    print("B4 QSP toy model (scipy) — final myelin M (higher = healthier):")
    print(f"  untreated: M(final) = {untreated['M'][-1]:.3f}   "
          f"peak cytokine C = {untreated['C'].max():.3f}")
    print(f"  treated:   M(final) = {treated['M'][-1]:.3f}   "
          f"peak cytokine C = {treated['C'].max():.3f}")
    better = treated["M"][-1] > untreated["M"][-1]
    print(f"  sanity: treatment preserves more myelin? {better}  "
          f"(toy model, validated=False)")
