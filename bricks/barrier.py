"""B6 — blood/CSF/CNS barrier brick (v0 TOY three-compartment PBPK).

Answers one question for the rest of the pipeline: of a dose given in blood,
what fraction reaches the CNS? For MS this is the difference between a therapy
that can act on a lesion and one that can only act peripherally.

Three compartments, linear transfer, integrated with scipy:

    dP/dt = -(k_pc + k_el) * P + k_cp * C
    dC/dt =  k_pc * P - (k_cp + k_cn) * C + k_nc * N
    dN/dt =  k_cn * C - k_nc * N

    P plasma   C csf   N cns parenchyma

**What is real:** the integration, the mass balance, the shape of the answer
(AUC in CNS over AUC in plasma), AND — now — the baseline CNS-penetration
MAGNITUDE. `k_pc` is calibrated so a large molecule reaches ~0.15% of blood
level in the CNS at baseline, matching the published therapeutic-antibody figure
of ~0.1-0.2% (Pardridge 2019, Expert Opin Investig Drugs,
doi:10.1080/13543784.2019.1627325 — via PubMed). Active-lesion BBB disruption
raises it (~0.6% at full disruption), which is directionally correct.

**What is NOT real:** the OTHER rate constants (k_cp, k_cn, k_nc, k_el) remain
illustrative, and the disruption scaling is directional, not fitted. Only the
baseline CNS-penetration fraction is grounded to an independent published value
(never to an MS outcome). Output stays `validated: False` and must not be read
as a full pharmacokinetic prediction.

Barrier state matters in MS specifically — the BBB is disrupted in active
lesions — so `bbb_disruption` scales the plasma-to-CSF rate. That coupling is
directionally right and quantitatively invented.

Writes state["cns_exposure"] = {"fraction": float, "validated": False, ...}.
Swap-in later: a real published PBPK model behind this same Stage.
"""

from __future__ import annotations

from typing import Any

import numpy as np

try:
    from scipy.integrate import solve_ivp
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover - scipy is pinned, but never hard-fail a brick
    _HAVE_SCIPY = False


# illustrative rate constants (1/h). NOT fitted to anything.
DEFAULTS = {
    "k_pc": 0.00035, # plasma -> csf — GROUNDED: calibrated so a large molecule reaches
                     # ~0.15% of blood level in CNS at baseline (Pardridge 2019, see docstring)
    "k_cp": 0.120,   # csf -> plasma   (efflux dominates influx)
    "k_cn": 0.050,   # csf -> cns parenchyma
    "k_nc": 0.080,   # cns -> csf
    "k_el": 0.200,   # plasma elimination
}


def simulate(dose: float = 1.0, hours: float = 48.0, bbb_disruption: float = 0.0,
             params: dict[str, float] | None = None, n_points: int = 400) -> dict[str, Any]:
    """Integrate the three-compartment model. Returns trajectories and AUC ratio."""
    p = {**DEFAULTS, **(params or {})}

    # Active MS lesions leak. Scale influx up to 4x at full disruption -
    # directionally correct, magnitude invented.
    disruption = float(np.clip(bbb_disruption, 0.0, 1.0))
    k_pc = p["k_pc"] * (1.0 + 3.0 * disruption)

    def rhs(_t, y):
        P, C, N = y
        dP = -(k_pc + p["k_el"]) * P + p["k_cp"] * C
        dC = k_pc * P - (p["k_cp"] + p["k_cn"]) * C + p["k_nc"] * N
        dN = p["k_cn"] * C - p["k_nc"] * N
        return [dP, dC, dN]

    t_eval = np.linspace(0.0, hours, n_points)
    y0 = [float(dose), 0.0, 0.0]

    if _HAVE_SCIPY:
        sol = solve_ivp(rhs, (0.0, hours), y0, t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-10)
        if not sol.success:
            raise RuntimeError(f"barrier ODE failed to integrate: {sol.message}")
        t, (P, C, N) = sol.t, sol.y
    else:  # explicit fallback so the brick still runs without scipy
        dt = hours / (n_points - 1)
        t = t_eval
        P, C, N = np.zeros(n_points), np.zeros(n_points), np.zeros(n_points)
        P[0] = dose
        for i in range(1, n_points):
            dP, dC, dN = rhs(0, [P[i - 1], C[i - 1], N[i - 1]])
            P[i], C[i], N[i] = P[i - 1] + dP * dt, C[i - 1] + dC * dt, N[i - 1] + dN * dt

    auc_p = float(np.trapezoid(P, t))
    auc_n = float(np.trapezoid(N, t))
    fraction = auc_n / auc_p if auc_p > 0 else 0.0

    return {
        "t": t, "plasma": P, "csf": C, "cns": N,
        "auc_plasma": auc_p, "auc_cns": auc_n,
        "fraction": float(fraction),
        "bbb_disruption": disruption,
    }


class BarrierStage:
    """Stage: compute CNS exposure for whatever intervention is in the state.

    An agent that does not need to reach the CNS (`cns_required: False`) is not
    penalised — its effective exposure is 1.0, because the barrier is not on its
    path. That distinction is the whole reason this brick exists: peripheral
    tolerance induction and a remyelination agent face completely different
    delivery problems, and a model that ignores that will rank them wrongly.
    """

    name = "B6 barrier/PBPK"
    requires: tuple[str, ...] = ("intervention",)

    def __init__(self, hours: float = 48.0, params: dict[str, float] | None = None) -> None:
        self.hours = hours
        self.params = params

    def run(self, state: dict) -> dict:
        interv = state.get("intervention") or {}
        dose = float(interv.get("dose", 1.0) or 0.0)
        cns_required = bool(interv.get("cns_required", False))
        disruption = float(state.get("bbb_disruption", 0.3))

        sim = simulate(dose=max(dose, 0.0), hours=self.hours,
                       bbb_disruption=disruption, params=self.params)

        effective = sim["fraction"] if cns_required else 1.0

        state["cns_exposure"] = {
            "fraction": sim["fraction"],
            "effective": float(effective),
            "cns_required": cns_required,
            "bbb_disruption": disruption,
            "auc_plasma": sim["auc_plasma"],
            "auc_cns": sim["auc_cns"],
            "engine": "scipy.solve_ivp" if _HAVE_SCIPY else "explicit-euler-fallback",
            "validated": False,
            "note": "toy 3-compartment PBPK; rate constants illustrative, not fitted",
        }
        state["cns_traj"] = {"t": sim["t"], "plasma": sim["plasma"],
                            "csf": sim["csf"], "cns": sim["cns"], "validated": False}
        return state

    __call__ = run
