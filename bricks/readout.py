"""B8 — clinical readout brick (v0). Maps sim -> clinical proxies.

Reads state["abm_damage"] (myelin damage trajectory) and state["cns_exposure"]
(barrier delivery), writes state["readout"].

THE MICRO->CLINICAL MAP IS OPEN RESEARCH. Mapping a simulated tissue-damage
trajectory to a lesion count or a relapse rate is unsolved; the numbers here are
invented scales, NOT clinical endpoints. This is the single easiest place in the
whole pipeline to accidentally claim a clinical result, so everything is flagged
validated=False and named a "proxy".

**Partly superseded for ARM-LEVEL scoring (blocker 2).** `MAX_RELAPSE` is still
an invented per-patient scale, but the per-patient number is not what the gate
needs: `bricks/sormani.py` supplies a PUBLISHED trial-level map from a lesion
rate ratio to a relapse rate ratio (Sormani & Bruzzi 2013, 31 trials, 18,901
patients, slope 0.52). Use `relapse_ratio_from_arms()` below to compare two arms.
`MAX_RELAPSE` cancels in any such ratio (BUILD_PLAN §8 blocker 2), so it was
never the thing to fix — the fix was to stop scoring an invented absolute and
start scoring a ratio through a map somebody fitted to real trials.

Barrier coupling (the reason B6 exists): the ABM applied the intervention at full
strength. In reality a therapy only helps to the extent it reaches its site of
action. So the benefit the ABM assumed -- proportional to the intervention's
`treat` -- is clawed back by the undelivered fraction (1 - effective). A
peripheral drug (effective=1.0) keeps all of it; a CNS-required drug that barely
crosses (effective~0.08) loses almost all of it. Ignore `effective` and B6 is
decoration.
"""

from __future__ import annotations

import numpy as np

from bricks.sormani import RelapsePrediction, lesion_ratio, predict_relapse_ratio

# Invented mapping scales -- NOT clinical calibration.
MAX_LESIONS = 40      # damage=1.0 -> 40 "lesions" (illustrative ceiling)
MAX_RELAPSE = 1.5     # damage=1.0 -> annualized relapse rate 1.5 (illustrative)


def relapse_ratio_from_arms(treated_lesions, comparator_lesions) -> RelapsePrediction:
    """Arm-level predicted relapse ratio, via the published Sormani map.

    `treated_lesions` and `comparator_lesions` are per-patient `lesion_proxy`
    values for the two arms — arrays or scalars. They are averaged here, because
    Sormani is a TRIAL-level relation between two arms and has no per-patient
    reading; applying it inside a single patient's readout would be using a
    between-trial regression as a within-person one.

    Returns the prediction object, not a float, so that its blind spot and its
    assumed intercept travel with the number. See bricks/sormani.py.
    """
    treated = float(np.mean(np.asarray(treated_lesions, dtype=float)))
    comparator = float(np.mean(np.asarray(comparator_lesions, dtype=float)))
    return predict_relapse_ratio(lesion_ratio(treated, comparator))


class ReadoutStage:
    name = "B8 clinical readout"
    requires: tuple[str, ...] = ("abm_damage", "cns_exposure")

    def run(self, state: dict) -> dict:
        damage = np.asarray(state["abm_damage"], dtype=float).ravel()
        final_damage = float(damage[-1]) if damage.size else 0.0

        cns = state.get("cns_exposure", {}) or {}
        effective = float(cns.get("effective", 1.0))
        interv = state.get("intervention", {}) or {}
        treat = float(interv.get("treat", 0.0))

        # Delivery gating: benefit realized only to the extent the drug arrives.
        undelivered_benefit = (1.0 - effective) * treat
        adjusted_damage = float(np.clip(final_damage + undelivered_benefit, 0.0, 1.0))

        # relapse_proxy is SCORED, not just displayed: backtest/clinical.py averages
        # it across the cohort and takes a ratio against the untreated arm. Rounding
        # it here quantized that ratio -- at treated-arm damage ~0.02 the per-patient
        # value collapsed onto {0.02, 0.03, 0.04, 0.05}, ~25% discretization error on
        # the gate's only quantitative number. Keep the exact value; round for display.
        state["readout"] = {
            "lesion_proxy": round(adjusted_damage * MAX_LESIONS),
            "relapse_proxy": adjusted_damage * MAX_RELAPSE,
            "sim_final_damage": round(final_damage, 4),
            "delivery_adjusted_damage": round(adjusted_damage, 4),
            "effective_exposure": round(effective, 4),
            "treat": treat,
        }
        state["readout_meta"] = {
            "validated": False,
            "note": "micro->clinical map is OPEN RESEARCH; lesion_proxy/relapse_proxy "
                    "are invented scales, not clinical endpoints. Barrier-gated by "
                    "cns_exposure['effective'].",
        }
        return state

    __call__ = run


if __name__ == "__main__":
    import numpy as _np
    peripheral = {"abm_damage": _np.linspace(0, 0.3, 20),
                  "cns_exposure": {"effective": 1.0}, "intervention": {"treat": 0.5}}
    cns_locked = {"abm_damage": _np.linspace(0, 0.3, 20),
                  "cns_exposure": {"effective": 0.08}, "intervention": {"treat": 0.5}}
    r_p = ReadoutStage().run(dict(peripheral))["readout"]
    r_c = ReadoutStage().run(dict(cns_locked))["readout"]
    print("B8 readout — same sim damage (0.30), same drug (treat=0.5), different delivery:")
    print(f"  peripheral (effective=1.00): lesion_proxy={r_p['lesion_proxy']} "
          f"relapse_proxy={r_p['relapse_proxy']:.2f}  adjusted={r_p['delivery_adjusted_damage']}")
    print(f"  CNS-locked (effective=0.08): lesion_proxy={r_c['lesion_proxy']} "
          f"relapse_proxy={r_c['relapse_proxy']:.2f}  adjusted={r_c['delivery_adjusted_damage']}")
    print("  -> barrier matters: the CNS-required drug that can't cross looks worse. "
          "(proxies invented, validated=False)")
