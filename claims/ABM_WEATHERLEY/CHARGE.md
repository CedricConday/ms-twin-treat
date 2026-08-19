# CHARGE — ground the population brick against the real Weatherley MS ABM
**Owner: cbk.** Claim it: `mkdir claims/ABM_WEATHERLEY` is already made — put `owner: cbk` in an OWNER file and go.

## The job
`bricks/abm.py` is a TOY grid (invented aggression/damage rates). Ground it against
the **one real open MS agent-based model**: `Georgia-Weatherley/MS_ABM_Weatherley`
(MIT, PLoS Comput Biol 2026 — T cells × BBB, macrophages, oligodendrocyte stress,
myelin agents, remyelination). Replace the toy rules/rates with rules/rates taken
from that published model, cited. Goal: the population brick's behavior is grounded
in real MS-ABM science, not made up.

## Hard constraints (do NOT break these)
1. **Keep the Stage interface exactly:** `run(state) -> state`, writing
   `state["abm_damage"]` (np.ndarray, damage over time) + `state["abm_meta"]`
   (dict with `validated: False`). The spine and `results/experiment.py` depend on it.
2. **Preserve the harm coupling** (added this session): the sim reads
   `intervention["treat"]` (lowers immune aggression) AND `intervention["immunogenic"]`
   (raises it — the harm mechanism). Effective aggression must still fall with treat
   and rise with immuno, or the clinical gate's APL arm breaks.
3. **Preserve** per-patient `state["seed"]`.
4. **validated=False** until the ABM reproduces a Weatherley result. Grounding the
   rates is not the same as validating the model — say so. Cite the Weatherley source
   in the docstring + `abm_meta`.
5. Do it as a grounded rewrite of `abm.py`, OR a new `bricks/abm_weatherley.py` behind
   the same interface — your call. If new, leave a NOTES file so vps wires it.

## Files
- **You own:** `bricks/abm.py` (+ any new abm module you add), and the ABM row in `GROUNDING.md`.
- **Do NOT touch** (vps's live thread): `bricks/barrier.py`, `bricks/intervention.py`,
  `bricks/grounding.py`, `backtest/clinical.py`. If you need a change there, drop a
  note in `claims/ABM_WEATHERLEY/NOTES`.
- Shared docs (README/RESULTS/GROUNDING): append/update only your brick's own row or
  section; don't rewrite vps's.

## Guardrails (unchanged)
- Author `cedric@condaydigital.com`, never a personal email.
- **No `Claude-Session:` trailer in commit messages** (the pre-push hook blocks it and
  it has to be stripped from history before every push — just leave it out).
- No datasets/weights committed. Selftest + clinical gate green before you commit:
  `python -m backtest.selftest`, `python -m backtest.clinical`.
- Built ≠ validated. Nothing here is evidence about MS.

## Definition of done
`bricks/abm.py` (or the new module) runs with rates grounded in the Weatherley model,
cited; the harm coupling + Stage interface intact; a test asserting treatment lowers
and immunogenic raises damage still passes; `GROUNDING.md` ABM row updated; committed
and (when the tree is clean) pushable.
