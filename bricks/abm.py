"""B5 — Population / agent-based model brick (v0 TOY).

A minimal spatial ABM: immune agents roam a grid of myelin and degrade it on
contact. Treatment lowers their effective aggression. This is a TOY (validated
=False) standing in for a real PhysiCell + MS_ABM_Weatherley port; the emergent
damage curve is real, the biology is illustrative.

Reads state["intervention"]["treat"] (0..1, fraction of aggression suppressed).
Writes state["abm_damage"] = np.ndarray (mean myelin damage per timestep).
"""

from __future__ import annotations

import numpy as np

try:
    import mesa
    from mesa.space import MultiGrid
    _HAVE_MESA = True
except Exception:  # pragma: no cover - fallback path
    _HAVE_MESA = False


if _HAVE_MESA:

    class _ImmuneAgent(mesa.Agent):
        def __init__(self, model, aggression: float):
            super().__init__(model)
            self.aggression = aggression

        def step(self):
            nbrs = self.model.grid.get_neighborhood(self.pos, moore=True, include_center=False)
            self.model.grid.move_agent(self, self.random.choice(nbrs))
            x, y = self.pos
            if self.model.myelin[x, y] > 0 and self.random.random() < self.aggression:
                self.model.myelin[x, y] = max(0.0, self.model.myelin[x, y] - 0.5)

    class MSABM(mesa.Model):
        def __init__(self, width=20, height=20, n_immune=40, aggression=0.5,
                     treat=0.0, seed=0):
            super().__init__(seed=seed)
            self.grid = MultiGrid(width, height, torus=True)
            self.myelin = np.ones((width, height))
            eff = aggression * (1.0 - float(treat))
            for _ in range(n_immune):
                agent = _ImmuneAgent(self, eff)
                self.grid.place_agent(agent, (self.random.randrange(width),
                                              self.random.randrange(height)))
            self.damage: list[float] = []

        def step(self):
            self.agents.shuffle_do("step")
            self.damage.append(float((1.0 - self.myelin).mean()))


def simulate(n_steps=60, treat=0.0, seed=0, **kw) -> np.ndarray:
    """Run the ABM and return the damage-over-time curve (mean myelin lost)."""
    if not _HAVE_MESA:
        return _simulate_numpy(n_steps=n_steps, treat=treat, seed=seed, **kw)
    model = MSABM(treat=treat, seed=seed, **kw)
    for _ in range(n_steps):
        model.step()
    return np.asarray(model.damage)


def _simulate_numpy(n_steps=60, treat=0.0, seed=0, width=20, height=20,
                    n_immune=40, aggression=0.5) -> np.ndarray:
    """Pure-numpy fallback ABM if mesa is unavailable — same output shape."""
    rng = np.random.default_rng(seed)
    myelin = np.ones((width, height))
    pos = np.column_stack([rng.integers(0, width, n_immune), rng.integers(0, height, n_immune)])
    eff = aggression * (1.0 - float(treat))
    damage = []
    for _ in range(n_steps):
        pos = (pos + rng.integers(-1, 2, pos.shape)) % [width, height]
        for x, y in pos:
            if myelin[x, y] > 0 and rng.random() < eff:
                myelin[x, y] = max(0.0, myelin[x, y] - 0.5)
        damage.append(float((1.0 - myelin).mean()))
    return np.asarray(damage)


class ABMBrick:
    """Stage: run the toy population model, honoring intervention.treat."""

    name = "abm:immune-vs-myelin-grid(mesa)" if _HAVE_MESA else "abm:immune-vs-myelin-grid(numpy)"

    def __init__(self, n_steps: int = 60, seed: int = 0) -> None:
        self.n_steps = n_steps
        self.seed = seed

    def __call__(self, state: dict) -> dict:
        treat = 0.0
        interv = state.get("intervention")
        if isinstance(interv, dict):
            treat = float(interv.get("treat", 0.0))
        damage = simulate(n_steps=self.n_steps, treat=treat, seed=self.seed)
        state["abm_damage"] = damage
        state["abm_meta"] = {"validated": False, "engine": self.name,
                             "final_damage": float(damage[-1]), "treat": treat}
        return state


if __name__ == "__main__":
    untreated = simulate(treat=0.0)
    treated = simulate(treat=0.8)
    print(f"B5 ABM ({'mesa' if _HAVE_MESA else 'numpy'}) — final mean myelin damage (lower=healthier):")
    print(f"  untreated: {untreated[-1]:.3f}")
    print(f"  treated:   {treated[-1]:.3f}")
    print(f"  sanity: treatment reduces damage? {treated[-1] < untreated[-1]}  (toy, validated=False)")
