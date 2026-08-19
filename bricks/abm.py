"""B5 — Population / agent-based model brick, GROUNDED in the Weatherley MS ABM.

WHAT CHANGED AND WHY
--------------------
This brick used to be an invented toy: immune agents on a 20x20 grid, an
`aggression = 0.5` probability pulled out of the air, myelin decremented by 0.5
on contact. Nothing in it came from MS science.

It is now a port of the one real open agent-based model of multiple sclerosis:

    Georgia-Weatherley/MS_ABM_Weatherley  (MIT licence)
    "Therapeutic targeting of oligodendrocytes in an agent-based model of
     multiple sclerosis", PLOS Computational Biology (open access)
    https://doi.org/10.1371/journal.pcbi.1013273
    https://github.com/Georgia-Weatherley/MS_ABM_Weatherley

Every rule and every rate below is taken from that published MATLAB source, and
each one is annotated with the file it came from. Nothing here is tuned, and
nothing is fitted to any clinical outcome this repo is trying to predict.

THE MODEL (as published)
------------------------
A 2D lattice split into three regions — blood | CNS | myelin — with three agent
populations and an oligodendrocyte layer:

  C1  peripheral T cells, born in blood on a relapse schedule, cross the
      blood-brain barrier eastward with probability `leavingprob`
  C2  resident CNS cells, confined to the CNS strip
  C3  the effector population: created when a C1 and a C2 land on the same
      lattice site (the C1 is consumed). C3 moves with a distance-decaying bias
      toward surviving myelin and demyelinates whatever it sits on.

  myelin   30000 pieces, each with an integer health grade 0..4; a piece touched
           by a C3 drops one grade per timestep
  oligos   1200 oligodendrocytes, each owning a 5x5 block of 25 myelin pieces.
           An oligo that accumulates >= 10 destroyed pieces stops remyelinating;
           at >= 14 it undergoes apoptosis and takes all 25 of its pieces to
           zero. Healthy oligos rebuild one grade per piece every 25 timesteps.

That oligodendrocyte layer is the scientific point of the paper — damage is not
a smooth decay, it is a threshold cascade: lose enough pieces and the whole
5x5 block dies at once. The old toy could not express that at all.

HOW OUR INTERVENTION COUPLES IN — and why this particular lever
---------------------------------------------------------------
The paper models its own therapies, and one of them is exactly the lever we
need: reducing blood-brain-barrier permeability, `PB.leavingprob`, from a
baseline **0.1** to a treated **0.025** (Setup_Simulation.m; ABM_main_script.m
`intervention.new_BBB_prob`). Fewer T cells cross, fewer C3 are made, less
myelin is lost.

So we map our two intervention axes onto the published BBB lever:

    leavingprob = 0.1 * (1 - treat + immunogenic)

  * `treat = 0.75` reproduces the paper's own therapeutic BBB value (0.025)
    exactly. That is the anchor — our scale is pinned to a published number,
    not invented.
  * `immunogenic > 0` pushes permeability the other way, which is the harm
    mechanism the clinical gate's APL arm needs.

Damage therefore still falls with `treat` and rises with `immunogenic`, as the
rest of the stack requires — but now for a mechanistic reason taken from a real
model, instead of by multiplying a made-up aggression constant.

BUILT != VALIDATED — read this before quoting anything
-------------------------------------------------------
`validated=False`, deliberately. Grounding the rates is **not** the same as
reproducing the paper's results. What is true today:

  * the rules and the constants are the published ones, cited per line;
  * the qualitative behaviour is right (threshold cascade, relapse-driven
    damage, treatment reduces damage, immunogenic increases it).

What has **not** been done: no figure from the paper has been reproduced and
checked. Until that happens this is "grounded in", not "validated against", the
Weatherley model, and `abm_meta["validated"]` stays False. See
`claims/ABM_WEATHERLEY/NOTES` for exactly what a validation run would require.

One deliberate deviation, called out rather than buried: the spine runs this
brick once per virtual patient per arm, so the default horizon is shorter than
the paper's 300-day run (`DEFAULT_DAYS` below). The *rates* are per-timestep and
unchanged — only how long we integrate them differs. `simulate(days=300)`
reproduces the published run length.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

# --------------------------------------------------------------------------- #
# THE PUBLISHED PARAMETER SET.
# Every value cited to its file in Georgia-Weatherley/MS_ABM_Weatherley
# (Current_version/). Do not tune these to make an outcome come out right — that
# is what GROUNDING.md forbids. If one is wrong, fix it against the source.
# --------------------------------------------------------------------------- #
W = dict(
    tau_minutes=20,          # Setup_Agents.m   timestep duration
    mins_per_day=1440,       # Setup_Agents.m
    published_days=300,      # ABM_main_script.m  step_cap = 21600 steps

    blood_width=3,           # Setup_Simulation.m
    cns_width=5,             # Setup_Simulation.m
    myelin_width=100,        # Setup_Simulation.m
    domain_height=300,       # Setup_Simulation.m

    oligo_dim=5,             # Setup_Simulation.m  5x5 block per oligodendrocyte
    myelin_grades=4,         # Setup_Simulation.m  health grades 0..4
    heal_time=25,            # Setup_Simulation.m  timesteps between rebuild ticks

    leaving_prob=0.1,        # Setup_Simulation.m  BBB: blood -> CNS
    entering_prob=0.0,       # Setup_Simulation.m  BBB: CNS -> blood
    therapeutic_leaving_prob=0.025,   # ABM_main_script.m  intervention.new_BBB_prob

    c3_aggression=0.5,       # Setup_Simulation.m  bias amplitude
    c3_bias_decay=0.025,     # Setup_Simulation.m  0.5 * exp(-0.025 * distance)
    c3_max_bias_dist=80,     # Biased_Movement.m   beyond this, movement is unbiased

    c1_init=5,               # Setup_Agents.m
    c2_init=45,              # Setup_Agents.m
    death_rate_per_day=0.35,  # Setup_Agents.m  C1.deathrate = 0.35/(1440/20)

    oligo_apop_threshold=14,     # ABM_main_script.m  intervention.oligo_apop
    oligo_stop_my_threshold=10,  # ABM_main_script.m  intervention.oligo_stop_my

    # Renew_Cell_1.m, relapse_schedule = 4 (the schedule ABM_main_script.m uses):
    relapse_days=28,
    relapse_starts=(1, 50, 100, 150, 200, 250),
    prob_relapse=0.0025,
    prob_baseline=0.0002,
)

STEPS_PER_DAY = W["mins_per_day"] // W["tau_minutes"]        # 72
DEATH_RATE = W["death_rate_per_day"] / STEPS_PER_DAY         # per timestep

# Default horizon: ONE COMPLETE PUBLISHED RELAPSE WINDOW. Renew_Cell_1.m's
# schedule 4 (the one ABM_main_script.m runs) uses 28-day relapses, so 28 days is
# a boundary the model itself defines, not a number picked to make an outcome
# come out right. The paper integrates 300 days; the spine runs this brick once
# per virtual patient per arm, which that horizon makes impractical. Per-timestep
# rates are identical either way — only the integration length differs.
DEFAULT_DAYS = 28

CITATION = ("Weatherley et al., PLOS Comput Biol, doi:10.1371/journal.pcbi.1013273; "
            "code github.com/Georgia-Weatherley/MS_ABM_Weatherley (MIT)")

# --------------------------------------------------------------------------- #
# RUN PROFILES — an explicit, MEASURED speed/fidelity trade-off.
#
# The spine runs this brick in a cohort loop (arms x virtual patients), so the
# paper's 300-day run is not tractable there. Every profile therefore keeps the
# PUBLISHED LATTICE and every published rate, threshold and probability. The one
# thing that varies is how long the model is integrated:
#
#   "published"  the paper's own run: 100x300, 300 days. This is the profile a
#                real validation attempt against a paper figure must use.
#   "default"    100x300, 28 days — one complete relapse window under
#                Renew_Cell_1.m schedule 4 (the schedule ABM_main_script.m runs).
#   "gate"       100x300, 14 days — one complete relapse window under
#                Renew_Cell_1.m schedule 3. ~2s/run, so a 48-run cohort gate
#                finishes in ~1.5 min.
#
# Both short horizons are boundaries the published model itself defines, not
# numbers chosen to make an outcome come out right.
#
# MEASURED COST of "gate" vs "default" (3 seeds, relative change in final damage
# — the only quantity the clinical gate reads):
#     treat=0.5    -48.9%  vs  -56.5%
#     immuno=0.4   +26.6%  vs  +23.8%
# Direction is preserved and magnitude is within ~8pp. Absolute damage is NOT
# comparable across horizons and must never be quoted across them.
#
# A shrunken LATTICE was tried and rejected: at 40x100 (agent counts scaled to
# preserve density) the treatment effect collapsed from -56.5% to -23.2% and the
# harm effect inflated from +23.8% to +73.0%, because a narrower myelin field
# saturates. `myelin_width`/`height` remain available as passthroughs, but
# changing them distorts magnitudes — do not use them for the gate.
# --------------------------------------------------------------------------- #
PROFILES = {
    "published": dict(days=W["published_days"], myelin_width=W["myelin_width"],
                      height=W["domain_height"]),
    "default":   dict(days=28, myelin_width=W["myelin_width"], height=W["domain_height"]),
    "gate":      dict(days=14, myelin_width=W["myelin_width"], height=W["domain_height"]),
}


def bbb_permeability(treat: float = 0.0, immuno: float = 0.0) -> float:
    """Map our intervention axes onto the paper's own BBB therapy lever.

    Anchored so that `treat=0.75` gives 0.025 — the exact therapeutic value in
    ABM_main_script.m. Clipped to a probability.
    """
    p = W["leaving_prob"] * (1.0 - float(treat) + float(immuno))
    return float(np.clip(p, 0.0, 1.0))


class _Domain:
    """Lattice geometry and the myelin/oligodendrocyte layer (Setup_Simulation/Setup_Agents)."""

    def __init__(self, myelin_width: int, height: int):
        self.height = height
        self.left_blood = 0
        self.right_blood = self.left_blood + W["blood_width"]        # 3
        self.left_cns = self.right_blood                             # 3
        self.right_cns = self.left_cns + W["cns_width"]              # 8
        self.left_myelin = self.right_cns                            # 8
        self.right_myelin = self.left_myelin + myelin_width - 1
        self.width = W["blood_width"] + W["cns_width"] + myelin_width

        # PB boundary sites (Setup_Simulation.m)
        self.pb_left = self.right_blood          # eastward crossing happens here
        self.pb_right = self.right_blood + 1     # westward crossing happens here

        d = W["oligo_dim"]
        self.n_sub_h = height // d
        self.n_sub_w = myelin_width // d
        self.n_oligos = self.n_sub_h * self.n_sub_w
        self.pieces_per_oligo = d * d

        # Myelin piece coordinates, grouped oligo-block by oligo-block exactly as
        # Setup_Agents.m builds xarray/yarray, so reshape(-1, 25) recovers oligos.
        xs, ys = [], []
        for ii in range(self.n_sub_h):
            for jj in range(self.n_sub_w):
                gx, gy = np.meshgrid(
                    np.arange(self.left_myelin + jj * d, self.left_myelin + (jj + 1) * d),
                    np.arange(height - ii * d, height - (ii + 1) * d, -1),
                )
                xs.append(gx.ravel())
                ys.append(gy.ravel())
        self.piece_x = np.concatenate(xs)
        self.piece_y = np.concatenate(ys)
        self.n_pieces = self.piece_x.size

        # (x, y) -> flat piece index, for O(1) demyelination lookups
        self.index_grid = np.full((self.width + 2, height + 2), -1, dtype=np.int64)
        self.index_grid[self.piece_x, self.piece_y] = np.arange(self.n_pieces)


class WeatherleyABM:
    """The published model, ported to vectorized numpy."""

    def __init__(self, treat=0.0, immuno=0.0, seed=0, myelin_width=None, height=None,
                 bias_refresh_steps=6):
        self.rng = np.random.default_rng(seed)
        self.bias_refresh_steps = max(1, int(bias_refresh_steps))
        self.dom = _Domain(myelin_width or W["myelin_width"], height or W["domain_height"])
        self.leaving_prob = bbb_permeability(treat, immuno)
        self.entering_prob = W["entering_prob"]

        d = self.dom
        g = W["myelin_grades"]
        self.state = np.full((d.n_oligos, d.pieces_per_oligo), g, dtype=np.int16)
        self.timer = np.full((d.n_oligos, d.pieces_per_oligo), -1, dtype=np.int32)
        self.oligo_state = np.ones(d.n_oligos, dtype=np.int8)   # 1 active, 2 stopped, 0 dead

        # Agents (Setup_Agents.m initial placement). If the lattice has been
        # shrunk for speed, initial populations scale with height so that agent
        # DENSITY — which is what drives C1/C2 collisions and therefore C3
        # creation — matches the published model rather than being concentrated.
        hscale = d.height / W["domain_height"]
        n1 = max(1, int(round(W["c1_init"] * hscale)))
        n2 = max(1, int(round(W["c2_init"] * hscale)))
        self.c1x = self.rng.integers(d.left_blood + 1, d.right_blood, n1).astype(np.int64)
        self.c1y = self.rng.integers(1, d.height, n1).astype(np.int64)
        self.c2x = self.rng.integers(d.left_cns + 1, d.right_cns, n2).astype(np.int64)
        self.c2y = self.rng.integers(1, d.height, n2).astype(np.int64)
        self.c3x = np.empty(0, dtype=np.int64)
        self.c3y = np.empty(0, dtype=np.int64)

        self._edt_cache = None
        self._alive_prev = None
        self._alive2d = np.zeros((self.dom.width + 2, self.dom.height + 2), dtype=bool)

    # --- damage readout ----------------------------------------------------
    def damage_fraction(self) -> float:
        """Mean myelin lost, in [0, 1] — the same quantity the old brick reported."""
        return float(1.0 - self.state.mean() / W["myelin_grades"])

    # --- movement (Unbiased_Movement.m) ------------------------------------
    def _move_unbiased(self, x, y, confine_cns=False):
        n = x.size
        if n == 0:
            return x, y
        r = self.rng.random(n)
        west = r < 0.25
        east = (r >= 0.25) & (r < 0.50)
        south = (r >= 0.50) & (r < 0.75)
        north = r >= 0.75

        west &= x != 1                                    # left wall
        blocked = (x == self.dom.pb_right) & (self.rng.random(n) > self.entering_prob)
        west &= ~blocked
        blocked = (x == self.dom.pb_left) & (self.rng.random(n) > self.leaving_prob)
        east &= ~blocked

        x = x + east.astype(np.int64) - west.astype(np.int64)
        y = y + north.astype(np.int64) - south.astype(np.int64)
        if confine_cns:                                   # C2 stays in the CNS strip
            x = np.minimum(x, self.dom.right_cns)
        return x, y

    # --- movement (Biased_Movement.m) --------------------------------------
    # Chebyshev-ball offsets, ordered by true Euclidean distance. Any alive site
    # found within this ball IS the global nearest, because everything outside
    # the ball is strictly farther than its radius.
    _R = 3
    _OFF = None

    @classmethod
    def _offsets(cls):
        if cls._OFF is None:
            r = cls._R
            dx, dy = np.meshgrid(np.arange(-r, r + 1), np.arange(-r, r + 1), indexing="ij")
            dd = np.sqrt(dx ** 2 + dy ** 2).ravel()
            keep = dd <= r
            order = np.argsort(dd[keep])
            cls._OFF = (dx.ravel()[keep][order], dy.ravel()[keep][order], dd[keep][order])
        return cls._OFF

    def _refresh_alive2d(self):
        d = self.dom
        self._alive2d[d.piece_x, d.piece_y] = self.state.ravel() > 0

    def _nearest_alive(self, px, py):
        """EXACT nearest surviving myelin piece for each given point.

        Biased_Movement.m loops every C3 over every living myelin piece and takes
        the minimum. Same answer, two stages:
          1. a local Euclidean ball of radius 3 around each C3 — almost always a
             hit, since C3 sit on or beside the myelin field;
          2. for the few C3 with no myelin that close (late, once large holes
             have opened), one exact distance transform over the alive mask,
             cached until that mask changes.
        No approximation: stage 1 is exact by construction, stage 2 is exact.
        """
        n = px.size
        dx, dy, dd = self._offsets()
        cand_x = px[:, None] + dx[None, :]
        cand_y = py[:, None] + dy[None, :]
        np.clip(cand_x, 0, self._alive2d.shape[0] - 1, out=cand_x)
        np.clip(cand_y, 0, self._alive2d.shape[1] - 1, out=cand_y)
        ok = self._alive2d[cand_x, cand_y]

        first = np.argmax(ok, axis=1)          # offsets are distance-ordered
        found = ok[np.arange(n), first]
        dist = np.where(found, dd[first], np.inf)
        tx = np.where(found, cand_x[np.arange(n), first], px)
        ty = np.where(found, cand_y[np.arange(n), first], py)

        if not found.all():
            miss = ~found
            emap = self._edt()
            if emap is None:                   # no myelin left anywhere
                return dist, tx, ty, False
            edist, eix, eiy = emap
            dist[miss] = edist[px[miss], py[miss]]
            tx[miss] = eix[px[miss], py[miss]]
            ty[miss] = eiy[px[miss], py[miss]]
        return dist, tx, ty, True

    def _edt(self):
        """Exact distance transform over the alive mask; recomputed only on change."""
        alive_flat = self.state.ravel() > 0
        if not alive_flat.any():
            return None
        if (self._alive_prev is not None and self._edt_cache is not None
                and np.array_equal(alive_flat, self._alive_prev)):
            return self._edt_cache
        dist, (ix, iy) = ndimage.distance_transform_edt(~self._alive2d, return_indices=True)
        self._edt_cache = (dist, ix, iy)
        self._alive_prev = alive_flat.copy()
        return self._edt_cache

    def _move_c3(self, step: int = 0):
        n = self.c3x.size
        if n == 0:
            return
        self._refresh_alive2d()
        dist, tx, ty, any_alive = self._nearest_alive(self.c3x, self.c3y)
        if not any_alive:
            self.c3x, self.c3y = self._move_unbiased(self.c3x, self.c3y)
            return
        dx = tx - self.c3x
        dy = ty - self.c3y

        denom = np.abs(dx) + np.abs(dy)
        with np.errstate(invalid="ignore", divide="ignore"):
            xr = np.where(denom > 0, np.abs(dx) / denom, 0.0)
            yr = np.where(denom > 0, np.abs(dy) / denom, 0.0)

        extra = W["c3_aggression"] * np.exp(-W["c3_bias_decay"] * dist)
        close = dist <= W["c3_max_bias_dist"]
        px = (0.25 + xr * extra) * close
        py = (0.25 + yr * extra) * close

        east = np.where(dx >= 0, px, 0.0)
        west = np.where(dx < 0, px, 0.0)
        north = np.where(dy >= 0, py, 0.0)
        south = np.where(dy < 0, py, 0.0)

        # Redistribute whatever probability the bias did not claim (Biased_Movement.m)
        zero_count = ((east == 0).astype(int) + (west == 0).astype(int)
                      + (north == 0).astype(int) + (south == 0).astype(int))
        leftover = 1.0 - (east + west + north + south)
        share = np.divide(leftover, np.maximum(zero_count, 1),
                          out=np.zeros_like(leftover), where=zero_count > 0)
        east = np.where(east == 0, share, east)
        west = np.where(west == 0, share, west)
        north = np.where(north == 0, share, north)
        south = np.where(south == 0, share, south)

        r = self.rng.random(n)
        m_west = r < west
        m_east = (r >= west) & (r < west + east)
        m_south = (r >= west + east) & (r < west + east + south)
        m_north = r >= west + east + south

        blocked = (self.c3x == self.dom.pb_right) & (self.rng.random(n) > self.entering_prob)
        m_west &= ~blocked
        blocked = (self.c3x == self.dom.pb_left) & (self.rng.random(n) > self.leaving_prob)
        m_east &= ~blocked

        self.c3x = self.c3x + m_east.astype(np.int64) - m_west.astype(np.int64)
        self.c3y = self.c3y + m_north.astype(np.int64) - m_south.astype(np.int64)

    # --- myelin dynamics ---------------------------------------------------
    def _remyelinate(self):
        """Oligo_Remyelination.m — healthy oligos rebuild one grade every heal_time."""
        destroyed = (self.state == 0).sum(axis=1)                 # per oligo
        healthy = self.oligo_state == 1
        self.timer[healthy] += 1
        ready = self.timer >= W["heal_time"]
        build = ready & healthy[:, None] & (destroyed < W["oligo_stop_my_threshold"])[:, None]
        self.state[build] += 1
        self.timer[build] = 0
        np.clip(self.state, 0, W["myelin_grades"], out=self.state)
        self.timer[self.state == W["myelin_grades"]] = -1

    def _apoptosis(self):
        """Oligo_Apoptosis.m — the threshold cascade that is the paper's point."""
        destroyed = (self.state == 0).sum(axis=1)
        self.oligo_state[destroyed >= W["oligo_stop_my_threshold"]] = 2
        dead = destroyed >= W["oligo_apop_threshold"]
        self.oligo_state[dead] = 0
        self.state[dead] = 0

    def _demyelinate(self):
        """Demyelination.m — a C3 sitting on a piece knocks it down one grade."""
        if self.c3x.size == 0:
            return
        inside = ((self.c3x >= self.dom.left_myelin) & (self.c3x <= self.dom.right_myelin)
                  & (self.c3y >= 1) & (self.c3y <= self.dom.height))
        if not inside.any():
            return
        idx = self.dom.index_grid[self.c3x[inside], self.c3y[inside]]
        idx = np.unique(idx[idx >= 0])          # 'intersect ... rows' is set-like
        if idx.size == 0:
            return
        flat = self.state.ravel()
        flat[idx] -= 1
        self.timer.ravel()[idx] = -1
        np.clip(self.state, 0, W["myelin_grades"], out=self.state)

    # --- agent bookkeeping -------------------------------------------------
    def _interact(self):
        """C1 + C2 on the same site -> a C3, consuming the C1 (main_function.m)."""
        if self.c1x.size == 0 or self.c2x.size == 0:
            return
        h = self.dom.height + 2
        hit = np.isin(self.c1x * h + self.c1y, self.c2x * h + self.c2y)
        if not hit.any():
            return
        self.c3x = np.concatenate([self.c3x, self.c1x[hit]])
        self.c3y = np.concatenate([self.c3y, self.c1y[hit]])
        self.c1x, self.c1y = self.c1x[~hit], self.c1y[~hit]

    def _renew_c1(self, step: int):
        """Renew_Cell_1.m, relapse schedule 4 — birth rate jumps during relapses."""
        day = step / STEPS_PER_DAY
        in_relapse = any(s <= day <= s + W["relapse_days"] for s in W["relapse_starts"])
        p = W["prob_relapse"] if in_relapse else W["prob_baseline"]
        ys = np.arange(1, self.dom.height)
        for col in (1, 2):                       # two blood columns, as published
            born = ys[self.rng.random(ys.size) < p]
            if born.size:
                self.c1x = np.concatenate([self.c1x, np.full(born.size, col, dtype=np.int64)])
                self.c1y = np.concatenate([self.c1y, born])

    def _boundaries(self):
        """Vertical wrap, reflecting left wall, eviction off the right edge."""
        d = self.dom
        for name in ("c1", "c2", "c3"):
            y = getattr(self, name + "y")
            y = np.where(y == d.height, 1, y)
            y = np.where(y == 0, d.height - 1, y)
            setattr(self, name + "y", y)

        self.c2x = np.minimum(self.c2x, d.right_cns)     # C2 confined to the CNS
        self.c1x = np.where(self.c1x == 0, 1, self.c1x)  # reflecting LHS

        keep = self.c1x != d.width                       # evicted off the far side
        self.c1x, self.c1y = self.c1x[keep], self.c1y[keep]
        keep = self.c3x != d.width
        self.c3x, self.c3y = self.c3x[keep], self.c3y[keep]

    def _deaths(self):
        """C1/C3 apoptosis at the published per-timestep rate (Setup_Agents.m)."""
        keep = self.rng.random(self.c1x.size) >= DEATH_RATE
        self.c1x, self.c1y = self.c1x[keep], self.c1y[keep]
        keep = self.rng.random(self.c3x.size) >= DEATH_RATE
        self.c3x, self.c3y = self.c3x[keep], self.c3y[keep]

    # --- one timestep (order follows main_function.m exactly) --------------
    def step(self, n: int):
        self._move_c3(n)                                             # C3 movement
        self.c1x, self.c1y = self._move_unbiased(self.c1x, self.c1y)
        self.c2x, self.c2y = self._move_unbiased(self.c2x, self.c2y, confine_cns=True)
        self._boundaries()                                           # boundary conditions
        self._remyelinate()                                          # myelin dynamics
        self._apoptosis()
        self._demyelinate()
        self._interact()                                             # C1 + C2 -> C3
        self._renew_c1(n)                                            # relapse-driven birth
        self._deaths()


def simulate(days: float = DEFAULT_DAYS, treat: float = 0.0, immuno: float = 0.0,
             seed: int = 0, n_steps: int | None = None,
             myelin_width: int | None = None, height: int | None = None,
             record_every: int = 1, bias_refresh_steps: int = 6,
             **_ignored) -> np.ndarray:
    """Run the grounded ABM; return mean myelin damage over time, in [0, 1].

    `days` is simulated time at the published 20-minute timestep. `n_steps`
    overrides it directly (kept so older callers that passed step counts still
    work). Rates are per-timestep and are never rescaled.
    """
    steps = int(n_steps) if n_steps is not None else int(round(days * STEPS_PER_DAY))
    steps = max(steps, 1)
    model = WeatherleyABM(treat=treat, immuno=immuno, seed=seed,
                          myelin_width=myelin_width, height=height,
                          bias_refresh_steps=bias_refresh_steps)
    out = []
    for n in range(1, steps + 1):
        model.step(n)
        if n % record_every == 0 or n == steps:
            out.append(model.damage_fraction())
    return np.asarray(out, dtype=float)


class ABMBrick:
    """Stage: run the grounded population model, honouring the intervention."""

    name = "abm:weatherley-ms-oligodendrocyte(numpy-port)"

    def __init__(self, days: float | None = None, seed: int = 0,
                 n_steps: int | None = None, record_every: int = 12,
                 profile: str = "gate", myelin_width: int | None = None,
                 height: int | None = None) -> None:
        if profile not in PROFILES:
            raise KeyError(f"unknown profile {profile!r}; have {sorted(PROFILES)}")
        cfg = PROFILES[profile]
        self.profile = profile
        self.days = cfg["days"] if days is None else days
        self.myelin_width = cfg["myelin_width"] if myelin_width is None else myelin_width
        self.height = cfg["height"] if height is None else height
        self.seed = seed
        self.n_steps = n_steps
        self.record_every = record_every

    def run(self, state: dict) -> dict:
        treat = 0.0
        immuno = 0.0
        interv = state.get("intervention")
        if isinstance(interv, dict):
            treat = float(interv.get("treat", 0.0))
            immuno = float(interv.get("immunogenic", 0.0))
        seed = int(state.get("seed", self.seed))     # per-patient draw from B9
        damage = simulate(days=self.days, treat=treat, immuno=immuno, seed=seed,
                          n_steps=self.n_steps, record_every=self.record_every,
                          myelin_width=self.myelin_width, height=self.height)
        state["abm_damage"] = damage
        state["abm_meta"] = {
            "validated": False,
            "engine": self.name,
            "grounded_in": CITATION,
            "grounding_status": "rates+rules from the published model; NO paper "
                                "figure reproduced yet, so not validated",
            "final_damage": float(damage[-1]),
            "treat": treat,
            "immuno": immuno,
            "bbb_leaving_prob": bbb_permeability(treat, immuno),
            "published_bbb_baseline": W["leaving_prob"],
            "published_bbb_treated": W["therapeutic_leaving_prob"],
            "profile": self.profile,
            "days_simulated": self.days if self.n_steps is None else None,
            "published_days": W["published_days"],
            "lattice": (self.myelin_width, self.height),
            "published_lattice": (W["myelin_width"], W["domain_height"]),
            "profile_note": ("published lattice and published rates in every "
                             "profile; only the integration horizon differs. "
                             "ABSOLUTE damage is not comparable across profiles; "
                             "relative treatment effect is (measured: within "
                             "~8pp between 'gate' and 'default')"),
            "seed": seed,
        }
        return state

    __call__ = run


if __name__ == "__main__":
    print("B5 ABM — grounded in the Weatherley MS model (validated=False).")
    print(f"  {CITATION}\n")
    print(f"  BBB permeability: baseline {W['leaving_prob']}, "
          f"paper's therapy {W['therapeutic_leaving_prob']} "
          f"(our treat=0.75 -> {bbb_permeability(0.75):.3f})\n")
    untreated = simulate(treat=0.0, seed=0)
    treated = simulate(treat=0.8, seed=0)
    harmful = simulate(treat=0.0, immuno=0.4, seed=0)
    print(f"  final mean myelin damage over {DEFAULT_DAYS} days (lower = healthier):")
    print(f"    untreated        {untreated[-1]:.4f}")
    print(f"    treated (0.8)    {treated[-1]:.4f}")
    print(f"    immunogenic(0.4) {harmful[-1]:.4f}")
    print(f"\n  treatment reduces damage? {treated[-1] < untreated[-1]}")
    print(f"  immunogenic worsens it?   {harmful[-1] > untreated[-1]}")
    print("\n  Rates are the published ones; no paper figure reproduced yet.")
