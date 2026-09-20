"""B4 — QSP brick, GROUNDED in the Vélez de Mendizábal MS cross-regulation model.

WHAT CHANGED AND WHY
--------------------
`bricks/qsp.py` is a 3-species toy (A, C, M) with hand-picked rates. BUILD_PLAN
§4 budgeted 20 minutes for it and said so; §8 blocker (1) records that the
August finding "no open MS QSP model exists" was too strong.

This module is a transcription of the one published, open-access, MS-specific
ODE model whose complete source the authors shipped:

    Vélez de Mendizábal N, Carneiro J, Solé RV, Goñi J, Bragard J,
    Martinez-Forero I, Martinez-Pasamar S, Sepulcre J, Torrealdea J, Bagnato F,
    Garcia-Ojalvo J, Villoslada P.
    "Modeling the effector - regulatory T cell cross-regulation reveals the
     intrinsic character of relapses in Multiple Sclerosis."
    BMC Systems Biology 2011;5:114.  PMID 21762505.
    doi:10.1186/1752-0509-5-114   PMC3155504

Two sources were used, and they disagree in three places (see DISCREPANCIES):

  (a) Table 1 of the paper — the parameter table as published.
  (b) Additional file 2, `1752-0509-5-114-S2.MDL` — the authors' own **Vensim
      model file**, i.e. the executable artifact behind the figures. Fetched
      from the Europe PMC supplementary package for PMC3155504.

Where they differ this module follows **(b)**, because (b) is what produced the
paper's figures — except where the paper's own body text states the intended
quantity outright, which it does once (DISCREPANCY 3, and there the MDL is
simply wrong). Every difference is recorded rather than resolved quietly.

NOTHING HERE IS TUNED. No rate below was chosen, fitted, or adjusted by this
repo. Every one is annotated with where it came from.

THE MODEL
---------
Six states. Two T-cell populations, each in a resting and an activated pool,
plus two damage compartments.

  Er  resting (naive/memory) effector T cells        cells
  E   activated effector T cells                     cells
  Rr  resting regulatory T cells                     cells
  R   activated regulatory T cells                   cells
  l   reversible damage                              arbitrary
  L   irreversible damage                            arbitrary

The cross-regulation — the scientific point of the paper — is two Hill terms
that run in opposite directions:

  Treg INHIBIT effector proliferation     kr^h / (kr^h + R^h)
  Treg PROMOTE effector death/egress      R^h  / (kr^h + R^h)
  Teff RECRUIT regulatory cells           E^h  / (ke^h + E^h)

So effectors call in their own regulators, and regulators both slow effector
growth and clear them. That negative feedback is what makes the system settle,
and what makes it *oscillate* rather than run away or flatline.

  dEr/dt = naive_E(t) - delta*Er - beta*Er + eta*E
  dRr/dt = naive_R(t) - delta*Rr - beta*Rr + eta*R
  dE/dt  = delta*Er + alpha_E*(kr^h/(kr^h+R^h))*E - gamma_E*(R^h/(kr^h+R^h))*E - eta*E
  dR/dt  = delta*Rr + alpha_R*(E^h/(ke^h+E^h))*R - gamma_R*R - eta*R
  dl/dt  = d1*(E/a)^n - d2*l - r*l
  dL/dt  = d2*l

RELAPSES ARE NOT SCHEDULED — and that is the whole reason to port this
-----------------------------------------------------------------------
The toy model this replaces had no relapses at all; damage was a smooth decay
and "a relapse" was whatever the readout scaled it into. Here, `naive_E` and
`naive_R` are a **stochastic input** (Additional file 2): at each step, with
probability `100*dt/365`, an influx of 100 cells arrives. That is the paper's
model of ordinary immune events — infections — arriving at ~100/year.

Each influx perturbs the Teff/Treg balance. The cross-regulation loop pulls it
back. Between the kick and the recovery the effector population spikes, and
`dl/dt` turns that spike into damage. The paper: "The presence of activated
Teff cell peaks produces both reversible and irreversible damage. The sum of
both represents the clinical relapses."

So a relapse is an EMERGENT event here, produced by a stochastic kick against a
regulated loop, not a scheduled one. That is what gives this brick something
the toy could not have: a candidate intervention can change *when and whether*
relapses happen, not just scale a damage number.

INTERVENTION — the candidate vocabulary (this repo's addition, marked as such)
------------------------------------------------------------------------------
Everything above is transcription. This paragraph is not, and is flagged so
nobody mistakes it for the paper.

The published model has no drug term. This module adds one in the only way that
does not invent biology: **an intervention is a set of multipliers on the
model's own named parameters.** Nothing new is introduced; a drug can only turn
the dials the model already has.

    InterventionPoint    multiplies   meaning for a candidate
    ------------------   ----------   ------------------------------------------
    alpha_E              alpha_E      antiproliferative on effector T cells
                                      (teriflunomide, dimethyl fumarate)
    gamma_E              gamma_E      raises effector death / egress / migration
                                      (alemtuzumab, ocrelizumab, natalizumab,
                                       S1P modulators)
    alpha_R              alpha_R      Treg support (>1) or Treg damage (<1)
    gamma_R              gamma_R      Treg clearance (>1 = removes regulation)
    delta                delta        antigen presentation / activation rate
    naive_E              naive influx of effectors (immunogenic agents raise it)
    ke                   ke           the Teff level at which Treg recruitment is
                                      half-maximal -- LOWER means effectors call in
                                      regulation sooner. This is the one point whose
                                      use is not just plausible but FITTED: see below.

`mult=1.0` is untreated, everywhere. A pure suppressive drug is
`{"alpha_E": 0.5}`; a depleting one is `{"gamma_E": 2.0}`.

**`ke` comes from the successor paper, with a number attached.**
Martinez-Pasamar et al. 2013 (*BMC Syst Biol* 7:34, PMC3651362) reuse these
exact equations ("the equations of the T-cell cross-regulation model as
described in [10]") and add no compartments and no B cells -- so there is
nothing structural to port from it. What it does add is a sensitivity analysis
against EAE flow-cytometry data, and one result is directly usable:

    "the dynamics of the antigen-specific T-cell subpopulation after anti-CD20
     therapy was reproduced by reducing the K_eff threshold below the healthy
     standard (<1,000 cells; e.g. 850 cells), independently of the alpha_reg
     parameter"

    "B-cell depletion therapy may influence the autoimmune process by preventing
     uncontrolled activation of T_eff without strengthening T_reg activation"

That is a mechanism for a B-cell-depleting drug in a model with no B cells,
fitted to mouse T-cell dynamics and not to any human relapse rate. It is why
`ke` exists as a dial here. Note what it is NOT: a licence to route every
B-lineage agent through `ke`. The result is specific to anti-CD20.

**This is what makes harm expressible mechanistically for the first time.** The
old two-constant scheme could only say "suppressive" or "immunogenic", so
lenercept — immunosuppressive AND harmful — had to reuse the immunogenic
constant (see bricks/grounding.py). Here it is one object: a drug that lowers
`alpha_E` *and* raises `gamma_R` suppresses the attack while stripping out the
regulation that contains it. Whether that nets to harm is then a PREDICTION of
the dynamics, not a constant anyone picked.

DISCREPANCIES between Table 1 and the authors' model file
----------------------------------------------------------
Recorded, not resolved by preference. For 1 and 2 this module uses the MDL
value and the Table 1 value is available via `params=`/`y0=`. For 3 the paper's
body text decides it, against the MDL.

  1. **d2 (irreversible damage rate).** Table 1: `0.02 day-1`.
     MDL: `d2 = 0.002`. Tenfold. This module uses 0.002.
  2. **Resting-pool initial conditions.** Table 1: `En0 = 0`, `Rn0 = 0`.
     MDL: `Resting Te = 7.5`, `Resting Tr = 2.4`. This module uses the MDL's,
     because starting both resting pools at 0 is inconsistent with the naive
     influx being the model's only source of new cells.
  3. **The damage drive — RESOLVED, and the MDL is the one that is wrong.**
     Table 1 does not list the exponent. The MDL defines `n = 2` and
     `eprima = (Activated Te^n)/A`, but the authors' own inline comment on that
     same variable reads `Activated Te^n/threshold^n`. Those differ by a factor
     of A = 22800, so it mattered. The paper's body text settles it, in the
     sentence introducing equations (7) and (8):

         "E' = (E/a)^2 is the second order effect of Teff cells on damage"

     So the intended quantity is `(E/a)^2`, the MDL's comment, and the MDL's
     *equation* is a transcription error in the authors' own file. This module
     therefore defaults to `damage_form="paper"`. The literal MDL equation is
     kept available as `damage_form="mdl_literal"` so the difference stays
     inspectable rather than becoming folklore.

     (The two differ by a constant factor on both damage compartments, so they
     cannot change a trajectory's shape or an arm-vs-arm ratio — which is what
     the clinical gate scores. It would change any absolute damage claim.)

  Also noted, not a discrepancy: Table 1 gives alpha_E and alpha_R as RANGES
  (`[1:2]` and `[0.25:2]`), being sweep ranges. The MDL's chosen operating
  point is `alpha_E = 2`, `alpha_R = 0.25`, and that is what this module uses.

WHAT THIS MODEL CANNOT REPRESENT: A DEPLETING THERAPY AS BENEFICIAL
--------------------------------------------------------------------
Measured 2026-09-20, and it is a structural property of the published model
rather than a defect in this port. Damage is driven by `(E/a)^2`, so it is set
by effector PEAK EXCURSIONS, not by the effector median. Over 48 histories at
730 days:

    arm                        med E    med PEAK E    damage
    untreated                   1126         52741     1.457
    alpha_E x0.5 (damp growth)  1020         30794     0.665
    gamma_E x1.5 (kill cells)   1077        192455    10.249

The median effector load is almost identical in all three. The PEAK differs
six-fold, and damage follows the peak.

The mechanism: effectors recruit their own regulators through `E^h/(ke^h+E^h)`.
Killing effectors lowers that recruitment, `R` falls, the proliferation brake
`kr^h/(kr^h+R^h)` releases, and the population rebounds into a LARGER excursion
than the one the killing removed. Damping proliferation (`alpha_E` down) instead
quiets the oscillator and does help.

So in this model **removing effector cells is self-defeating**, whatever channel
does the removing. This was tested directly: an additive,
regulation-independent loss term `-depletion_rate * E` was added as an explicit
extension and it did NOT fix the sign (median damage 1.46 -> 4.69 -> 11.3 ->
39.8 as the term rose, with in-regime runs falling from 64 to 26). The extension
was removed rather than tuned, which is what its own comment said to do.

Consequence for any screen over this model: the entire depleting / sequestering
/ trafficking-blocking class -- natalizumab, fingolimod, ponesimod, alemtuzumab,
atacicept -- cannot come out beneficial here. That is not a calibration gap. It
is the reason the grounded direction gate scores 5/13 against the ABM path's
9/14 (`backtest/clinical_velez.py`).

BUILT != VALIDATED
------------------
`validated=False`, deliberately, exactly as bricks/abm.py is after the
Weatherley port. Transcribing a published model is not the same as reproducing
its results. What is true today: the equations and rates are the authors'; the
integrator is this repo's; no figure of the paper has been reproduced. Moving
that flag requires a run that reproduces a published figure, and that has not
been done.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------- #
# Parameters. Every value from Additional file 2 (the authors' Vensim model),
# cross-checked against Table 1. `T1:` marks the Table 1 value where it differs.
# --------------------------------------------------------------------------- #

VELEZ_PARAMS: dict[str, float] = {
    # --- activation / turnover -------------------------------------------- #
    "delta":   1.0,     # antigen presentation, day^-1        MDL AntigenPresentation=1; T1 delta=1 (ref 35)
    "beta":    0.01,    # anergy of resting cells, day^-1     MDL TeAnergy/TrAnergy = Resting*0.01; T1 beta=0.01
    "eta":     0.01,    # activated -> memory, day^-1         MDL TeMemory/TrMemory = Activated*0.01; T1 eta=0.01
    "alpha_E": 2.0,     # max Teff proliferation, day^-1      MDL Max Te Proliferation Rate=2; T1 [1:2] (ref 36)
    "alpha_R": 0.25,    # max Treg activation+prolif, day^-1  MDL Max Tr activation and proliferation Rate=0.25; T1 [0.25:2]
    "gamma_E": 0.2,     # max Teff death/anergy/migration     MDL Max Te Death Anergy Migration Rate=0.2; T1 0.2 (ref 37)
    "gamma_R": 0.2,     # Treg death/anergy/migration         MDL Tr Death Anergy Migration Rate=0.2; T1 0.2
    # --- cross-regulation -------------------------------------------------- #
    "ke":      1000.0,  # Teff half-maximal for Treg recruitment, cells   MDL ke=1000; T1 Ke=1000
    "kr":      200.0,   # Treg half-maximal for Teff control, cells       MDL kr=200;  T1 Kr=200
    "h":       5.0,     # Hill coefficient                                MDL h=5;     T1 h=5
    # --- damage ------------------------------------------------------------ #
    "d1":      1.0,     # reversible damage rate, day^-1      MDL d1=1;    T1 d1=1
    "d2":      0.002,   # reversible -> irreversible, day^-1  MDL d2=0.002; T1 says 0.02  <-- DISCREPANCY 1
    "r":       0.1,     # natural recovery, day^-1            MDL r=0.1;   T1 r=0.1
    "A":       22800.0, # damage threshold, cells             MDL A=22800; T1 a=22800
    "n":       2.0,     # damage exponent                     MDL n=2;     T1 does not list it
    # --- stochastic immune events ------------------------------------------ #
    "naive_events_per_year": 100.0,  # MDL: P(step) = 100*TIME STEP/365
    "naive_influx_cells":    100.0,  # MDL: amount 100/TIME STEP for one step = 100 cells
}

# Initial conditions. MDL INTEG(...) second argument.
VELEZ_Y0: dict[str, float] = {
    "Er": 7.5,     # MDL Resting Te INTEG(..., 7.5);   T1 says En0 = 0   <-- DISCREPANCY 2
    "Rr": 2.4,     # MDL Resting Tr INTEG(..., 2.4);   T1 says Rn0 = 0   <-- DISCREPANCY 2
    "E":  1000.0,  # MDL Activated Te INTEG(..., 1000); T1 Ea0 = 1000
    "R":  200.0,   # MDL Activated Tr INTEG(..., 200);  T1 Ra0 = 200
    "l":  0.0,     # MDL RevesibleDamage INTEG(..., 0); T1 l0 = 0
    "L":  0.0,     # MDL IrreversibleDamage INTEG(..., 0); T1 L0 = 0
}

# Simulation control, MDL .Control section.
VELEZ_DT = 0.1          # TIME STEP = 0.1 Day
VELEZ_T_END = 1825.0    # FINAL TIME = 1825 Day (5 years)

# The parameters an intervention is allowed to touch. A candidate that wants to
# act anywhere else needs a model that represents that place — which is the
# honest failure mode, and the reason this list is explicit and short.
INTERVENTION_POINTS: tuple[str, ...] = (
    "alpha_E", "gamma_E", "alpha_R", "gamma_R", "delta", "naive_E", "ke",
)


@dataclass(frozen=True)
class MechanismProfile:
    """A candidate expressed as multipliers on the published model's own dials.

    1.0 means "untouched". This is the whole drug representation: no new
    biology, only the rates Vélez de Mendizábal already defined.

    `label` is free text for reporting. `source` should say where the numbers
    came from — a citation for a real drug, or "hypothetical" for a screened
    candidate. It is carried into the output so a trajectory can always be
    traced back to what produced it.
    """

    label: str = "untreated"
    alpha_E: float = 1.0
    gamma_E: float = 1.0
    alpha_R: float = 1.0
    gamma_R: float = 1.0
    delta: float = 1.0
    naive_E: float = 1.0
    ke: float = 1.0
    source: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        for point in INTERVENTION_POINTS:
            v = getattr(self, point)
            if not np.isfinite(v) or v < 0.0:
                raise ValueError(f"{self.label}: {point} must be finite and >= 0, got {v}")

    @property
    def is_untreated(self) -> bool:
        return all(getattr(self, p) == 1.0 for p in INTERVENTION_POINTS)

    def as_dict(self) -> dict:
        return {"label": self.label, "source": self.source,
                **{p: getattr(self, p) for p in INTERVENTION_POINTS}}


UNTREATED_PROFILE = MechanismProfile(label="untreated", source="control arm")


def _hill(x: float, k: float, h: float) -> float:
    """x^h / (k^h + x^h). Guarded at both ends.

    x=0 so that h<1 cannot produce a 0/0, and large x because this model has no
    carrying capacity on E: a run that has left its regime can reach values
    where x**h overflows a float before any finiteness check downstream sees it.
    The limit is 1, so returning it is exact, not a fudge.
    """
    if x <= 0.0:
        return 0.0
    if x > k * 1e12:
        return 1.0
    xh = x ** h
    return xh / (k ** h + xh)


def simulate(
    profile: MechanismProfile = UNTREATED_PROFILE,
    params: dict | None = None,
    *,
    t_end: float = VELEZ_T_END,
    dt: float = VELEZ_DT,
    y0: dict | None = None,
    seed: int | None = None,
    damage_form: str = "paper",
    regime_ceiling_mult: float = 100.0,
) -> dict:
    """Integrate the Vélez de Mendizábal 2011 model under one intervention.

    Fixed-step explicit Euler at `dt`, because that is what the authors' Vensim
    model does (TIME STEP = 0.1, Euler) and because the naive-cell influx is a
    per-step discrete event — an adaptive solver would step over it. Fidelity to
    the published integration beats a nicer integrator here.

    `damage_form="paper"` uses `(E/a)^n`, the quantity the paper's body text
    defines. `"mdl_literal"` uses `E^n/a`, the authors' Vensim equation as
    written. See DISCREPANCIES 3.

    `seed` controls the stochastic immune events only. Two arms run on the same
    seed see the same infection history. That sounds like it should be what
    makes an arm-vs-arm ratio meaningful, and it is NOT: measured, pairing
    narrows the ratio's IQR by only ~1.1x
    (`scripts/measure_qsp_variance.py`). What makes a comparison valid here is
    cohort size -- this model's damage spans ~90,000-fold across histories, so
    anything comparing arms on a handful of runs is wrong by construction. See
    the SEEDS comment in backtest/lomo.py.

    `regime_ceiling_mult` guards the model's one structural weakness. There is
    no carrying capacity on E — the *only* thing bounding the effector
    population is the Treg feedback. Strip enough regulation (a large `gamma_R`
    or a small `alpha_R`) and `dE/dt -> (alpha_E - eta)*E`, which is unbounded
    exponential growth. That is a real property of the published model, not a
    porting error, but the resulting numbers are not "extreme harm" — they are
    the model outside the regime anything was ever said about. So a run whose E
    exceeds `regime_ceiling_mult * a` is flagged `in_regime: False`, integration
    HALTS at that point, and its damage must not be ranked against an in-regime
    arm. The value is not clamped to something plausible: silently bounding it
    would hide exactly the thing a screen needs to see.
    """
    p = {**VELEZ_PARAMS, **(params or {})}
    if damage_form not in ("paper", "mdl_literal"):
        raise ValueError(
            f"damage_form must be 'paper' or 'mdl_literal', got {damage_form!r}")

    # Intervention multipliers applied to the published rates. This is the only
    # place a drug touches the model.
    delta = p["delta"] * profile.delta
    alpha_E = p["alpha_E"] * profile.alpha_E
    gamma_E = p["gamma_E"] * profile.gamma_E
    alpha_R = p["alpha_R"] * profile.alpha_R
    gamma_R = p["gamma_R"] * profile.gamma_R

    beta, eta, h = p["beta"], p["eta"], p["h"]
    ke, kr = p["ke"] * profile.ke, p["kr"]
    d1, d2, r, A, n = p["d1"], p["d2"], p["r"], p["A"], p["n"]

    y = {**VELEZ_Y0, **(y0 or {})}
    # `l` and `L` are the paper's own symbols for reversible and irreversible
    # damage (equations 7 and 8). E741 is silenced rather than renaming them,
    # because matching the source's notation is the point of a transcription.
    Er, Rr, E, R, l, L = y["Er"], y["Rr"], y["E"], y["R"], y["l"], y["L"]  # noqa: E741

    n_steps = int(round(t_end / dt))
    rng = np.random.default_rng(seed)
    regime_ceiling = regime_ceiling_mult * A
    left_regime_at: float | None = None

    # MDL: IF RANDOM UNIFORM(0,1) < (100*TIME STEP/365) THEN 100/TIME STEP ELSE 0
    p_event = p["naive_events_per_year"] * dt / 365.0
    influx_rate = p["naive_influx_cells"] / dt   # a rate held for exactly one step

    ts = np.empty(n_steps + 1)
    arr = {k: np.empty(n_steps + 1) for k in ("Er", "Rr", "E", "R", "l", "L")}

    def _record(i, t):
        ts[i] = t
        arr["Er"][i], arr["Rr"][i] = Er, Rr
        arr["E"][i], arr["R"][i] = E, R
        arr["l"][i], arr["L"][i] = l, L

    _record(0, 0.0)

    for i in range(1, n_steps + 1):
        # --- stochastic immune events (infections), MDL Naive Te/Tr Input --- #
        naive_E_in = influx_rate * profile.naive_E if rng.random() < p_event else 0.0
        naive_R_in = influx_rate if rng.random() < p_event else 0.0

        # --- cross-regulation Hill terms --------------------------------- #
        treg_inhibits_prolif = 1.0 - _hill(R, kr, h)   # MDL kr^h/(kr^h+R^h)
        treg_drives_death = _hill(R, kr, h)            # MDL R^h/(kr^h+R^h)
        teff_recruits_treg = _hill(E, ke, h)           # MDL E^h/(ke^h+E^h)

        # --- derivatives -------------------------------------------------- #
        dEr = naive_E_in - delta * Er - beta * Er + eta * E
        dRr = naive_R_in - delta * Rr - beta * Rr + eta * R
        dE = (delta * Er
              + alpha_E * treg_inhibits_prolif * E
              - gamma_E * treg_drives_death * E
              - eta * E)
        dR = (delta * Rr
              + alpha_R * teff_recruits_treg * R
              - gamma_R * R
              - eta * R)

        drive = (E / A) ** n if damage_form == "paper" else (E ** n) / A
        dl = d1 * drive - d2 * l - r * l
        dL = d2 * l

        Er = max(0.0, Er + dt * dEr)
        Rr = max(0.0, Rr + dt * dRr)
        E = max(0.0, E + dt * dE)
        R = max(0.0, R + dt * dR)
        l = max(0.0, l + dt * dl)  # noqa: E741  -- paper's symbol, see above
        L = max(0.0, L + dt * dL)

        out_of_regime = E > regime_ceiling or not np.isfinite(E + R + l + L)
        if out_of_regime and left_regime_at is None:
            left_regime_at = i * dt

        if out_of_regime:
            # Integration STOPS here. Past the regime ceiling the effector
            # population is growing without bound, so every further step is
            # arithmetic about a state the model was never claimed to describe
            # -- and continuing costs a float overflow inside the Hill terms.
            # The remaining samples are held at the last value so the array
            # shape stays contractual; `in_regime=False` is what a consumer
            # must branch on, never the damage number.
            _record(i, i * dt)
            for key in arr:
                arr[key][i:] = arr[key][i]
            ts[i:] = np.arange(i, n_steps + 1) * dt
            break

        _record(i, i * dt)

    return {
        "t": ts,
        **arr,
        "total_damage": arr["l"] + arr["L"],
        "profile": profile.as_dict(),
        "damage_form": damage_form,
        "seed": seed,
        "in_regime": left_regime_at is None,
        "left_regime_at": left_regime_at,
        "validated": False,
        "note": ("Vélez de Mendizábal 2011 (BMC Syst Biol 5:114) transcribed from the "
                 "authors' Additional file 2 Vensim model; rates not tuned. Ported, "
                 "NOT reproduced — no figure of the paper has been checked."),
    }


def relapse_events(traj: dict, *, threshold_mult: float = 2.0,
                   min_separation_days: float = 30.0,
                   baseline: float | None = None) -> list[float]:
    """Times (days) of emergent effector peaks — this model's relapses.

    A relapse is a local maximum of activated effector cells exceeding
    `threshold_mult * baseline`, with peaks closer than `min_separation_days`
    merged.

    **`baseline` is the load-bearing argument and it must be the UNTREATED
    arm's median E when comparing arms.** Left None it falls back to the
    trajectory's own median, which is only correct for describing a single run:
    a drug that lowers effector levels also lowers its own median, so a
    self-referenced threshold slides down with the treatment and reports an
    unchanged relapse count for an arm whose effector load fell tenfold. Use
    `baseline_from(untreated_traj)` and pass it to every arm in the comparison.

    **The detector is this repo's, not the paper's.** Vélez de Mendizábal
    identify relapses visually from the damage curve ("the presence of activated
    Teff cell peaks produces both reversible and irreversible damage"); they
    publish no numeric criterion. So the threshold is a reading convention over
    the paper's trajectory, and its constants are stated rather than hidden.
    Anything reporting a relapse RATE inherits that choice and must say so.
    """
    E = np.asarray(traj["E"], dtype=float)
    t = np.asarray(traj["t"], dtype=float)
    if E.size < 3:
        return []
    ref = float(np.median(E)) if baseline is None else float(baseline)
    if ref <= 0.0:
        return []
    cut = ref * threshold_mult

    peaks: list[float] = []
    for i in range(1, E.size - 1):
        if E[i] > cut and E[i] >= E[i - 1] and E[i] > E[i + 1]:
            if peaks and (t[i] - peaks[-1]) < min_separation_days:
                continue
            peaks.append(float(t[i]))
    return peaks


def baseline_from(untreated_traj: dict) -> float:
    """The reference effector level every arm in a comparison must be scored on.

    Always the UNTREATED arm's median E. See `relapse_events`.
    """
    return float(np.median(np.asarray(untreated_traj["E"], dtype=float)))


class VelezQSPBrick:
    """Stage: integrate the grounded MS QSP model, honoring state["intervention"].

    Writes `state["qsp_traj"]` (same key the toy brick wrote, so the spine and
    the run report are unchanged) plus `state["qsp_damage"]`, the total-damage
    trajectory, which is the quantity a readout should consume.

    **Per-patient variation is currently the infection history only.**
    `bricks/vpop.py` samples per-patient values for the TOY model's parameters
    (`r_CA`, `k_dmg`), which do not exist in this model. They are dropped and
    listed in `qsp_traj["ignored_params"]` rather than silently absorbed. So a
    cohort run through this brick varies by `seed` and by nothing else, and
    vpop's Allen-Rieger plausibility filter — built on toy parameters — does not
    constrain it. Porting that filter onto Vélez parameters is open work; until
    then "virtual patient" means "one stochastic infection history".

    Reads `state["mechanism_profile"]` (a MechanismProfile or a plain dict of
    multipliers) if present. Falls back to `state["intervention"]`'s scalar
    `treat`/`immunogenic` so the existing pipeline keeps running during the
    migration: `treat` becomes an alpha_E multiplier of (1 - treat) and
    `immunogenic` a naive_E multiplier of (1 + immunogenic). **That fallback is
    a shim, not science** — it reproduces the old brick's single-axis behaviour
    inside the new model so nothing breaks mid-port, and it should disappear
    once every arm carries a real profile.
    """

    name = "qsp:velez2011-teff-treg-cross-regulation"
    requires: tuple[str, ...] = ()

    def __init__(self, params: dict | None = None, t_end: float = VELEZ_T_END,
                 dt: float = VELEZ_DT, damage_form: str = "paper") -> None:
        self.params = params
        self.t_end = t_end
        self.dt = dt
        self.damage_form = damage_form

    def _profile_from_state(self, state: dict) -> MechanismProfile:
        prof = state.get("mechanism_profile")
        if isinstance(prof, MechanismProfile):
            return prof
        if isinstance(prof, dict):
            known = {k: v for k, v in prof.items()
                     if k in INTERVENTION_POINTS or k in ("label", "source")}
            return MechanismProfile(**known)

        interv = state.get("intervention")
        if isinstance(interv, dict):
            treat = float(interv.get("treat", 0.0))
            immuno = float(interv.get("immunogenic", 0.0))
            if treat or immuno:
                return MechanismProfile(
                    label=str(interv.get("name", "shim")),
                    alpha_E=max(0.0, 1.0 - treat),
                    naive_E=1.0 + immuno,
                    source="SHIM from scalar treat/immunogenic — not a mechanism",
                )
        return UNTREATED_PROFILE

    def run(self, state: dict) -> dict:
        profile = self._profile_from_state(state)

        # state["qsp_params"] is per-patient variation from bricks/vpop.py, and
        # vpop still speaks the TOY model's parameter names (r_CA, k_dmg). Merging
        # those into VELEZ_PARAMS would add dead keys that change nothing, so a
        # cohort would silently collapse to one identical patient with no error
        # anywhere. Unknown names are dropped and REPORTED instead.
        incoming = {**(self.params or {}), **(state.get("qsp_params") or {})}
        known = {k: v for k, v in incoming.items() if k in VELEZ_PARAMS}
        ignored = sorted(set(incoming) - set(known))

        traj = simulate(
            profile,
            params=known,
            t_end=self.t_end,
            dt=self.dt,
            seed=state.get("seed"),
            damage_form=self.damage_form,
        )
        traj["ignored_params"] = ignored
        state["qsp_traj"] = traj
        state["qsp_damage"] = traj["total_damage"]
        state["qsp_relapses"] = relapse_events(traj)
        return state

    __call__ = run


if __name__ == "__main__":
    print("B4 QSP — Vélez de Mendizábal 2011 (BMC Syst Biol 5:114), ported.\n")

    # --- the paper's own Figure 3 result, as a reproduction check ----------- #
    # "by decreasing the maximum activation and proliferation rate of Treg
    #  (alpha_R) ... immune homeostasis was lost and spontaneous immune
    #  responses emerged" -- healthy (3A/3C) vs autoimmune (3B/3D).
    print("Figure 3 check — alpha_R is the paper's health/autoimmunity axis.")
    fig3_ref = baseline_from(simulate(UNTREATED_PROFILE, seed=1))
    print(f"{'alpha_R':>10} {'mean E':>12} {'relapses':>9} {'total damage':>14} {'regime':>8}")
    for mult in (8.0, 4.0, 2.0, 1.0):
        tr = simulate(MechanismProfile(label=f"alpha_R x{mult}", alpha_R=mult,
                                       source="Figure 3 sweep"), seed=1)
        print(f"{VELEZ_PARAMS['alpha_R'] * mult:>10.2f} {tr['E'].mean():>12.1f} "
              f"{len(relapse_events(tr, baseline=fig3_ref)):>9d} {tr['total_damage'][-1]:>14.3f} "
              f"{'ok' if tr['in_regime'] else 'OUT':>8}")
    print("  -> higher alpha_R should mean fewer effector peaks and less damage;\n"
          "     the MDL's operating point (0.25) is the autoimmune configuration.\n")

    # --- what the intervention vocabulary buys ----------------------------- #
    print("Intervention points, same infection history (seed=1), 5 years:")
    arms = (
        UNTREATED_PROFILE,
        MechanismProfile(label="antiproliferative (alpha_E x0.9)", alpha_E=0.9,
                         source="illustrative, not a real drug"),
        MechanismProfile(label="depleting (gamma_E x1.5)", gamma_E=1.5,
                         source="illustrative, not a real drug"),
        MechanismProfile(label="Treg support (alpha_R x2)", alpha_R=2.0,
                         source="illustrative, not a real drug"),
        MechanismProfile(label="suppress + strip regulation", alpha_E=0.9, gamma_R=1.5,
                         source="illustrative, not a real drug"),
    )
    ref = baseline_from(simulate(UNTREATED_PROFILE, seed=1))
    for prof in arms:
        tr = simulate(prof, seed=1)
        flag = "" if tr["in_regime"] else "   <-- OUT OF REGIME, damage not comparable"
        print(f"  {prof.label:<34} damage = {tr['total_damage'][-1]:12.3f}   "
              f"relapses = {len(relapse_events(tr, baseline=ref)):3d}{flag}")
    print("\n  -> the last arm suppresses effectors AND strips regulation. The model has")
    print("     no carrying capacity on E, so that runs away rather than 'harming':")
    print("     a real limit of the published model, reported instead of ranked.")
    print("     (ported, NOT reproduced — validated=False)")
