"""Probe A — the two-equation myelin/inflammation model of Jenner et al. 2026.

    Jenner AL, Weatherley G, Frascoli F. "A minimal mathematical model of
    demyelination and inflammation in multiple sclerosis."
    J R Soc Interface 2026;23:20260364. PMID 42481015. Open access.
    Sections 2-6 vendored at docs/research/jenner2026/sections_2_to_6.md.

Equation (2.1), transcribed exactly, time unit one month:

    dM/dt = r M (1 - M) - phi M I / (M + eta)
    dI/dt = phi M I / (M + eta) - delta I

    M     proportion of viable healthy myelin, in (0, 1]
    I     level of inflammation (the paper's proxy: serum neurofilament light)
    r     remyelination rate, "between 0.3 and 1.0"           (section 2)
    delta decay rate of inflammation, "chosen around unity"   (section 2)
    eta   half-effect myelin level, "from zero to unity"      (section 2)
    phi   disease strength, "similar to the rate of infection" (section 2)

WHY THIS MODEL IS BEING SCORED
------------------------------
docs/FOUR_DAY_PLAN.md, day 1. The Velez port's damage is peak-driven and its
effectors recruit their own brake, so removing effectors raises damage
(docs/GATE_GAP_ANALYSIS.md G3). Here damage is load-driven by construction:
the only thing that demyelinates is the term phi M I / (M + eta), so anything
that lowers the inflammatory attack lowers demyelination, monotonically, at
the level of the equations. Whether that survives the dynamics (limit cycles,
and a disease equilibrium whose inflammation is NOT monotone in phi, section
4) is what the probe measures. Nothing in this file is tuned; the four rates
are the paper's, and the untreated operating point is the paper's own
labelled RRMS example.

THE UNTREATED PATIENT
---------------------
Figure 2D: r = 0.5, delta = 0.2, eta = 0.5, phi = 0.7, "stable period
oscillations with a period of t = 30, capturing RRMS patients". That is the
one parameter set the paper itself labels as a relapsing-remitting patient,
so it is the untreated arm. Its landmarks, from equations (3.1) and (3.7):

    branch point  phi* = delta (1 + eta)            = 0.30
    Hopf          phi_H = delta (1 + eta) / (1 - eta) = 0.60

so phi = 0.7 sits on the limit-cycle side of the Hopf, as the figure says.
The paper's other canonical set (figure 3: r = 0.5, eta = 0.7, delta = 1,
phi* = 1.7, phi_H = 5.67) is kept as `FIG3_PARAMS` for the reproduction test
and is not an arm.

THE READOUTS, AND WHY THE HORIZON IS LONG
-----------------------------------------
Two, both defined against the untreated arm run on the same schedule:

  lesion ratio    integral of I over the scoring window, treated / untreated.
                  This is the quantity that goes through bricks/sormani, the
                  same seat total damage occupies on the Velez path.
  relapse proxy   peaks of I above 2x the UNTREATED arm's median I, at least
                  one month apart, per year. That is `qsp_velez.relapse_events`
                  applied to I, with its constants unchanged, so the two models
                  are read by one convention. It is reported, not scored.

The model is deterministic and periodic. A trial samples patients at random
phases of their cycle, so the population mean of a 1-2 year trial equals the
time average over whole cycles; a 24-month window on ONE deterministic run
would instead measure the phase the drug happened to start at. So every arm
is run through a common untreated burn-in (`BURN_IN_MONTHS`, an established
patient on the attractor), then treated for `SCORE_MONTHS`, long enough to
average over several cycles of the slowest arm. Both constants are stated,
neither is fitted, and both are carried into the output.

INTEGRATOR
----------
Fixed-step classical Runge-Kutta (RK4), `DT` months. Fixed-step is what the
plan asks for and what makes the run reproducible to the bit; RK4 rather than
Euler because a limit cycle integrated with forward Euler at any affordable
step drifts in amplitude, and the reproduction test measures amplitude. The
step is checked in tests/test_qsp_minimal.py by halving it.

WHAT THIS MODEL CANNOT EXPRESS — read before trusting any ranking
----------------------------------------------------------------
No compartment, no cell type, no regulator, no antigen step. A drug is a
multiplier on r, phi, eta or delta. So depletion, sequestration, transit
block, anti-proliferation and cytokine blockade are all "less phi" and are
indistinguishable here; bricks/profiles_minimal.py says so on every arm. The
probe's question is only whether the SIGN table of the exam comes out right
once damage is load-driven. Ranking within a dial is gap G4's job.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from bricks.qsp_velez import relapse_events

# --------------------------------------------------------------------------- #
# Parameters, all from the paper.
# --------------------------------------------------------------------------- #
RATES: tuple[str, ...] = ("r", "phi", "eta", "delta")

# Figure 2 caption: r = 0.5, delta = 0.2, eta = 0.5; panel D phi = 0.7,
# "capturing RRMS patients". The untreated arm.
FIG2D_PARAMS: dict[str, float] = {"r": 0.5, "phi": 0.7, "eta": 0.5, "delta": 0.2}
FIG2D_PERIOD_MONTHS = 30.0          # figure 2 caption, panel D

# Figure 3 caption: r = 0.5, eta = 0.7, delta = 1; BP at phi* = 1.7, HB at
# phi_HB = 5.67; periods "from a minimum of about T = 10 months to T = 20".
FIG3_PARAMS: dict[str, float] = {"r": 0.5, "eta": 0.7, "delta": 1.0}
FIG3_PHI_STAR = 1.7
FIG3_PHI_HOPF = 5.67

# Section 2 ranges. Enforced on the untreated arm, reported (not enforced) on
# treated arms: a drug that takes a rate outside the paper's range is a drug
# the paper's analysis does not cover, and the run says so in `in_range`.
RANGES: dict[str, tuple[float, float]] = {
    "r": (0.3, 1.0),        # "a value of r between 0.3 and 1.0 is chosen"
    "eta": (0.0, 1.0),      # "can only assume values from zero to unity"
    "delta": (0.0, math.inf),   # "chosen around unity"; delta = 0 excluded, section 4
    "phi": (0.0, math.inf),     # positive; the disease strength being swept
}

UNTREATED_PARAMS = FIG2D_PARAMS

# Schedule, in months. See the docstring for why the window is long.
DT = 0.01                    # months; halving it moves the period by < 0.1%
BURN_IN_MONTHS = 120.0       # ten years untreated: four figure-2D cycles
SCORE_MONTHS = 240.0         # twenty years treated: eight figure-2D cycles
Y0 = (1.0, 0.01)             # healthy myelin, a small inflammatory seed

# The relapse detector's constants are qsp_velez's. Its separation is in DAYS.
DAYS_PER_MONTH = 365.25 / 12.0


# --------------------------------------------------------------------------- #
# Profiles: a drug is a multiplier on each of the four rates.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MinimalProfile:
    """Multipliers on r, phi, eta, delta. 1.0 means untouched."""

    label: str = "untreated"
    r: float = 1.0
    phi: float = 1.0
    eta: float = 1.0
    delta: float = 1.0
    source: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        for rate in RATES:
            v = getattr(self, rate)
            if not np.isfinite(v) or v < 0.0:
                raise ValueError(f"{self.label}: {rate} must be finite and >= 0, got {v}")

    @property
    def is_untreated(self) -> bool:
        return all(getattr(self, p) == 1.0 for p in RATES)

    def as_dict(self) -> dict:
        return {"label": self.label, "source": self.source,
                **{p: getattr(self, p) for p in RATES}}


UNTREATED_PROFILE = MinimalProfile(label="untreated", source="control arm")


def params_for(profile: MinimalProfile, base: dict[str, float] | None = None) -> dict[str, float]:
    base = UNTREATED_PARAMS if base is None else base
    return {k: base[k] * getattr(profile, k) for k in RATES}


# --------------------------------------------------------------------------- #
# Analysis, equations (3.1) and (3.7). Used by the tests, not by the scorer.
# --------------------------------------------------------------------------- #
def phi_star(p: dict[str, float]) -> float:
    """Branch point: the healthy state loses stability. phi* = delta (1 + eta)."""
    return p["delta"] * (1.0 + p["eta"])


def phi_hopf(p: dict[str, float]) -> float:
    """Hopf locus, equation (3.7): phi_H = delta (1 + eta) / (1 - eta), eta < 1."""
    if not p["eta"] < 1.0:
        return math.inf
    return p["delta"] * (1.0 + p["eta"]) / (1.0 - p["eta"])


def disease_equilibrium(p: dict[str, float]) -> tuple[float, float]:
    """(M*, I*) of equation (3.1); biologically meaningful only for phi > phi*."""
    u = p["phi"] - p["delta"]
    m = p["delta"] * p["eta"] / u
    i = p["r"] * p["eta"] * (u - p["delta"] * p["eta"]) / (u * u)
    return m, i


# --------------------------------------------------------------------------- #
# Integration.
# --------------------------------------------------------------------------- #
def _rhs(m: float, i: float, r: float, phi: float, eta: float, delta: float) -> tuple[float, float]:
    attack = phi * m * i / (m + eta)
    return r * m * (1.0 - m) - attack, attack - delta * i


def integrate(p: dict[str, float], y0: tuple[float, float], t_end: float,
              dt: float = DT) -> dict:
    """RK4 from y0 for t_end months. Returns t, M, I as arrays, plus in_regime."""
    r, phi, eta, delta = p["r"], p["phi"], p["eta"], p["delta"]
    n = int(round(t_end / dt))
    t = np.empty(n + 1)
    M = np.empty(n + 1)
    infl = np.empty(n + 1)
    m, i = float(y0[0]), float(y0[1])
    t[0], M[0], infl[0] = 0.0, m, i
    ok = True
    for k in range(1, n + 1):
        k1m, k1i = _rhs(m, i, r, phi, eta, delta)
        k2m, k2i = _rhs(m + 0.5 * dt * k1m, i + 0.5 * dt * k1i, r, phi, eta, delta)
        k3m, k3i = _rhs(m + 0.5 * dt * k2m, i + 0.5 * dt * k2i, r, phi, eta, delta)
        k4m, k4i = _rhs(m + dt * k3m, i + dt * k3i, r, phi, eta, delta)
        m += dt * (k1m + 2.0 * k2m + 2.0 * k3m + k4m) / 6.0
        i += dt * (k1i + 2.0 * k2i + 2.0 * k3i + k4i) / 6.0
        if not (math.isfinite(m) and math.isfinite(i)) or m <= 0.0 or m > 1.0 + 1e-9 or i < 0.0:
            ok = False
            M[k:], infl[k:] = np.nan, np.nan
            t[k:] = t[k - 1] + dt * np.arange(1, n - k + 2)
            break
        t[k], M[k], infl[k] = k * dt, m, i
    return {"t": t, "M": M, "I": infl, "in_regime": ok}


def simulate(profile: MinimalProfile = UNTREATED_PROFILE, *,
             burn_in: float = BURN_IN_MONTHS, t_end: float = SCORE_MONTHS,
             dt: float = DT, base: dict[str, float] | None = None) -> dict:
    """Untreated burn-in, then the treated window. Time in months from treatment start.

    The burn-in is always run with the untreated parameters, whatever the
    profile: the drug is given to an established patient, not to the initial
    condition. The returned arrays cover the treated window only.
    """
    base = UNTREATED_PARAMS if base is None else base
    pre = integrate(base, Y0, burn_in, dt) if burn_in > 0 else None
    if pre is not None and not pre["in_regime"]:
        raise RuntimeError("untreated burn-in left the model's regime")
    y0 = (float(pre["M"][-1]), float(pre["I"][-1])) if pre is not None else Y0
    p = params_for(profile, base)
    traj = integrate(p, y0, t_end, dt)
    traj["params"] = p
    traj["profile"] = profile.as_dict()
    traj["in_range"] = {k: bool(RANGES[k][0] < p[k] <= RANGES[k][1]) for k in RATES}
    traj["integrated_I"] = float(np.trapezoid(traj["I"], traj["t"])) if traj["in_regime"] else float("nan")
    traj["mean_M"] = float(np.mean(traj["M"])) if traj["in_regime"] else float("nan")
    return traj


# --------------------------------------------------------------------------- #
# Readouts.
# --------------------------------------------------------------------------- #
def baseline_from(untreated_traj: dict) -> float:
    """The untreated arm's median I: the reference every arm's peaks are counted against."""
    return float(np.median(untreated_traj["I"]))


def relapses_per_year(traj: dict, baseline: float) -> float:
    """Peaks of I above 2x the untreated median, >= 30 days apart, per year.

    `qsp_velez.relapse_events` with its constants untouched; it reads the key
    "E" as the variable to detect peaks on and "t" in days.
    """
    if not traj["in_regime"]:
        return float("nan")
    peaks = relapse_events({"t": traj["t"] * DAYS_PER_MONTH, "E": traj["I"]}, baseline=baseline)
    years = float(traj["t"][-1]) / 12.0
    return len(peaks) / years if years > 0 else float("nan")


def oscillation_period(traj: dict, tail_fraction: float = 0.5) -> float:
    """Mean spacing of I maxima over the last `tail_fraction` of the run; nan if < 2 maxima."""
    t, infl = traj["t"], traj["I"]
    start = int(len(t) * (1.0 - tail_fraction))
    seg = infl[start:]
    idx = [k for k in range(1, len(seg) - 1) if seg[k] >= seg[k - 1] and seg[k] > seg[k + 1]]
    if len(idx) < 2:
        return float("nan")
    times = t[start:][idx]
    return float(np.mean(np.diff(times)))


def oscillation_amplitude(traj: dict, tail_fraction: float = 0.25) -> float:
    """max - min of I over the tail of the run. ~0 at a stable equilibrium."""
    start = int(len(traj["t"]) * (1.0 - tail_fraction))
    seg = traj["I"][start:]
    return float(np.max(seg) - np.min(seg))


def compare(profile: MinimalProfile, untreated: dict | None = None) -> dict:
    """Lesion ratio and relapse proxy of one arm against the untreated arm."""
    untreated = simulate(UNTREATED_PROFILE) if untreated is None else untreated
    treated = simulate(profile)
    base = baseline_from(untreated)
    ratio = (treated["integrated_I"] / untreated["integrated_I"]
             if treated["in_regime"] and untreated["integrated_I"] > 0 else float("nan"))
    return {
        "arm": profile.label,
        "in_regime": treated["in_regime"],
        "in_range": treated["in_range"],
        "lesion_ratio": ratio,
        "relapses_per_year": relapses_per_year(treated, base),
        "untreated_relapses_per_year": relapses_per_year(untreated, base),
        "mean_M": treated["mean_M"],
        "untreated_mean_M": untreated["mean_M"],
        "params": treated["params"],
    }


if __name__ == "__main__":
    u = simulate(UNTREATED_PROFILE)
    print("untreated (figure 2D):", u["params"])
    print(f"  period {oscillation_period(u):.1f} months (paper: {FIG2D_PERIOD_MONTHS:g})")
    print(f"  relapses/yr {relapses_per_year(u, baseline_from(u)):.2f}   mean M {u['mean_M']:.3f}")
    for name, mult in (("phi x0.5", 0.5), ("phi x1.5", 1.5)):
        c = compare(MinimalProfile(label=name, phi=mult), u)
        print(f"  {name}: lesion ratio {c['lesion_ratio']:.3f}  relapses/yr {c['relapses_per_year']:.2f}")
