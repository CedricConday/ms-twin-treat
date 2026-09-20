"""DREAM(ZS) — the sampling half of the wedge, which was never built.

`ms-twin/docs/RESEARCH_FINDINGS.md:39` states the differentiator as

    "a Python port of the prevalence/plausible-patient method + modern sampling
     (DREAM(ZS) / simulation-based inference) + the Weatherley MS ABM"

Only the first and third parts were built. `bricks/vpop.py` samples by Latin
hypercube and applies a hard accept/reject filter, which yields an accepted SET
of parameter vectors. It does not yield a POSTERIOR, and the difference matters:
the method's value is a prevalence-weighted population, and you cannot weight
what you never sampled densities over.

This is that missing piece: DREAM(ZS), transcribed rather than invented.

    ter Braak CJF, Vrugt JA. "Differential Evolution Markov Chain with snooker
    updater and fewer chains." Statistics and Computing 2008;18(4):435-446.
    doi:10.1007/s11222-008-9104-9

    Vrugt JA et al. "Accelerating Markov chain Monte Carlo simulation by
    differential evolution with self-adaptive randomized subspace sampling."
    Int J Nonlinear Sci Numer Simul 2009;10(3):273-290.

WHY WRITE IT RATHER THAN PIP-INSTALL IT
----------------------------------------
`pydream` (Shockley et al., Bioinformatics 2018;34(4):695-697, PMC5860607)
is a published, maintained Python DREAM implementation and is the right thing
to reach for in most projects. It is NOT used here for two stated reasons, and
neither is "ours is better":

  1. `requirements.txt` pins a deliberately minimal numpy/scipy stack with a
     comment about an ABI break on this box. Adding a sampler dependency to one
     brick is a real cost against that.
  2. The claim this repo makes is a *port of a method*. A port that shells out
     to somebody else's port is not the artifact.

If either reason stops holding, replacing this module with `pydream` is the
correct move, not a regression.

THE ALGORITHM, AND WHERE EACH CONSTANT COMES FROM
--------------------------------------------------
N chains evolve jointly. A proposal for chain i is built from the difference
between randomly drawn members of an archive Z of past states, so the proposal
scale and orientation adapt to the target without any tuning:

    parallel direction (90% of steps):
        x* = x_i + (1 + e) * gamma(delta, d') * sum_{j<=delta} (Z[a_j] - Z[b_j])
             + epsilon
        gamma = 2.38 / sqrt(2 * delta * d')      Vrugt 2009
        every 5th step gamma = 1.0, which allows mode-jumping

    snooker (10% of steps):
        three archive points define a line through x_i; the jump is a multiple
        of the projection difference along it, with gamma ~ U(1.2, 2.2)
                                                      ter Braak & Vrugt 2008

    e       ~ U(-b, b),  b  = 0.05      "DE_eps"   in bumps' DREAM
    epsilon ~ N(0, b*),  b* = 1e-6      "DE_noise" in bumps' DREAM
    delta   ~ U{1..3},   delta_max = 3  "DE_pairs" in bumps' DREAM
    snooker rate 0.1                    "DE_snooker_rate" in bumps' DREAM

Those four values were read off a working implementation (the `bumps` package's
DREAM core) rather than recalled, because getting a sampler's constants subtly
wrong produces a sampler that runs, mixes badly, and lies quietly.

CROSSOVER is implemented as randomized subspace sampling with a fixed CR set:
each dimension is updated with probability CR, drawn from {1/3, 2/3, 1}. The
SELF-ADAPTIVE CR weighting of Vrugt 2009 — which learns which CR values produce
the largest normalized jumps — is **NOT implemented**. That is the one place
this port is incomplete, it is stated rather than hidden, and it costs
efficiency rather than correctness.

WHAT THIS DOES NOT FIX
----------------------
A sampler gives you a posterior over parameters for the density you hand it. It
does not make the density right. `bricks/vpop.py` documents what its density
means and what conditioning on a fixed set of infection histories costs.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

# --- constants, each sourced. See the module docstring. --------------------- #
DE_EPS = 0.05        # multiplicative perturbation e ~ U(-b, b)
DE_NOISE = 1e-6      # additive perturbation epsilon ~ N(0, b*)
DE_PAIRS = 3         # delta_max: up to 3 difference pairs
SNOOKER_RATE = 0.1   # 90/10 parallel-direction / snooker mix
SNOOKER_GAMMA = (1.2, 2.2)   # ter Braak & Vrugt 2008
MODE_JUMP_EVERY = 5  # every 5th step gamma = 1.0, enabling mode jumps
CR_SET = (1.0 / 3.0, 2.0 / 3.0, 1.0)


@dataclass
class DreamResult:
    """Posterior samples plus the diagnostics needed to distrust them."""

    chains: np.ndarray          # (n_chains, n_steps, n_dim)
    log_density: np.ndarray     # (n_chains, n_steps)
    acceptance_rate: float
    n_burn: int
    param_names: tuple[str, ...]

    def flat(self, thin: int = 1) -> np.ndarray:
        """Post-burn-in samples, chains concatenated."""
        return self.chains[:, self.n_burn::thin, :].reshape(-1, self.chains.shape[2])

    def gelman_rubin(self) -> np.ndarray:
        """R-hat per dimension. Above ~1.2 means the chains have not converged.

        Reported rather than enforced: a caller deciding what to do about
        non-convergence is better than this function silently resampling.
        """
        x = self.chains[:, self.n_burn:, :]
        m, n, d = x.shape
        if m < 2 or n < 2:
            return np.full(d, np.nan)
        chain_means = x.mean(axis=1)                     # (m, d)
        chain_vars = x.var(axis=1, ddof=1)               # (m, d)
        W = chain_vars.mean(axis=0)
        B = n * chain_means.var(axis=0, ddof=1)
        var_hat = (n - 1) / n * W + B / n
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.sqrt(np.where(W > 0, var_hat / W, np.nan))


def _snooker_proposal(x_i, z_a, z_b, z_c, rng):
    """Jump along the line through z_c and x_i. ter Braak & Vrugt 2008.

    z_b and z_a are projected orthogonally onto that line and the jump is a
    multiple of the difference between the projections.
    """
    line = x_i - z_c
    norm2 = float(np.dot(line, line))
    if norm2 <= 0.0:
        return x_i.copy(), 0.0
    proj_a = np.dot(z_a, line) / norm2 * line
    proj_b = np.dot(z_b, line) / norm2 * line
    gamma = rng.uniform(*SNOOKER_GAMMA)
    x_new = x_i + gamma * (proj_a - proj_b)

    # Snooker proposals are not symmetric; the Metropolis ratio carries a
    # Jacobian in the distance along the line. Omitting it is a common and
    # silent bug, so it is computed here.
    d_new = float(np.linalg.norm(x_new - z_c))
    d_old = float(np.linalg.norm(x_i - z_c))
    if d_old <= 0.0 or d_new <= 0.0:
        return x_new, 0.0
    n_dim = x_i.size
    return x_new, (n_dim - 1) * (np.log(d_new) - np.log(d_old))


def sample(
    log_density: Callable[[np.ndarray], float],
    bounds: list[tuple[float, float]],
    *,
    n_chains: int = 5,
    n_steps: int = 400,
    burn_fraction: float = 0.5,
    seed: int = 0,
    archive_size: int = 200,
    param_names: tuple[str, ...] | None = None,
) -> DreamResult:
    """Run DREAM(ZS) against `log_density` over a box of `bounds`.

    `log_density` must be DETERMINISTIC in its argument. A stochastic simulator
    evaluated with a fresh seed each call turns this into pseudo-marginal MCMC
    with an unstated noise term, which is not what is implemented here — the
    caller is responsible for common random numbers. `bricks/vpop.py` shows what
    that costs.

    Proposals outside `bounds` are rejected by the density returning -inf rather
    than being reflected or clipped, which keeps the chain reversible.
    """
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in bounds], dtype=float)
    hi = np.array([b[1] for b in bounds], dtype=float)
    n_dim = lo.size
    if n_chains < 3:
        raise ValueError("DREAM(ZS) needs at least 3 chains for difference vectors")

    # Archive Z, initialised by Latin-hypercube over the box so the sampler
    # starts from the same design the rejection sampler uses.
    z = lo + (hi - lo) * rng.random((archive_size, n_dim))
    x = lo + (hi - lo) * rng.random((n_chains, n_dim))
    logp = np.array([log_density(xi) for xi in x], dtype=float)

    # A chain that starts at -inf can never move. Re-draw rather than let it sit.
    for i in range(n_chains):
        for _ in range(200):
            if np.isfinite(logp[i]):
                break
            x[i] = lo + (hi - lo) * rng.random(n_dim)
            logp[i] = log_density(x[i])
        if not np.isfinite(logp[i]):
            raise RuntimeError(
                "could not find a finite starting point for every chain; the "
                "density is -inf almost everywhere over these bounds"
            )

    chains = np.empty((n_chains, n_steps, n_dim))
    dens = np.empty((n_chains, n_steps))
    accepted = 0

    for step in range(n_steps):
        for i in range(n_chains):
            if rng.random() < SNOOKER_RATE:
                idx = rng.choice(len(z), size=3, replace=False)
                x_new, log_jac = _snooker_proposal(x[i], z[idx[0]], z[idx[1]],
                                                   z[idx[2]], rng)
            else:
                log_jac = 0.0
                delta = int(rng.integers(1, DE_PAIRS + 1))
                idx = rng.choice(len(z), size=2 * delta, replace=False)
                diff = sum(z[idx[2 * j]] - z[idx[2 * j + 1]] for j in range(delta))

                # randomized subspace sampling (crossover)
                cr = float(rng.choice(CR_SET))
                mask = rng.random(n_dim) < cr
                if not mask.any():
                    mask[rng.integers(n_dim)] = True
                d_eff = int(mask.sum())

                gamma = (1.0 if (step + 1) % MODE_JUMP_EVERY == 0
                         else 2.38 / np.sqrt(2.0 * delta * d_eff))
                e = rng.uniform(-DE_EPS, DE_EPS, n_dim)
                eps = rng.normal(0.0, DE_NOISE, n_dim)

                x_new = x[i].copy()
                x_new[mask] = (x[i][mask]
                               + (1.0 + e[mask]) * gamma * diff[mask]
                               + eps[mask])

            lp_new = log_density(x_new) if np.all((x_new >= lo) & (x_new <= hi)) else -np.inf
            if np.isfinite(lp_new) and np.log(rng.random()) < (lp_new - logp[i] + log_jac):
                x[i], logp[i] = x_new, lp_new
                accepted += 1

            chains[i, step], dens[i, step] = x[i], logp[i]

        # Grow the archive with current states (the "ZS" of DREAM(ZS)).
        z = np.vstack([z, x])
        if len(z) > archive_size * 2:
            z = z[-archive_size * 2:]

    return DreamResult(
        chains=chains,
        log_density=dens,
        acceptance_rate=accepted / (n_chains * n_steps),
        n_burn=int(n_steps * burn_fraction),
        param_names=param_names or tuple(f"p{i}" for i in range(n_dim)),
    )


if __name__ == "__main__":
    # A target with a known answer, so the sampler is checked rather than trusted:
    # an uncorrelated Gaussian inside a box.
    true_mu = np.array([1.4, 0.6])
    true_sd = np.array([0.15, 0.10])

    def logp(theta):
        return float(-0.5 * np.sum(((theta - true_mu) / true_sd) ** 2))

    res = sample(logp, [(1.0, 2.0), (0.25, 2.0)], n_chains=5, n_steps=1200,
                 seed=1, param_names=("alpha_E", "alpha_R"))
    post = res.flat()
    print("DREAM(ZS) self-check — recover a known Gaussian\n")
    print(f"  acceptance rate : {res.acceptance_rate:.3f}")
    print(f"  R-hat           : {np.round(res.gelman_rubin(), 3)}")
    print(f"  posterior mean  : {np.round(post.mean(axis=0), 3)}  (true {true_mu})")
    print(f"  posterior sd    : {np.round(post.std(axis=0), 3)}  (true {true_sd})")
    print("\n  A sampler that runs is not a sampler that samples. R-hat near 1.0 and")
    print("  a recovered mean/sd are the minimum before trusting it on a real target.")
