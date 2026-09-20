"""Tests for the DREAM(ZS) port and the posterior virtual population.

A sampler that runs is not a sampler that samples, and the failure is silent:
badly-mixed chains return numbers that look exactly like good ones. So the first
group checks it recovers distributions whose answers are known, and the last
group checks that the vpop wrapper REFUSES to hand back a cohort from chains
that did not converge — the first version of it returned one at R-hat 1.34 and
nothing in the output said so.
"""

from __future__ import annotations

import numpy as np
import pytest

from bricks import dream, vpop

# --------------------------------------------------------------------------- #
# constants, each sourced. See the module docstring in bricks/dream.py.
# --------------------------------------------------------------------------- #

def test_constants_match_a_working_implementation():
    """Read off the `bumps` DREAM core rather than recalled.

    Getting a sampler's constants subtly wrong produces one that runs, mixes
    badly, and lies quietly — so they are pinned.
    """
    assert dream.DE_EPS == 0.05          # bumps DE_eps
    assert dream.DE_NOISE == 1e-6        # bumps DE_noise
    assert dream.DE_PAIRS == 3           # bumps DE_pairs
    assert dream.SNOOKER_RATE == 0.1     # bumps DE_snooker_rate
    assert dream.SNOOKER_GAMMA == (1.2, 2.2)   # ter Braak & Vrugt 2008


# --------------------------------------------------------------------------- #
# does it actually sample?
# --------------------------------------------------------------------------- #

def test_recovers_a_known_gaussian():
    mu, sd = np.array([1.4, 0.6]), np.array([0.15, 0.10])

    def logp(theta):
        return float(-0.5 * np.sum(((theta - mu) / sd) ** 2))

    res = dream.sample(logp, [(1.0, 2.0), (0.25, 2.0)], n_chains=5,
                       n_steps=1200, seed=1)
    post = res.flat()
    assert np.allclose(post.mean(axis=0), mu, atol=0.05)
    assert np.allclose(post.std(axis=0), sd, rtol=0.35)
    assert np.all(res.gelman_rubin() < 1.2)


def test_recovers_a_uniform_box():
    """A flat density must give a flat posterior, not pile up in the middle."""
    res = dream.sample(lambda t: 0.0, [(0.0, 1.0), (0.0, 1.0)],
                       n_chains=5, n_steps=1500, seed=3)
    post = res.flat()
    assert np.allclose(post.mean(axis=0), 0.5, atol=0.06)
    # A uniform on [0,1] has sd 1/sqrt(12) ~ 0.289.
    assert np.allclose(post.std(axis=0), 0.289, rtol=0.25)


def test_acceptance_rate_is_in_a_workable_band():
    mu = np.array([0.5, 0.5])

    def logp(theta):
        return float(-0.5 * np.sum(((theta - mu) / 0.2) ** 2))

    res = dream.sample(logp, [(0.0, 1.0), (0.0, 1.0)], n_chains=5,
                       n_steps=600, seed=5)
    assert 0.1 < res.acceptance_rate < 0.7, (
        f"acceptance {res.acceptance_rate:.3f} — near 0 means the proposal is too "
        "wide, near 1 means it is too narrow; either way the chain is not exploring")


def test_samples_stay_inside_the_bounds():
    res = dream.sample(lambda t: 0.0, [(0.3, 0.7), (1.0, 2.0)],
                       n_chains=4, n_steps=300, seed=7)
    post = res.flat()
    assert np.all(post[:, 0] >= 0.3) and np.all(post[:, 0] <= 0.7)
    assert np.all(post[:, 1] >= 1.0) and np.all(post[:, 1] <= 2.0)


def test_gelman_rubin_flags_chains_that_have_not_mixed():
    """A bimodal target with narrow wells strands chains; R-hat must say so."""
    def logp(theta):
        a = -0.5 * np.sum(((theta - np.array([0.1, 0.1])) / 0.01) ** 2)
        b = -0.5 * np.sum(((theta - np.array([0.9, 0.9])) / 0.01) ** 2)
        return float(np.logaddexp(a, b))

    res = dream.sample(logp, [(0.0, 1.0), (0.0, 1.0)], n_chains=4,
                       n_steps=200, seed=11)
    assert np.any(res.gelman_rubin() > 1.2) or res.acceptance_rate < 0.05


def test_reproducible_for_a_given_seed():
    f = lambda t: float(-np.sum(t ** 2))  # noqa: E731
    a = dream.sample(f, [(-1.0, 1.0)] * 2, n_chains=3, n_steps=120, seed=2)
    b = dream.sample(f, [(-1.0, 1.0)] * 2, n_chains=3, n_steps=120, seed=2)
    assert np.array_equal(a.chains, b.chains)


def test_too_few_chains_is_refused():
    with pytest.raises(ValueError, match="at least 3 chains"):
        dream.sample(lambda t: 0.0, [(0.0, 1.0)], n_chains=2, n_steps=10)


def test_an_everywhere_impossible_density_raises():
    with pytest.raises(RuntimeError, match="-inf almost everywhere"):
        dream.sample(lambda t: -np.inf, [(0.0, 1.0)] * 2, n_chains=3, n_steps=10)


# --------------------------------------------------------------------------- #
# the plausibility density
# --------------------------------------------------------------------------- #

def test_plausibility_shoulder_is_in_log_damage_space():
    """The bug this was first written with.

    With a shoulder of 0.5 in LINEAR damage units, a healthy configuration at
    damage 0.0086 against a lower edge of 0.05 is 0.041 away — a penalty of
    -0.003, i.e. none — and the posterior filled with configurations the hard
    filter rejects. Damage spans ~0.0086 to 1000 here, so the tolerance has to
    be multiplicative.
    """
    inside = vpop.plausibility_log_density(np.array([1.5, 0.3]))
    healthy = vpop.plausibility_log_density(np.array([1.5, 1.9]))
    assert inside == 0.0
    assert healthy < -4.0, (
        f"a healthy configuration should be strongly disfavoured, got {healthy:.3f}")


def test_plausibility_is_monotone_in_alpha_R():
    vals = [vpop.plausibility_log_density(np.array([1.5, a]))
            for a in (0.3, 0.6, 1.0, 1.5, 1.9)]
    assert vals == sorted(vals, reverse=True)


def test_density_uses_fixed_seeds_so_it_is_deterministic():
    """MCMC on a density that resamples its noise is pseudo-marginal MCMC with an
    unstated error term. Common random numbers are what make this valid."""
    theta = np.array([1.4, 0.4])
    assert vpop.plausibility_log_density(theta) == vpop.plausibility_log_density(theta)
    assert len(vpop.DENSITY_SEEDS) >= 4


# --------------------------------------------------------------------------- #
# the convergence guard
# --------------------------------------------------------------------------- #

def test_cohort_is_refused_when_chains_have_not_converged():
    """The guard that matters. An unconverged cohort is a walk, and it is
    indistinguishable from a good one by inspection."""
    with pytest.raises(RuntimeError, match="have not converged"):
        vpop.sample_vpop_dream(n=4, seed=1, n_chains=3, n_steps=12)


def test_the_guard_can_be_opted_out_of_explicitly():
    cohort, res = vpop.sample_vpop_dream(n=3, seed=1, n_chains=3, n_steps=14,
                                         require_convergence=False)
    assert len(cohort) == 3
    assert "r_hat" in cohort[0]["vpop_meta"]


def test_r_hat_threshold_is_the_conventional_one():
    assert vpop.R_HAT_MAX == 1.2
