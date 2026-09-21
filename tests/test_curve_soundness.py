"""The cached response curve is per-MECHANISM, and one property has to hold for
that to be sound.

`backtest/lomo.build_table` tabulates one column per mechanism pattern and
computes it from a single representative arm — `groups[pattern][0]`. Every other
arm sharing that pattern then reads the representative's curve. That is the
saving which makes the table affordable: five gamma_E arms are one column, not
five.

It is sound only while the arms sharing a pattern move their dials in the same
DIRECTION, and nothing enforces that. `bricks/profiles.touched_points` returns
point NAMES only, so an arm that raised `alpha_E` and an arm that lowered it
would land in the same group, and whichever sorted first would supply the curve
for both.

Today every group is homogeneous, so this is a latent hazard rather than a bug.
It is worth a test rather than a comment because of what it would look like if it
ever fired: no crash, no missing key, just an arm predicted with the wrong sign
of its own mechanism — which is exactly the class of failure this repo spent a
night discovering it cannot see.

WHY THIS IS THE RIGHT GUARD AND AN ARM-SET FINGERPRINT IS NOT
--------------------------------------------------------------
`results/lomo_certificate.json` needed a fingerprint because it is a measurement
OF AN EXAM: change the arms and the number stops describing anything current.
`results/mechanism_curve*.json` is not that. It is a measurement of the MODEL —
pattern by potency — and wiring an arm onto an existing dial does not invalidate
it, because the curve was never per-arm. Wiring an arm onto a NEW dial fails
loudly already, with a KeyError on the lookup.

The one case that is silently wrong is a direction collision, and a fingerprint
would not catch it: it would report "the arm set changed", prompt an eight-minute
rebuild, and the rebuild would pool the same two arms the same wrong way. The
fingerprint would make someone pay for the answer without giving it to them.
"""

from __future__ import annotations

from backtest.lomo import mechanism_groups
from bricks.profiles import PROFILES, touched_points


def _direction(arm: str, point: str) -> str:
    return "down" if getattr(PROFILES[arm], point) < 1.0 else "up"


def test_every_mechanism_group_moves_its_dials_the_same_way():
    """If this fails, the cached curve is pooling arms it must not pool, and the
    arms on the losing side of the collision are being predicted with the wrong
    sign. Split the group by direction before rebuilding anything."""
    offenders = []
    for pattern, arms in mechanism_groups().items():
        if len(arms) < 2:
            continue
        for point in pattern:
            directions = {_direction(a, point) for a in arms}
            if len(directions) > 1:
                offenders.append(
                    f"{'|'.join(pattern)} on {point}: "
                    + ", ".join(f"{a}={_direction(a, point)}" for a in arms))
    assert not offenders, (
        "mechanism groups mix directions, so the shared response curve is unsound:\n  "
        + "\n  ".join(offenders))


def test_the_representative_is_deterministic():
    """build_table uses groups[pattern][0]. If group membership order depended on
    anything unstable, the curve's identity would change between rebuilds without
    the arm set changing."""
    first = {p: arms[0] for p, arms in mechanism_groups().items()}
    second = {p: arms[0] for p, arms in mechanism_groups().items()}
    assert first == second


def test_every_group_representative_is_a_member_of_its_group():
    for pattern, arms in mechanism_groups().items():
        assert arms[0] in arms
        assert touched_points(PROFILES[arms[0]]) == pattern
