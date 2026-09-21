"""Generate the candidate space, screen it, and publish the survivors.

This is the half of `screen/` that produces an artifact. `kill_filter.py`
decides; this enumerates, runs it over the whole space, and writes the result
somewhere a person can read and disagree with it.

WHAT IT REFUSES TO DO, AND HOW THAT IS ENFORCED
-------------------------------------------------
It does not order survivors by predicted benefit. `backtest/lomo.py` measures
whether this model generalises to a mechanism it has never seen and returns
45.9pp MAE against a 12.3pp predict-the-mean null, so a predicted effect size
for a NEW mechanism is not information — and every survivor here is, by
construction, a mechanism pattern no existing drug occupies.

So survivors are listed in **alphabetical order by label**, stated in the output
itself, and `rank_candidates()` in `kill_filter.py` still raises. Sorting by
`best_damage` would be ranking with extra steps, and the sorted list is the
artifact somebody would screenshot.

`best_damage` IS reported per survivor, because hiding it would be worse: it is
the number a reader needs to challenge a verdict. It is labelled as what it is —
a simulated median over eight infection histories, under a model whose
out-of-sample error exceeds its null.

SURVIVING IS NOT PASSING
------------------------
A survivor has only avoided the failures this model can see, and the model is
blind to the entire depleting / sequestering / trafficking class (BUILD_PLAN
§8.4). Read the caveats block in the output before quoting any line of it.

Run:  PYTHONPATH=. python -m screen.report
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

from bricks.profiles import PROFILES, touched_points
from bricks.qsp_velez import INTERVENTION_POINTS
from screen.kill_filter import (
    NOISE_FLOOR,
    REACHABILITY_PROBES,
    SCREEN_SEEDS,
    enumerate_candidates,
    screen,
)

ROOT = Path(__file__).resolve().parent.parent
JSON_OUT = ROOT / "results" / "screen.json"
MD_OUT = ROOT / "docs" / "SCREEN_RESULTS.md"

# The largest relapse reduction any real arm in docs/TRIAL_ANCHORS.md reports.
# Survivors claiming to beat it are the report's most misleading number, so the
# comparison is COMPUTED into the output rather than left to the reader.
BEST_REAL_ARM_PCT = -68.0        # natalizumab, AFFIRM, PMID 16510744

CAVEATS = [
    "Nothing here is evidence about multiple sclerosis. validated=False throughout.",
    "Survivors are listed ALPHABETICALLY, never by predicted benefit. Ranking is "
    "gated on backtest/lomo.py beating its null; it does not (45.9pp vs 12.3pp).",
    "The model is blind to the depleting / sequestering / trafficking class — "
    "natalizumab, fingolimod, ponesimod, alemtuzumab and atacicept cannot come "
    "out beneficial in it at any potency (BUILD_PLAN §8.4). A candidate whose "
    "novelty is in trafficking is invisible to this screen.",
    "A surviving candidate has only avoided the failures this model can see.",
]


def _implausibility(survivors: list[dict]) -> dict:
    """How many survivors claim to beat the best drug ever tested in MS.

    This is the number most likely to be misread off the survivor table, so it
    is computed and stated rather than left implicit. A model that cannot
    reproduce a -30% trial effect (see backtest/clinical.py) reporting -99% for
    an invented dial pair is telling you about itself, not about the candidate.
    """
    rel = sorted((s["best_damage"] / s["untreated_damage"] - 1) * 100
                 for s in survivors
                 if s.get("best_damage") and s.get("untreated_damage"))
    if not rel:
        return {}
    return {
        "n": len(rel),
        "min_pct": rel[0],
        "max_pct": rel[-1],
        "beating_best_real_arm": sum(x < BEST_REAL_ARM_PCT for x in rel),
        "beyond_90pct": sum(x < -90.0 for x in rel),
        "best_real_arm_pct": BEST_REAL_ARM_PCT,
    }


def run(max_points: int = 2, potency: float = 0.5,
        carrying_capacity: float | None = None) -> dict:
    candidates = enumerate_candidates(max_points=max_points, potency=potency)
    results = screen(candidates, carrying_capacity=carrying_capacity)

    by_reason: dict[str, int] = {}
    for r in results:
        key = r.killed_by.name if r.killed_by else "SURVIVED"
        by_reason[key] = by_reason.get(key, 0) + 1

    def row(r) -> dict:
        return {
            "label": r.profile.label,
            "verdict": r.killed_by.name if r.killed_by else "SURVIVED",
            "reason": r.killed_by.value if r.killed_by else r.detail,
            "detail": r.detail,
            "median_damage": r.median_damage,
            "best_damage": r.best_damage,
            "best_potency": r.best_potency,
            "untreated_damage": r.untreated_damage,
            "like_existing": list(r.like_existing),
        }

    survivors = sorted((r for r in results if r.survived),
                       key=lambda r: r.profile.label)
    killed = sorted((r for r in results if not r.survived),
                    key=lambda r: (r.killed_by.name, r.profile.label))

    # The DEGENERATE filter compares candidates against the dial patterns real
    # arms occupy, so this report is only valid for the arm set it ran against.
    # Recording that set is what makes staleness detectable instead of silent:
    # adding one arm on a previously free dial flips candidates to DEGENERATE.
    occupied = sorted({"|".join(touched_points(prof))
                       for name, prof in PROFILES.items() if name != "untreated"})

    imp = _implausibility([row(r) for r in survivors])
    return {
        "occupied_patterns": occupied,
        "n_arms_in_library": len(PROFILES) - 1,
        "implausibility": imp,
        "generated": date.today().isoformat(),
        "carrying_capacity": carrying_capacity,
        "n_candidates": len(candidates),
        "max_points": max_points,
        "enumerated_potency": potency,
        "intervention_points": list(INTERVENTION_POINTS),
        "screen_seeds": list(SCREEN_SEEDS),
        "reachability_probes": list(REACHABILITY_PROBES),
        "noise_floor": NOISE_FLOOR,
        "counts": by_reason,
        "survivors": [row(r) for r in survivors],
        "killed": [row(r) for r in killed],
        "ordering": "survivors alphabetical by label; NOT by predicted benefit",
        "caveats": CAVEATS,
        "validated": False,
    }


def to_markdown(rep: dict) -> str:
    k = rep.get("carrying_capacity")
    banner = ([] if k is None else [
        f"> **RUN UNDER THE CAPACITY EXTENSION (K = {k:g}), NOT THE TRANSCRIPTION.**",
        "> `bricks/qsp_velez.py`'s carrying capacity is an extension to Vélez de",
        "> Mendizábal 2011, off by default. These results describe a DIFFERENT MODEL",
        "> from `SCREEN_RESULTS.md`. The extension does not lift the ranking gate",
        "> (`backtest/lomo_capacity.py`: 45.6pp against a 12.3pp null); it changes",
        "> which candidates are screenable at all.",
        "",
    ])
    lines = [
        "# Screen results — candidates the model cannot rule out",
        "",
        *banner,
        f"Generated {rep['generated']} by `PYTHONPATH=. python -m screen.report`. "
        f"{rep['n_candidates']} candidates: every combination of up to "
        f"{rep['max_points']} of the model's {len(rep['intervention_points'])} named "
        f"intervention points, each turned both ways at potency "
        f"{rep['enumerated_potency']}.",
        "",
        "**This is a kill list with a remainder, not a ranking.** " + CAVEATS[1],
        "",
        "## Counts",
        "",
        "| verdict | n |",
        "|---|---|",
    ]
    for k, v in sorted(rep["counts"].items()):
        lines.append(f"| {k} | {v} |")

    lines += ["", f"## Survivors ({len(rep['survivors'])})", "",
              "Alphabetical. The damage column is a simulated median over "
              f"{len(rep['screen_seeds'])} infection histories, shown so a reader can "
              "challenge the verdict — not as a score.", "",
              "| candidate | best simulated damage | vs untreated | at potency |",
              "|---|---|---|---|"]
    for s in rep["survivors"]:
        unt = s["untreated_damage"]
        rel = f"{(s['best_damage'] / unt - 1) * 100:+.0f}%" if unt and s["best_damage"] else "—"
        pot = s["best_potency"] if s["best_potency"] is not None else "as enumerated"
        bd = f"{s['best_damage']:.4f}" if s["best_damage"] is not None else "—"
        lines.append(f"| `{s['label']}` | {bd} | {rel} | {pot} |")

    lines += ["", "## Killed", "", "| candidate | verdict | why |", "|---|---|---|"]
    for k in rep["killed"]:
        lines.append(f"| `{k['label']}` | {k['verdict']} | {k['detail']} |")

    lines += ["", "## Read this before quoting any line above", ""]
    imp = rep.get("implausibility") or {}
    if imp:
        lines += [
            f"- **{imp['beating_best_real_arm']} of {imp['n']} survivors claim a larger "
            f"effect than the best drug ever tested in MS** (natalizumab, "
            f"{imp['best_real_arm_pct']:.0f}% relapse reduction in AFFIRM), and "
            f"{imp['beyond_90pct']} of them claim better than -90%. Survivor effects "
            f"here run {imp['min_pct']:.0f}% to {imp['max_pct']:.0f}%. **This is not "
            "credible and it is not meant to be read as a prediction.** The same model "
            "cannot reproduce a -30% effect for interferon beta on a real arm. The "
            "damage column measures the model, not the candidate.",
        ]
    lines += [f"- {c}" for c in CAVEATS]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--carrying-capacity", type=float, default=None,
                    help="run the screen under the qsp_velez EXTENSION instead of "
                         "the transcription; writes to a separate, labelled file")
    args = ap.parse_args(argv)
    k = args.carrying_capacity

    rep = run(carrying_capacity=k)
    json_out = JSON_OUT if k is None else JSON_OUT.with_name(f"screen_K{k:g}.json")
    md_out = MD_OUT if k is None else MD_OUT.with_name(f"SCREEN_RESULTS_K{k:g}.md")
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(rep, indent=2))
    md_out.parent.mkdir(parents=True, exist_ok=True)
    md_out.write_text(to_markdown(rep))

    print(f"{rep['n_candidates']} candidates screened\n")
    for k, v in sorted(rep["counts"].items()):
        print(f"  {k:<16} {v:>4}")
    print(f"\n  survivors, alphabetical (NOT ranked — {CAVEATS[1]}):")
    for s in rep["survivors"]:
        print(f"    {s['label']}")
    imp = rep.get("implausibility") or {}
    if imp:
        print(f"\n  {imp['beating_best_real_arm']} of {imp['n']} survivors claim a "
              f"bigger effect than natalizumab ({imp['best_real_arm_pct']:.0f}%), "
              f"{imp['beyond_90pct']} claim better than -90%. Not credible; the damage "
              "column measures the model, not the candidate.")
    print(f"\nwritten: {json_out}\n         {md_out}")
    for c in CAVEATS:
        print(f"  - {c}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
