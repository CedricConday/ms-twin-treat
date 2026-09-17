# CHARGE — ground the intervention class strengths

**Brick:** `bricks/grounding.py`
**Priority:** HIGH — `GROUNDING.md` row 1. This is the one that unlocks the magnitude column.
**Status:** DONE 2026-09-06 — see `NOTES` in this directory and the dated section in
`GROUNDING.md`. Derived value 0.78; direction gate held 4/4; the magnitude gap WIDENED
(-47% -> -77%) and was written up rather than tuned.

---

## The task

`SUPPRESSIVE_STRENGTH = 0.5` and `IMMUNOGENIC_STRENGTH = 0.4` are the last two
hand-set numbers in the grounded path. They are class-level, never per-drug, and
that was the right first move — but they were set from mechanism *reasoning*, not
from data. Replace at least `SUPPRESSIVE_STRENGTH` with a value derived from an
independent measurement.

The intended source, per `GROUNDING.md`: **the IFN-beta effect magnitude in the Kang
2018 data**, which is already loaded, already scored, and has nothing to do with
relapse rates.

---

## THE TRAP — read this before touching a number

The stack currently predicts **~-47%** relapse reduction where the trials report
**-27..-33%**. That gap is the reason this task exists.

**Do not tune the strengths to close it.**

The clinical gate exists to test whether parameters set *without* the answers can
reproduce the answers. Any value chosen because it moves -47% toward -30% has been
fitted to the exam. If the honest derivation lands on -47%, **the deliverable is
-47% plus the write-up of why** — that is a result, not a failure.

From `GROUNDING.md`, and it governs here:

> Ground every parameter against an INDEPENDENT source. Never fit it to the
> clinical outcome the gate is trying to predict.

---

## Deliverable

1. A derivation of `SUPPRESSIVE_STRENGTH` from the Kang IFN-beta magnitude, with the
   arithmetic written down — enough that a reader can redo it without the code.
2. The number changed in `bricks/grounding.py`, with the source cited in the comment
   the way `barrier.py` cites Pardridge 2019.
3. A test in the citation-guard style asserting the derived value, so a silent drift
   fails the suite.
4. Re-run `python -m backtest.clinical`. Record BOTH the direction gate (currently
   4/4) and the magnitude column (currently 0/2) before and after.
5. Append the outcome to `GROUNDING.md` under a dated heading — including the case
   where the gap does not close.

## Acceptance

- Direction gate stays **4/4**. If grounding breaks it, that is reportable, not
  fixable-by-tuning — write it up and stop.
- `IMMUNOGENIC_STRENGTH` stays untouched unless an independent source for it is
  found. Only one arm (APL CGP77116) exercises it; nothing in Kang speaks to it.
  Say so plainly rather than inventing a derivation.
- No change anywhere in `backtest/`. The ruler does not move to fit the result.

## Do not

- Do not fit per-drug. The two-constant class rule is the honesty property; keeping
  it is worth more than a better magnitude number.
- Do not touch the ABM lattice or profiles. `test_profiles_all_use_the_published_lattice`
  is there for a reason.
- Do not report the pipeline as validated. It is not, and a test enforces that.

## Context to read first

- `GROUNDING.md` — the whole file, especially "What grounding the ABM immediately exposed"
- `bricks/grounding.py` — the module docstring states the honest bound already
- `backtest/clinical.py` — how the four arms are scored
- `data/kang.py` — where the IFN-beta magnitude comes from
