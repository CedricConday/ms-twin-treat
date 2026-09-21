# Research quality and limitations

This repository is research software. It is not a clinical system, not a
regulatory artefact, and not evidence about multiple sclerosis.

Its central discipline is that **software correctness and biological validity
are different claims**, and the repo is required to say which one it is making.

## Evidence labels

The codebase already carries these at runtime; this is what they mean.

| Label | Where it appears | What it claims |
|---|---|---|
| **Data-grounded** | `GROUNDING.md` status table | A parameter is derived from an independent measured source, with the derivation in the repo (e.g. `scripts/derive_suppressive_strength.py`). |
| **Validated (software)** | test suite, `backtest/selftest.py` | The code behaves as specified. Nothing more. |
| **Illustrative** | `validated=False` on a brick's output | A mechanistic toy. Directionally reasoned, parameters not fitted to patient data. |
| **Stand-in** | `STANDIN` key in the pipeline state | An explicit placeholder that exists to keep the scales composable. |

A passing test demonstrates that code behaves as specified. It does not
demonstrate that a model is biologically correct, clinically predictive, safe,
or therapeutically useful. `spine/run_demo.py` prints this distinction at the
end of every run on purpose.

## The one rule that keeps grounding honest

From `GROUNDING.md`, restated here because it is the load-bearing one:

> Ground every parameter against an **independent** source. Never fit it to the
> clinical outcome the gate is trying to predict.

Fitting a parameter to the outcome makes the clinical gate circular. A gate that
restates what it was told is not a test.

## Review checklist

Before treating a result as informative, record:

1. dataset and exact source identifier;
2. preprocessing and filtering choices;
3. random seeds, package pins, and the commit SHA;
4. whether the test was prospective, held out, or used to design the model;
5. the comparison against the identity and mean-shift nulls, both reported;
6. the measurement-noise ceiling, and which folds were excluded as unreliable;
7. which outputs are illustrative or stand-ins;
8. uncertainty and failure cases;
9. **which exam a cached number was measured on** — see below.

Point 6 is not optional here. The Kang backtest carries a per-fold reliability
number precisely because one cell type (Megakaryocytes, 63/69 cells) scores
noise for every model including the nulls, and silently drags every aggregate.

## Staleness is the failure mode this repository actually has

Most of the discipline above guards against claiming too much from a measurement.
The failure that has actually occurred here, repeatedly, is subtler: a
measurement that was correct when taken, still being read after the thing it
measured changed.

On one night — 2026-09-20 into 2026-09-21 — four instances were found in two
independent work streams:

- a cached mechanism-holdout run measured over twelve quantified arms, still
  being read by the decision gate after fifteen were wired;
- a screen's survivor list computed against the dials the arms occupied *before*
  three more arms were wired onto them, which changes which candidates count as
  duplicates of an existing drug;
- a cached capacity comparison with no record of the arm set it ran against;
- a hardcoded summary string in `backtest/clinical_velez.py` printing a
  measurement from a previous arm set beside a freshly computed one.

- and, in the same night's clean-up, **five stale figures that were fixed,
  committed, and reported as fixed — while sitting unpushed.** The working tree
  was clean and the fixes were real; they were simply not in the repository. The
  other work stream, scanning a fresh clone of the remote, reported the same
  figures still stale, which from where it stood was indistinguishable from the
  fix having been wrong.

Plus a citation that resolved to the wrong paper (one digit out, with the row's
own numbers correct), and a model comparison quoted as a separation when the gap
sat inside the stated noise floor.

That last one has the same shape as all the caches and is worth naming as its
own rule: **committed is not published.** A cache is correct on disk and quoted
against a different exam; a fix is correct in a working tree and quoted as being
in the repo. In both cases an artifact that is correct *somewhere* is reported as
correct *everywhere*, and nothing local can tell the difference — a clean `git
status` looks identical whether or not you have pushed.

**None of these were caught by a failing test.** Every one was found by a person
reading, and mechanised afterwards. The unpushed-fixes case was found by reading
from a *clean clone of the remote* rather than from the tree the work was done
in — and that is the general point, not a fact about clones. **Reading your own
work from the place you did it cannot show you what you left behind.** It applies
to a cache you did not regenerate, to a commit you did not push, and to a
document you wrote and did not re-read. That ordering is the honest lesson, and it is
the opposite of the reassuring one: the checks in this repo do not find staleness,
they *keep it found* once a human has found it once. `scripts/verify_gate_docs.py`,
`scripts/verify_anchors.py` and `scripts/state_of_build.py` all exist because
someone read something and noticed, not because a test went red.

### The rule that came out of it

**A cached measurement carries the exam it was taken on, and is refused when that
exam has changed.** Not the date, not the commit — those were present on every
stale artifact above and told nobody anything. The arm set itself, by name rather
than by count, because swapping one arm for another leaves the count unchanged
and changes the exam completely.

Caches under `results/` that this applies to are listed in
`docs/REPRODUCIBILITY.md`, with what regenerates each.

### The limit of that rule, which matters as much as the rule

**A staleness guard tells you an input changed. It does not tell you a result is
wrong, and it does not tell you a result is right.** Those are three different
claims and the guard only makes the first.

The case that forced this distinction: `results/mechanism_curve*.json` was
nominated for the same arm-set fingerprint, and it does not need one. It is a
measurement of the *model* — one column per mechanism pattern per potency — not
of the exam, so wiring an arm onto an existing dial leaves it exactly valid, and
wiring one onto a new dial already fails loudly with a missing key. A fingerprint
there would fire on changes that invalidate nothing, charge an eight-minute
rebuild, and hand back a byte-identical table.

Worse, it would not catch that file's one genuinely silent failure. Each pattern's
column is computed from a single representative arm and shared by every arm on
that pattern, which is sound only while those arms move their dials in the same
direction — a property nothing enforced, because the grouping key carries point
names and not signs. An arm wired onto an existing dial pointing the other way
*does* change the arm set, so the fingerprint fires, and the rebuild pools the
same two arms the same wrong way. Loud, expensive, and still wrong.

`tests/test_curve_soundness.py` enforces the property instead. If it ever fires,
**split the group by direction before rebuilding anything** — the other cache
failures in this repo all taught the reflex of reaching for a regeneration, and
here a regeneration reproduces the fault.

The general lesson is that paying a cost creates the feeling of having checked
something. A guard that fires on the wrong signal is not neutral; it buys false
confidence at real expense.

### What this means when reading a number in this repo

A figure quoted in prose is the least trustworthy form a number takes here. It
was true of some exam, and the exam is not usually stated. Prefer a number that
a named command regenerates, and when one is quoted, check what it was measured
on before relying on it.

## What the clinical gate is and is not

The clinical gate is a **capability milestone**. Every arm in it is one whose
outcome informed the setup, so direction-correctness is necessary and not
sufficient. Historical outcomes used to design or calibrate a test cannot also
be counted as an untouched prospective confirmation of it.

Real viability needs data-grounded parameters, out-of-sample arms, and validated
magnitudes. `GROUNDING.md` tracks how far that has actually got.
