# Contributing

The thing worth preserving here is the distinction between **software
correctness** and **biological validity**. A contribution that blurs it costs
more than it adds, however good the code is.

Read [`docs/QUALITY.md`](docs/QUALITY.md) first, and
[`GROUNDING.md`](GROUNDING.md) if you are touching a parameter.

## Before opening a pull request

```bash
python -m ruff check .
python -m pytest
```

Both are the gates CI runs. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)
for the environment — the pins matter, and a system Python will give you
different numbers.

## New model components

Document, in the module itself:

- the assumption being made, and where it came from;
- the provenance of every parameter — measured, published, ported, or reasoned;
- whether the output is data-grounded, illustrative, or a stand-in, using the
  labels the pipeline already enforces (`validated=False`, `STANDIN`).

A brick whose output is a toy must say so in its output, not only in a comment.
`spine/pipeline.py` enforces this rather than trusting convention, and that is
deliberate.

## New parameters

Ground against an **independent** source. Never fit a parameter to the clinical
outcome the gate is trying to predict — that makes the gate circular. If no
independent source exists, say so and leave the value labelled ungrounded;
`IMMUNOGENIC_STRENGTH` is the worked example of an honest gap.

Where a constant is derived, the derivation belongs in the repo as a script, and
a test should assert that the constant is still what the derivation produces.
`tests/test_grounding.py` is the pattern.

## New benchmarks

State the data source, the split strategy, the nulls, and any route by which
information could leak from the held-out fold into the model. Report the noise
ceiling alongside the score. A benchmark without a null is not a benchmark.

## What not to do

Do not describe a passing unit test, a backtest on data used during development,
or a directional toy result as clinical evidence. Do not quote a number in the
README that the code does not currently produce.
