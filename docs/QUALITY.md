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
8. uncertainty and failure cases.

Point 6 is not optional here. The Kang backtest carries a per-fold reliability
number precisely because one cell type (Megakaryocytes, 63/69 cells) scores
noise for every model including the nulls, and silently drags every aggregate.

## What the clinical gate is and is not

The clinical gate is a **capability milestone**. Every arm in it is one whose
outcome informed the setup, so direction-correctness is necessary and not
sufficient. Historical outcomes used to design or calibrate a test cannot also
be counted as an untouched prospective confirmation of it.

Real viability needs data-grounded parameters, out-of-sample arms, and validated
magnitudes. `GROUNDING.md` tracks how far that has actually got.
