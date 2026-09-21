# Recoverability: can an independent source find the potency the model needs?

Regenerate with:

```
PYTHONPATH=. python3 -m gate.ceiling
```

Measured 2026-09-21 on the 15-quantified-arm exam. If a number here disagrees
with what that command prints, the command is right and this page is stale.
`scripts/verify_gate_docs.py` checks that mechanically.

## The question

`bricks/profiles.py` says which parameters each drug moves and in which
direction, from pharmacology. It leaves **how much** as one shared stub, because
a hand-picked number per drug is indistinguishable from fitting, and fitting to
a drug's own relapse number is circular.

So the magnitude has to come from somewhere independent. `backtest/potency.py`
supplies one: each drug's potency fitted to its trial's **MRI lesion ratio**,
never to its relapse number. The question this page answers is whether that
independent estimate lands where the model would need it to.

The reference it is measured against is the **oracle** potency from
`gate/ceiling.py` — the potency fitted in-sample to the arm's own relapse
number. That is cheating on purpose: it is the best any source could possibly
supply, so it is the right yardstick for one that is not cheating.

## Result: it does not

| arm | oracle *s* | MRI-fitted *s* | ratio | predicted at MRI *s* | error | null error |
|---|---|---|---|---|---|---|
| ocrelizumab | 0.30 | 0.60 | 2.00 | −73.8% | 27.8pp | 0.6pp |
| dimethyl fumarate | 0.74 | 0.80 | 1.08 | −62.1% | 9.1pp | 8.1pp |
| glatiramer acetate | 0.04 | 0.05 | 1.25 | −38.4% | 9.4pp | 17.6pp |
| teriflunomide | 0.48 | 0.75 | 1.55 | −54.2% | 22.7pp | 14.9pp |
| IFN-beta | 0.48 | 0.75 | 1.56 | −54.2% | 24.2pp | 16.5pp |
| ofatumumab | 0.33 | 0.80 | 2.46 | −84.9% | 34.9pp | 4.9pp |

**21.4pp against a 10.4pp predict-the-mean null.** The independent channel
loses, on the only six arms where it can be scored at all.

It is also biased, not merely noisy: every ratio is above 1, median 1.55. The
MRI channel systematically reads a drug as more potent than the potency that
reproduces its relapse outcome.

## Two different kinds of "unreachable"

Seven of the fifteen quantified arms have no oracle potency, and they fail for
two distinct reasons that should not be pooled:

**No dial can produce the effect (six arms)** — natalizumab, fingolimod,
ponesimod, daclizumab, cladribine and, since it was wired on 2026-09-21,
ozanimod. `gamma_E` and a lowered `alpha_R` both *raise*
damage in this model, so no potency produces a benefit at any strength. These
are the same five `backtest/potency.py` reports OUT OF RANGE. A direction
failure; no data source of any quality touches it.

**No data exists on this route (alemtuzumab)** — CARE-MS I reports a
lesion-free proportion and a lesion *volume* change, neither convertible to the
count ratio the MRI fitter needs, so alemtuzumab never enters `potency.py` at
all. It is absent from that module's OUT OF RANGE list for that reason, not
because it is reachable. The oracle fits on the arm's own relapse number, which
alemtuzumab does have, and finds the grid edge — so it is **unreachable on both
routes at once**, for two unrelated reasons.

That is why this page counts one more unreachable arm than `backtest/potency.py`
reports OUT OF RANGE. Both are correct about different things.

## The ocrelizumab cross-check, corrected

`backtest/potency.py` calls the ocrelizumab discrepancy "the best estimate this
repo has of how far the potency layer can be trusted" and puts it at 2.1x. That
figure is in **multiplier** space (`ke` = 0.85 from the EAE fit against 0.40 from
OPERA's Gd ratio). The model's fitted quantity is **potency**, `s = 1 − ke`, and
the same two numbers in that space are:

| estimate | `ke` multiplier | potency *s* |
|---|---|---|
| Martinez-Pasamar 2013, EAE flow cytometry | 0.85 | 0.15 |
| MRI channel, OPERA Gd-enhancing ratio | 0.40 | 0.60 |
| oracle — what reproduces OPERA's ARR | 0.70 | 0.30 |

In potency space the two independent estimates are **4.0x apart**, not 2.1x. The
repo has been quoting the smaller of two numbers that describe the same
disagreement, because they were computed in different units.

It also resolves an apparent convergence that is not one. The 2.00 in the table
above (MRI ÷ oracle, potency space) and the 2.1 in `potency.py` (EAE ÷ MRI,
multiplier space) are **not two independent routes landing on the same factor**.
They are different ratios in different parameterisations that happen to be
numerically close, and both share the MRI estimate as one of their two terms.

What the three estimates do show, read consistently, is cleaner than a
convergence: the potency that reproduces the trial sits **almost exactly halfway
between the two independent sources on a log scale** — 2.00x above the EAE fit,
2.00x below the MRI fit. Neither independent source finds it, and they miss in
opposite directions.

## What this decides

`BUILD_PLAN` blockers (4) and (6) — per-drug potency, and the patient-level data
that would improve it — are **closed routes on this model, not deferred ones.**
More or better potency data cannot lift a gate whose independent channel misses
the potency it would need, on the half of the arm set the model can reach at all.

The remaining blocker is the model's form: no CNS compartment and no trafficking,
so depletion, sequestration and transit blockade collapse onto one dial that
cannot turn downward. `bricks/profiles.py` names the only way that lump has ever
broken — someone measured what the drug does to a parameter the model already
has — and notes that reasoning harder about the mechanism has never broken it.

Nothing on this page is evidence about multiple sclerosis. Trial numbers are
real and cited (`docs/TRIAL_ANCHORS.md`); simulation numbers are proxies from a
ported toy model, with the Sormani map applied outside what it was fitted on.
