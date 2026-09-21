"""gate/ — the accept/reject device: does a simulated therapy pass, or not?

`screen/` kills doomed candidates. `backtest/` measures whether the model can
predict an arm it was not fitted on. Neither one answers the question a person
outside this repo actually asks:

    "You simulated a therapy. Is it a pass?"

This package answers it, and the answer today is **never PASS** — not because
the verdict is hard-coded, but because PASS is gated on a measurement the repo
does not currently make. That distinction is the whole design.

THE DEVICE
----------
    decide(profile) -> Verdict(KILL | ABSTAIN | PASS, reasons, evidence)

KILL     the candidate is doomed for a reason that needs no predicted effect
         size: out of regime, unreachable at any potency, or degenerate with a
         drug that already exists. Delegated to `screen.kill_filter`, which
         could already say this.

ABSTAIN  the candidate survived the kill filters, and the model is not entitled
         to say anything further about it. This is the honest verdict for
         everything that is not killed, and it is what the device returns today
         for every survivor.

PASS     the candidate survived AND the predictor that would score it has been
         shown, out of sample, to carry drug-specific information. Reachable
         only when an `EvidenceCertificate` says so.

WHY PASS IS UNREACHABLE TODAY, IN ONE NUMBER
---------------------------------------------
Both of this repo's out-of-sample scorers lose to predict-the-mean, and not
narrowly -- the paired bootstrap CI on each lies entirely ABOVE zero, meaning
the predictor is reliably worse than the null rather than close to it:

    backtest/loo.py    hold out an ARM        30.0pp MAE vs 11.3pp null
                                              95% CI on the paired gap [+8.9, +29.3]
    backtest/lomo.py   hold out a MECHANISM   45.4pp MAE vs 11.8pp null
                                              95% CI on the paired gap [+20.3, +42.8]

Both re-measured 2026-09-21 on the 15-quantified-arm exam, LOMO via
`results/lomo_certificate.json`, which now records the arm set it was taken on
and is refused if that set has changed.
The carrying-capacity extension to the QSP model (ec4ee5b) was the live candidate
for lifting this and does not: 45.6pp against the same null, measured by the
session that owns the master checkout. It fixes the depletion SIGN and not the
magnitude.

Do not read those numbers from here. `gate.evidence.certify()` recomputes them,
and every figure in `BUILD_PLAN.md` older than 2026-09-20 predates the arm set
growing to 17 arms / 12 quantified -- stale by one exam.

A scoring rule built on a predictor that loses to the mean of its training arms
is a rule about the training arms' average, not about the therapy in front of
it. So a PASS emitted today would be a statement about the null. The device
refuses to make it.

WHAT THIS IS NOT
----------------
It is not a clinical decision, not a regulatory artefact, and not evidence about
multiple sclerosis. A PASS from this device -- if the certificate ever unlocks
one -- would mean "this candidate survived the filters this model can see, and
the model that scored it beat its null on held-out data". That is a statement
about a simulation. Every verdict carries `validated = False`.

THE POINT OF BUILDING IT NOW
-----------------------------
The criterion is written down BEFORE the measurement clears it. Writing an
acceptance rule after seeing which candidates you like is how a screen launders
a preference into a result. `gate/criterion.py` is dated, frozen, and cites the
standard it comes from -- which is this repo's own, already stated in
`screen/__init__.py` ("Until that inverts, a predicted effect size is not
information") and `docs/QUALITY.md` item 5 (a null beside every benchmark).
Nothing here is a new threshold invented tonight.
"""

# Deliberately no eager re-export, for the reason `screen/__init__.py` gives:
# importing the submodules here makes `python -m gate.evidence` emit a
# RuntimeWarning about the module already being in sys.modules, and a warning
# nobody can act on trains people to ignore warnings. Import from the submodule:
#
#     from gate.device import decide, Verdict, VerdictKind
#     from gate.evidence import certify
#     from gate.criterion import CRITERION
