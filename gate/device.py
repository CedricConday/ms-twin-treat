"""decide(candidate) -> KILL | ABSTAIN | PASS.

The device. It composes two things that already exist and adds one that did
not: the authority to refuse.

    screen.kill_filter      can say "doomed", and needs no predicted effect
    gate.evidence           says whether the predictor may be listened to
    -> this module          turns those into one verdict a person can act on

ORDER OF OPERATIONS, AND WHY
-----------------------------
1. KILL first. A doomed candidate is doomed whatever the predictor's
   certificate says, and killing is the one judgement this model can make on
   its own evidence. Running it first also means the expensive certificate is
   never computed for a candidate that cannot survive.

2. Then ask the certificate. If it does not unlock PASS, every survivor gets
   ABSTAIN with the blocking reason attached -- the reason, not a shrug, so the
   reader can see exactly which measurement would change the answer.

3. Only with an unlocking certificate does a survivor become a PASS candidate.

WHAT `decide` WILL RETURN TODAY, EVERY TIME, FOR EVERY INPUT
--------------------------------------------------------------
KILL or ABSTAIN. Never PASS. Both scorers lose to predict-the-mean, so
`EvidenceCertificate.unlocks_pass` is False and clause 2 catches everything.

This is not a stub and the PASS branch is not dead code: hand `decide` a
certificate whose scorers clear the criterion and it returns PASS. That path is
exercised by `tests/test_gate.py` with a synthetic certificate, precisely so
that the day the real measurement inverts, the device changes its answer with
no edit -- nobody has to remember to turn PASS on, and equally nobody can turn
it on without a certificate.

THE REFUSAL THAT MATTERS MOST
------------------------------
`decide` never returns a score, a rank, or an effect size, and there is no
parameter that makes it. `screen.rank_candidates()` already raises for this
reason. A pass/fail device that also whispered "and this one looks strongest"
would hand back exactly the artifact the LOMO number says is not information.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from backtest.lomo import load as load_table
from bricks.qsp_velez import MechanismProfile
from gate.criterion import CRITERION, AcceptanceCriterion
from gate.evidence import EvidenceCertificate, certify
from gate.provenance import ModelProvenance, combined, from_simulator, from_table
from screen.kill_filter import ScreenResult, screen


class VerdictKind(Enum):
    KILL = "killed: this candidate is doomed for a reason the model can see"
    ABSTAIN = "no verdict: the model is not entitled to judge this candidate"
    PASS = "passed: survived the filters, scored by a predictor that beat its null"


@dataclass(frozen=True)
class Verdict:
    """One candidate's answer, with everything needed to disagree with it.

    `model` is REQUIRED and has no default. Two models this repo cannot tell
    apart -- the transcription and the K = 2000 extension, 45.4pp and 48.0pp,
    a gap inside the noise floor -- disagree about 11 of 40 survivors, so a
    verdict that does not name its model is not reproducible. It is a field rather than a docstring warning because
    this repo has spent a night discovering that docstring warnings do not hold.
    """

    kind: VerdictKind
    candidate: str
    model: ModelProvenance
    reasons: list[str] = field(default_factory=list)
    screen_result: ScreenResult | None = None
    certificate: EvidenceCertificate | None = None
    criterion: AcceptanceCriterion = CRITERION

    @property
    def validated(self) -> bool:
        """Always False. A verdict is about a simulation, not about a patient."""
        return False

    @property
    def is_pass(self) -> bool:
        return self.kind is VerdictKind.PASS

    def line(self) -> str:
        head = f"{self.candidate:<28} {self.kind.name:<8}"
        if self.reasons:
            return f"{head} {self.reasons[0]}"
        return head

    def explain(self) -> str:
        lines = [f"{self.candidate}: {self.kind.name} under {self.model.label()}",
                 f"  {self.kind.value}",
                 f"  model: {self.model.model}"]
        for r in self.reasons:
            lines.append(f"  - {r}")
        lines.append(f"  validated={self.validated} (a verdict about a simulation)")
        return "\n".join(lines)


def _provenance() -> ModelProvenance:
    """The one model this verdict is about, or a refusal if the halves disagree.

    Read from the artifacts that did the work: the simulator's own signature for
    the kill filters, and the cached table's recorded parameters for the
    certificate. Raises `ModelMismatch` rather than picking one.
    """
    return combined(from_simulator(), from_table(load_table()))


def decide(profile: MechanismProfile,
           certificate: EvidenceCertificate | None = None,
           criterion: AcceptanceCriterion = CRITERION,
           target_gene: dict[str, str] | None = None) -> Verdict:
    """The device. Pass a certificate to avoid re-measuring it per candidate."""
    result = screen([profile], target_gene=target_gene)[0]
    return _verdict(result, certificate, criterion, _provenance())


def _verdict(result: ScreenResult,
             certificate: EvidenceCertificate | None,
             criterion: AcceptanceCriterion,
             model: ModelProvenance) -> Verdict:
    """Turn one screen result into a verdict. Shared by `decide` and `decide_all`."""
    profile = result.profile

    if not result.survived:
        reasons = [f"{result.killed_by.name}: {result.killed_by.value}"]
        if result.detail:
            reasons.append(result.detail)
        if result.like_existing:
            reasons.append("indistinguishable from: " + ", ".join(result.like_existing))
        return Verdict(kind=VerdictKind.KILL, candidate=profile.label, model=model,
                       reasons=reasons, screen_result=result, criterion=criterion)

    cert = certificate if certificate is not None else certify(criterion)

    if not cert.unlocks_pass:
        reasons = ["survived every kill filter, which is not the same as working"]
        reasons += [f"PASS unavailable — {r}" for r in cert.blocking_reasons()]
        return Verdict(kind=VerdictKind.ABSTAIN, candidate=profile.label, model=model,
                       reasons=reasons, screen_result=result, certificate=cert,
                       criterion=criterion)

    reasons = ["survived every kill filter"]
    for sc in cert.scorers.values():
        reasons.append(f"predictor evidence — {sc.line().strip()}")
    if result.regulatory_liability is not None:
        reasons.append(f"regulatory liability flag (soft, not a veto): "
                       f"Treg:effector = {result.regulatory_liability:.2f}")
    return Verdict(kind=VerdictKind.PASS, candidate=profile.label, model=model,
                   reasons=reasons, screen_result=result, certificate=cert,
                   criterion=criterion)


def decide_all(profiles: list[MechanismProfile],
               criterion: AcceptanceCriterion = CRITERION,
               target_gene: dict[str, str] | None = None) -> list[Verdict]:
    """Decide a batch against ONE certificate and ONE screen pass.

    Cheaper than looping `decide`, which re-derives the untreated baseline per
    candidate. Deliberately returns the verdicts in the order given -- sorting
    them would be ranking by the back door.
    """
    model = _provenance()
    cert = certify(criterion)
    results = screen(profiles, target_gene=target_gene)
    return [_verdict(r, cert, criterion, model) for r in results]


def main() -> int:
    from screen.kill_filter import enumerate_candidates

    candidates = enumerate_candidates(max_points=1)
    verdicts = decide_all(candidates)

    print("THE ACCEPT/REJECT DEVICE — one verdict per candidate, in input order\n")
    print(f"  model: {verdicts[0].model.line()}\n")
    for v in verdicts:
        print("  " + v.line())

    counts = {k.name: sum(1 for v in verdicts if v.kind is k) for k in VerdictKind}
    print(f"\n  {len(verdicts)} candidates under {verdicts[0].model.label()}: "
          + ", ".join(f"{n} {k}" for k, n in counts.items() if n))
    print("  Survival is model-relative: the transcription and the K=2000 extension")
    print("  score 45.4pp and 48.0pp — a gap inside the noise floor — and disagree on")
    print("  11 of 40 survivors. Quote the label above with any verdict from this run.")

    cert = verdicts[0].certificate if verdicts and verdicts[0].certificate else certify()
    print()
    print(cert.report())
    if not cert.unlocks_pass:
        print("\n  No input can return PASS while that is False. The device is working;")
        print("  the predictor underneath it has not earned the verdict yet.")
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Simulation")
    print("  numbers are proxies from toy models. Nothing here is evidence about")
    print("  multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
