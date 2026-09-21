"""Which model produced this verdict? Read it off the artifact, never declare it.

A candidate's survival is not a model-independent fact. The published
transcription and the carrying-capacity extension at K = 2000 score 45.4pp and
48.0pp on the current 15-arm exam. That 2.6pp gap sits inside this repo's own
~10% noise floor for these figures (`backtest/lomo.py`), so no measurement here
can prefer one model over the other -- and they disagree about 11 of 40
survivors. Nine patterns survive the transcription and are killed
by the extension, two go the other way.

So "this candidate survived" means nothing without the model that said so, and
the model is not named by its family. K = 50000 and K = 2000 are both "the
capacity extension" and behave oppositely on the same dial: +109% and -3%. A
verdict tagged "extension" without its K is as unreproducible as an untagged one.

TWO PATHS, AND THEY CAN DISAGREE WITHOUT ANYONE NOTICING
---------------------------------------------------------
`gate.decide` reads from two places that are built differently:

    the KILL filters   simulate live through `bricks/qsp_velez.simulate`, at
                       whatever parameters that call passes
    the certificate    reads CACHED response tables built by an earlier run,
                       which record their own parameters in the file

Nothing has ever checked that those two describe the same model. They do today
-- both are the published transcription -- but a verdict that killed a candidate
under one model and withheld PASS under another would be incoherent, and would
look exactly like a working verdict. `combined()` compares them and refuses to
produce a provenance when they differ.

WHY THIS IS READ AND NOT PASSED IN
-----------------------------------
A `model=` argument is a label someone types, and a label someone types is the
thing that goes stale. `from_table` reads `carrying_capacity` out of the table
file itself; its absence means no cap, which is the published model. The
simulator's provenance comes from the default of the keyword actually used. Both
are derived from the artifact that did the work.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass, field

from bricks import qsp_velez

# The published model this repo transcribes. An extension is anything that adds
# a term the paper does not have, which today means the carrying capacity.
PUBLISHED = "Velez de Mendizabal 2011 (transcription)"
EXTENDED = "Velez de Mendizabal 2011 + carrying-capacity extension"


@dataclass(frozen=True)
class ModelProvenance:
    """The model AND its parameters, in a form that can be printed beside a verdict."""

    model: str
    parameters: dict = field(default_factory=dict)
    is_extension: bool = False
    source: str = ""

    #: Short names for parameters that appear in a label. A parameter with no
    #: short name is printed in full rather than dropped -- an unnamed parameter
    #: silently vanishing from a label is the failure this class exists to stop.
    SHORT = {"carrying_capacity": "K"}

    def label(self) -> str:
        """Short, sortable, and unambiguous: 'velez2011' or 'velez2011+K=2000'."""
        base = "velez2011"
        if not self.parameters:
            return base
        extras = ",".join(
            f"{self.SHORT.get(k, k)}={v:g}" if isinstance(v, float)
            else f"{self.SHORT.get(k, k)}={v}"
            for k, v in sorted(self.parameters.items()))
        return f"{base}+{extras}"

    def line(self) -> str:
        return f"{self.label():<20} {self.model}" + (f"   [{self.source}]" if self.source else "")


def from_table(table: dict, source: str = "cached response table") -> ModelProvenance:
    """Provenance of a cached response table, read from the file's own metadata.

    A table with no `carrying_capacity` key was built without a cap, which is the
    published model. That is the convention `results/mechanism_curve.json` (no
    key) and `results/mechanism_curve_K2000.json` (`carrying_capacity: 2000.0`)
    already use.
    """
    k = table.get("carrying_capacity")
    if k is None:
        return ModelProvenance(model=PUBLISHED, parameters={}, is_extension=False,
                               source=source)
    return ModelProvenance(model=EXTENDED, parameters={"carrying_capacity": float(k)},
                           is_extension=True, source=source)


def from_simulator(**call_kwargs) -> ModelProvenance:
    """Provenance of a LIVE simulation, from the keyword actually in force.

    Reads `simulate`'s signature rather than assuming its default, so that
    changing the default in `bricks/qsp_velez.py` cannot silently make every
    verdict misdescribe itself.
    """
    sig = inspect.signature(qsp_velez.simulate)
    default_k = sig.parameters["carrying_capacity"].default
    k = call_kwargs.get("carrying_capacity", default_k)
    if k is None:
        return ModelProvenance(model=PUBLISHED, parameters={}, is_extension=False,
                               source="live simulation")
    return ModelProvenance(model=EXTENDED, parameters={"carrying_capacity": float(k)},
                           is_extension=True, source="live simulation")


class ModelMismatch(RuntimeError):
    """The kill filters and the certificate describe different models."""


def combined(screen_prov: ModelProvenance, cert_prov: ModelProvenance) -> ModelProvenance:
    """One provenance for a verdict, or a refusal if the two halves disagree.

    Raising is the right behaviour rather than picking one or recording both: a
    verdict assembled from two models is not a verdict about either, and the
    repo has spent long enough removing numbers that were quietly about
    something other than what they said.
    """
    if (screen_prov.model, screen_prov.parameters) != (cert_prov.model, cert_prov.parameters):
        raise ModelMismatch(
            "the kill filters and the evidence certificate describe different models: "
            f"{screen_prov.label()} ({screen_prov.source}) vs "
            f"{cert_prov.label()} ({cert_prov.source}). A verdict assembled from both "
            "is about neither. Rebuild the cached table under the model the filters "
            "run, or run the filters under the table's model.")
    return ModelProvenance(model=screen_prov.model, parameters=dict(screen_prov.parameters),
                           is_extension=screen_prov.is_extension,
                           source=f"{screen_prov.source} + {cert_prov.source}")
