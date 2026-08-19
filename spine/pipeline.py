"""The spine: a multi-scale pipeline that runs stages in order over a shared state.

Design notes, because the shape of this file constrains every brick that plugs
into it:

`MultiScaleState` is a plain dict, not a class. Each stage reads the keys it
needs and writes the keys it produces. That is deliberately loose — the science
is unsettled, and a rigid schema would need renegotiating every time a brick
learns something. The cost of looseness is that a typo becomes a missing key
rather than a type error, so `Pipeline.run` reports exactly which keys each
stage added and every stage declares what it requires up front.

The STANDIN discipline is enforced here rather than left to convention. Most of
this pipeline is scaffolding around components nobody has validated, and a
pipeline that runs is not a model that works. `standin()` marks a value as
unvalidated, and the run report counts them, so "it ran" can never be quietly
read as "it works".
"""

from __future__ import annotations

from typing import Any, Iterable, Protocol, runtime_checkable

MultiScaleState = dict[str, Any]

STANDIN_FLAG = "STANDIN"


# --------------------------------------------------------------------------- #
# stand-in marking
# --------------------------------------------------------------------------- #

def standin(value: Any, reason: str = "") -> dict:
    """Wrap a value that is NOT validated science.

    Every brick that cannot yet compute its key for real writes its output
    through this. The wrapper travels with the data, so a downstream consumer
    can see what it is standing on, and the run report can count what is real.
    """
    out = {"value": value, STANDIN_FLAG: True}
    if reason:
        out["reason"] = reason
    return out


def is_standin(value: Any) -> bool:
    """True if a value is not validated science.

    Two conventions are in use and both count, because counting only one of
    them would let the run report call a toy model a result:

      1. `standin(x)`            -> {"value": x, "STANDIN": True}
         A placeholder for a brick that is not built.

      2. `{..., "validated": False}`
         The bricks' own convention (see bricks/qsp.py, abm.py, grn.py). The
         computation is real; the parameterisation is illustrative. Output
         travels as a normal dict so downstream code can read its keys
         directly.

    Form 2 is preferred for a brick that actually computes something, because
    a consumer can do `state["intervention"]["treat"]` without unwrapping.
    """
    if not isinstance(value, dict):
        return False
    return value.get(STANDIN_FLAG) is True or value.get("validated") is False


def unwrap(value: Any) -> Any:
    """Read a value whether or not it is wrapped by `standin()`.

    Only form 1 is a wrapper. Form 2 IS the value, so it passes through.
    """
    if isinstance(value, dict) and value.get(STANDIN_FLAG) is True and "value" in value:
        return value["value"]
    return value


# --------------------------------------------------------------------------- #
# stages
# --------------------------------------------------------------------------- #

@runtime_checkable
class Stage(Protocol):
    """One scale of the model. Canonical interface is BUILD_PLAN.md §3(b):

        class Stage(Protocol):
            name: str
            def run(self, state: dict) -> dict: ...

    `requires` is optional and additive — it lets the pipeline fail at the
    boundary with a useful message rather than deep inside a brick's numerics,
    where a missing upstream key surfaces as an unrelated exception. A brick
    that does not declare it simply gets no pre-check.
    """

    name: str

    def run(self, state: MultiScaleState) -> MultiScaleState: ...


def invoke(stage: Any, state: MultiScaleState) -> MultiScaleState:
    """Call a stage by its declared interface.

    BUILD_PLAN §3(b) specifies `.run(state)`. Plain callables are also accepted
    so a lambda or a partial can stand in during development without needing a
    class. `.run()` wins when both exist.
    """
    runner = getattr(stage, "run", None)
    if callable(runner):
        return runner(state)
    if callable(stage):
        return stage(state)
    raise StageError(
        f"stage {getattr(stage, 'name', stage)!r} has neither .run(state) nor __call__"
    )


class StageError(RuntimeError):
    pass


class MissingRequirement(StageError):
    pass


# --------------------------------------------------------------------------- #
# pipeline
# --------------------------------------------------------------------------- #

class Pipeline:
    """Runs stages in order, threading one state through all of them.

    Stages mutate-and-return the state. The runner records what each stage
    added so a run is legible without stepping through a debugger.
    """

    def __init__(self, stages: Iterable[Stage], name: str = "ms-twin") -> None:
        self.stages = list(stages)
        self.name = name
        self.trace: list[dict] = []

    def _check(self, stage: Stage, state: MultiScaleState) -> None:
        missing = [k for k in getattr(stage, "requires", ()) if k not in state]
        if missing:
            raise MissingRequirement(
                f"stage {stage.name!r} requires {missing} which no earlier stage wrote. "
                f"state has: {sorted(state)}"
            )

    def run(self, state: MultiScaleState | None = None, verbose: bool = True) -> MultiScaleState:
        state = dict(state or {})
        self.trace = []
        for i, stage in enumerate(self.stages, 1):
            before = set(state)
            self._check(stage, state)
            result = invoke(stage, state)
            if not isinstance(result, dict):
                raise StageError(
                    f"stage {stage.name!r} returned {type(result).__name__}, expected the state dict"
                )
            state = result
            added = sorted(set(state) - before)
            flagged = [k for k in added if is_standin(state[k])]
            self.trace.append({"stage": stage.name, "added": added, "standin": flagged})
            if verbose:
                marks = "".join(" *" if k in flagged else "  " for k in added)
                keys = ", ".join(f"{k}{'*' if k in flagged else ''}" for k in added) or "(nothing)"
                print(f"  [{i}/{len(self.stages)}] {stage.name:<28} -> {keys}")
        return state

    def run_cohort(self, states: Iterable[MultiScaleState], verbose: bool = False) -> list[MultiScaleState]:
        """Run one state per virtual patient. The VPop brick produces these."""
        out = []
        for i, s in enumerate(states, 1):
            if verbose:
                print(f"  -- cohort member {i}")
            out.append(self.run(s, verbose=verbose))
        return out

    # ----------------------------------------------------------------- report
    @staticmethod
    def _unvalidated(state: MultiScaleState, key: str) -> bool:
        """Is this key unvalidated, by any of the conventions in use?

        A brick may return a bare array and put its provenance in a sibling
        `<key>_meta` — bricks/abm.py writes `abm_damage` (ndarray) alongside
        `abm_meta = {"validated": False}`. Checking only the value itself would
        report a toy model's output as a validated-path key, which is precisely
        the overclaim BUILD_PLAN §5 forbids.
        """
        if is_standin(state[key]):
            return True
        # Bricks name their meta by BRICK prefix, not by key:
        #   abm_damage -> abm_meta,  cell_delta -> cell_meta,  grn_edges -> grn_meta
        # so check both forms. Checking only `<key>_meta` silently reports a toy
        # model's raw ndarray output as validated.
        for meta_key in (f"{key}_meta", f"{key.split('_')[0]}_meta"):
            meta = state.get(meta_key)
            if isinstance(meta, dict) and meta.get("validated") is False:
                return True
        return False

    def report(self, state: MultiScaleState) -> str:
        """Honest summary. Counts stand-ins so 'it ran' cannot be read as 'it works'."""
        keys = [k for k in sorted(state) if not k.endswith("_meta")]
        stand = [k for k in keys if self._unvalidated(state, k)]
        real = [k for k in keys if k not in stand]
        lines = [
            f"pipeline: {self.name}",
            f"stages run: {len(self.trace)}",
            f"state keys: {len(keys)}",
            "",
            f"  validated-path keys ({len(real)}): {', '.join(real) or '(none)'}",
            f"  STANDIN keys      ({len(stand)}): {', '.join(stand) or '(none)'}",
        ]
        if stand:
            lines += [
                "",
                "  NOTE: STANDIN keys are placeholders, not results. This pipeline",
                "  running end to end demonstrates that the scales compose. It does",
                "  NOT demonstrate that any of them are correct.",
            ]
        return "\n".join(lines)


# --------------------------------------------------------------------------- #
# helper for building trivial stages (used by run_demo before bricks land)
# --------------------------------------------------------------------------- #

class PassThroughStage:
    """A labelled placeholder that writes one STANDIN key and nothing else.

    Exists so the whole spine can be wired and run before any brick is written.
    Every one of these is a slot waiting to be replaced by a real brick.
    """

    def __init__(self, name: str, produces: str, value: Any = None,
                 requires: tuple[str, ...] = (), reason: str = "brick not built yet") -> None:
        self.name = name
        self.produces = produces
        self.requires = requires
        self._value = value
        self._reason = reason

    def run(self, state: MultiScaleState) -> MultiScaleState:
        state[self.produces] = standin(self._value, self._reason)
        return state

    # convenience so a stage is usable either way
    __call__ = run
