# HANDOFF — cbk → vps

**Written 2026-08-19. `cbk` is off the bricks. You own everything that is left.**

Reason: `cbk` works over SSH on a slow link, so every file crosses the wire
twice. You are on the box. You are roughly 60x faster at this and the numbers
say so — you finished four bricks while `cbk` finished two. `cbk` moves to
review: reading your commits, checking claims against what the code does, and
guarding the honesty rules.

Take everything below. Claim it the normal way (`mkdir claims/<B>`), then build.

---

## 1. WHAT IS DONE (do not rebuild)

| Brick | File | Owner | State |
|---|---|---|---|
| SPINE | `spine/pipeline.py` | cbk | done |
| WIRE | `spine/run_demo.py` | cbk | done — **now yours to keep current** |
| B2 cell | `bricks/cell_transfer.py` | vps | done |
| B3 GRN | `bricks/grn.py` | vps | done |
| B4 QSP | `bricks/qsp.py` | vps | done |
| B5 ABM | `bricks/abm.py` | vps | done |
| B6 barrier | `bricks/barrier.py` | cbk | done |
| B7 intervention | `bricks/intervention.py` | cbk | done |
| **B8 readout** | `bricks/readout.py` | **YOURS** | **not built** |
| **B9 VPop** | `bricks/vpop.py` | **YOURS** | **not built** |

`python spine/run_demo.py` runs end to end right now, 7 stages, no dataset
needed. `--with-data` adds B2/B3 against the Kang matrix. `--arm` switches
treatment arm. Keep that script running at every commit — it is the demo.

---

## 2. THREE INTEGRATION BUGS ALREADY FOUND AND FIXED

These are written down because each one was silent — nothing crashed, the
pipeline just quietly lied. Same class of bug will happen again in B8/B9.

**(a) Stage interface mismatch.** `cbk` first wrote the spine calling
`stage(state)`; your `BUILD_PLAN` §3(b) specifies `.run(state)`. Your bricks
would not have plugged in. Spine now calls `.run()` first and falls back to
`__call__`. **Write `.run(state)`.**

**(b) Intervention was invisible to your models.** `run_demo` originally wrote
`state["intervention"]` through `standin()`, which wraps as
`{"value": {...}, "STANDIN": True}`. Your QSP and ABM read
`state["intervention"]["treat"]` **directly** — the lookup missed, the default
applied, and every run was silently untreated with no error anywhere. B7 now
writes a **plain dict** with `validated: False` inside it. **Any brick whose
output is read by key must do the same — never wrap it.**

**(c) Toy output reported as validated.** `abm_damage` is a bare ndarray; its
`validated: False` lives in the sibling `abm_meta`. The report counted it as a
validated-path key, i.e. the demo overclaiming — exactly what §5 forbids.
`Pipeline._unvalidated` now checks the value, `<key>_meta`, **and**
`<brick-prefix>_meta`. **Keep writing `<prefix>_meta` with `validated: False`.**

---

## 3. B8 — CLINICAL READOUT (what "good" looks like)

Reads `abm_damage` and `cns_exposure`, writes `readout`.

- Map the damage trajectory to a lesion-count proxy and a relapse-rate proxy.
- Gate the treatment effect by `cns_exposure["effective"]` — that is the whole
  reason B6 exists. A peripheral agent gets `effective = 1.0`; a CNS-required
  agent gets the barrier fraction (~0.08). If B8 ignores it, B6 is decoration.
- Write `readout_meta = {"validated": False, ...}`.
- **The micro→clinical map is open research.** Label it hard. This is the single
  easiest place in the whole pipeline to accidentally claim a clinical result.

## 4. B9 — VPOP (the wedge, build this even if B8 slips)

Currently faked by `make_cohort()` in `run_demo.py`, which only varies seed and
`bbb_disruption`. That is not a virtual population and the docstring says so.

- Latin-hypercube or rejection sampling of parameter sets against plausibility
  bounds (Allen–Rieger flavour). MAPEL prevalence-weighting as a documented TODO.
- Return `list[dict]` that `Pipeline.run_cohort()` can consume directly.
- Replace `make_cohort()` in `run_demo.py` when it lands.
- **This is the differentiator** — per `ms-twin/docs/RESEARCH_FINDINGS.md` no
  open Python implementation of the plausible-patient method exists on GitHub at
  all. If time runs short, B9 beats B8.

---

## 5. UNCHANGED GUARDRAILS

- Author **always** `cedric@condaydigital.com`. Never the gmail.
- **Push held until Cedric says go.** Local commits only.
- Every unvalidated value flagged `validated: False`, in code and in output.
- Selftest + Kang green before each commit: `python -m backtest.selftest`,
  `python -m backtest.run_kang` (module form — NOT `python backtest/selftest.py`).
- **Built ≠ validated.** Nothing in this repo is evidence about MS.

---

## 6. NUMBERS — INDEPENDENTLY REPRODUCED

`cbk` re-ran the cold-boot checks from scratch. All three claims hold:

```
identity      delta_pearson = 0.0000    (claimed 0.00)  OK
mean-shift    delta_pearson = 0.8498    (claimed 0.85)  OK
Megakaryocytes              = 0.476     (claimed 0.48)  OK
HARNESS VALIDATED
```

These are safe to say out loud, including to an interviewer. Nothing else in
this repo is.
