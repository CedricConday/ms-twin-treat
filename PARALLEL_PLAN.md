# PARALLEL_PLAN — two instances, one repo, no collisions

**Read this before writing any code. Both instances work in `~/repos/ms-twin-treat/`.**

Instances:
- **`vps`** — Claude Code on vps1 (has the box, does the heavy installs)
- **`cbk`** — Claude Code on the Chromebook, working over SSH into the same repo

Supersedes nothing in `BUILD_PLAN.md` except the build order in §4 — see §3 below.
All guardrails in `BUILD_PLAN.md` §5 still apply, unchanged.

---

## 1. THE CLAIM PROTOCOL (non-negotiable, do this first every time)

Claiming is done with `mkdir`, because `mkdir` is atomic on POSIX — exactly one
caller can create a given directory. Two instances racing for the same brick,
one wins, the other gets an error and moves on. No coordination service, no
polling, no lost work.

**Before you write a single line for a brick:**

```bash
cd ~/repos/ms-twin-treat
mkdir claims/B4 2>/dev/null && {
    echo "owner: vps"            >  claims/B4/OWNER   # or cbk
    echo "started: $(date -u +%FT%TZ)" >> claims/B4/OWNER
    echo "claimed B4"
} || {
    echo "B4 ALREADY CLAIMED BY: $(cat claims/B4/OWNER 2>/dev/null)"
    echo "pick another brick"
}
```

If the claim fails, **you do not touch that brick.** Pick the next unclaimed one
from your column in §2.

**When a brick is finished:**

```bash
echo "done: $(date -u +%FT%TZ)" >> claims/B4/OWNER
git add -A && git commit -m "feat(b4): <what it does>"
```

**To see the board at any time:**

```bash
for d in claims/*/; do printf "%-6s %s\n" "$(basename $d)" "$(head -1 $d/OWNER)"; done
```

Claims are committed to git, so `git log claims/` is the audit trail of who
built what and when.

---

## 2. OWNERSHIP — default split

Split by dependency on heavy installs, not by count. `vps` owns everything that
needs a big pip install; `cbk` owns everything that is pure Python + scipy, so
neither instance ever waits on the other's toolchain.

| Brick | File (owner writes ONLY this) | Owner | Notes |
|---|---|---|---|
| **SPINE** | `spine/pipeline.py` | **cbk** | Stage protocol + MultiScaleState + runner |
| **WIRE** | `spine/run_demo.py` | **cbk** | all 8 stages as STANDIN pass-throughs — **do this 2nd** |
| B2 Cell | `bricks/cell_transfer.py` | **vps** | ridge cross-cell-type transfer; scGPT later |
| B3 GRN | `bricks/grn.py` | **vps** | pyscenic/arboreto; MI-graph fallback |
| B4 QSP | `bricks/qsp.py` | **vps** | tellurium toy ODE |
| B5 ABM | `bricks/abm.py` | **vps** | mesa immune-vs-myelin grid |
| B6 Barrier | `bricks/barrier.py` | **cbk** | scipy 3-compartment ODE, no R |
| B7 Intervention | `bricks/intervention.py` | **cbk** | param object modifying QSP/ABM inputs |
| B8 Readout | `bricks/readout.py` | **cbk** | sim trajectory → clinical proxy |
| B9 VPop | `bricks/vpop.py` | **cbk** | LHS/rejection sampling — **the wedge** |

**`vps` starts its installs at t=0 in the background** (`torch`, `pyscenic`,
`tellurium`, `mesa`) and builds while they run. `cbk` needs no installs beyond
scipy/numpy, so it produces the skeleton immediately.

**File ownership is absolute.** Nobody edits a file in another instance's row.
If you need something changed there, note it in `claims/<brick>/NOTES` and let
the owner make the change.

---

## 3. BUILD ORDER — the one change to BUILD_PLAN §4

`BUILD_PLAN.md` §4 puts `run_demo.py` last. That means nothing runs end to end
until the very end, so if the clock stops early there are good bricks and no
demo. A prefix of that plan is not a demonstration.

**Reordered:**

1. `cbk` — **SPINE** (`spine/pipeline.py`) — ~20 min
2. `cbk` — **WIRE** (`spine/run_demo.py`) with all 8 stages as trivial
   pass-throughs that write a labelled `STANDIN` key — ~20 min
3. Everything else, in any order, by both instances in parallel

After step 2, `python spine/run_demo.py` runs end to end **and keeps running
after every subsequent commit.** Each brick *replaces* a stand-in rather than
*adding* a missing piece. Whenever you stop, there is a working demo.

`vps` does not wait for steps 1–2 — start installs and build brick internals
against the contract in §4; they slot in when the spine lands.

---

## 4. THE CONTRACT (freeze this; do not renegotiate mid-build)

`MultiScaleState` is a plain `dict[str, Any]` passed stage to stage. Each stage
reads the keys it needs and writes the keys it produces. A stage that cannot
compute a key yet writes a value clearly marked `STANDIN`.

```python
# spine/pipeline.py
from typing import Any, Protocol

MultiScaleState = dict[str, Any]

class Stage(Protocol):
    name: str
    def __call__(self, state: MultiScaleState) -> MultiScaleState: ...
```

Reserved keys — **owner of the brick owns the key**:

| Key | Written by | Type |
|---|---|---|
| `intervention` | B7 | `dict` — dose, target, on/off |
| `cell_delta` | B2 | `dict[cell_type, np.ndarray]` — predicted expression delta |
| `grn_edges` | B3 | `list[tuple[str, str, float]]` — TF, target, weight |
| `qsp_traj` | B4 | `dict` with `t: np.ndarray`, `y: np.ndarray` |
| `abm_damage` | B5 | `np.ndarray` — damage over time |
| `cns_exposure` | B6 | `float` — fraction crossing the barrier, 0–1 |
| `readout` | B8 | `dict` — `lesion_proxy`, `relapse_proxy` |
| `vpop` | B9 | `list[dict]` — sampled parameter sets |

Every STANDIN value must be labelled in **both** the data and the printed
output, e.g. `state["cns_exposure"] = {"value": 0.12, "STANDIN": True}`.

---

## 5. GUARDRAILS (from BUILD_PLAN §5 — unchanged, restated so nobody misses them)

- Author identity **always** `cedric@condaydigital.com`. **NEVER** the gmail.
- **Push stays held until Cedric says go.** Local commits only. No remote.
- No datasets, weights or secrets in git.
- **Built ≠ validated.** Every stand-in labelled `STANDIN` in code and output.
  The demo running is not the model working. Under-promise to the MS community.
- Commit after each brick, so a pause never costs more than one brick.
- Harness selftest + Kang run must stay green. Re-run before you commit.

---

## 6. DEFINITION OF DONE

`python spine/run_demo.py` runs a small virtual cohort through all stages
without error, printing the state accumulating keys at each brick, with every
unvalidated component labelled `STANDIN`. Selftest + Kang still green. Committed
brick by brick. Push still held.

---

## 7. IF YOU ARE READING THIS AFTER A RESET

1. `cd ~/repos/ms-twin-treat && git log --oneline` — see what exists.
2. Run the board command in §1 — see what is claimed and by whom.
3. Run `BUILD_PLAN.md` §1 cold-boot checks — confirm harness + Kang green.
4. Claim the next unclaimed brick **in your own column** and build it.

Do not rebuild something with a `done:` line in its `claims/*/OWNER`.
