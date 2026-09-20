"""Run the leave-one-mechanism-out scorer and cache it for `gate.evidence`.

LOMO simulates every arm on every paired seed for each fold, which costs minutes.
The device must not pay that per candidate, and it must not silently skip it
either, so the measurement is cached to a file with the date and commit it was
taken at, and `gate/evidence.py` reads that file. No cache means UNMEASURED,
which can never unlock PASS.

Re-run this whenever anything under `bricks/` that LOMO touches changes. A stale
cache is the failure mode here, which is why the commit SHA is recorded: compare
it against HEAD before trusting a certificate that cites it.

Run:  PYTHONPATH=. python3 scripts/cache_lomo_certificate.py
"""

from __future__ import annotations

import json
import subprocess
from datetime import date
from pathlib import Path

from backtest.lomo import run_lomo

OUT = Path(__file__).resolve().parent.parent / "results" / "lomo_certificate.json"


def _commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True,
                              cwd=OUT.parent.parent).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def main() -> int:
    r = run_lomo()
    payload = {
        "measured_on": date.today().isoformat(),
        "commit": _commit(),
        "mae": r["mae"],
        "null_mae": r["null_mae"],
        "n_groups": r["n_groups"],
        "folds": r["folds"],
        "note": "written by scripts/cache_lomo_certificate.py; read by gate/evidence.py. "
                "Simulation numbers are proxies from a ported toy model.",
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(OUT.parent.parent)}: "
          f"MAE {r['mae']:.1f}pp vs null {r['null_mae']:.1f}pp over {len(r['folds'])} folds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
