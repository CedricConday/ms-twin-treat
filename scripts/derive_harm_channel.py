"""Re-derive the harm channel's expression table from the Human Protein Atlas.

`bricks/harm_channel.py` carries transcribed nTPM values. This repo's rule is
that a number in a constant ships with a way to reproduce it (see
`derive_suppressive_strength.py`), so:

    python scripts/derive_harm_channel.py            # print the table
    python scripts/derive_harm_channel.py --check    # compare against the module

Downloads ~2.8 MB from proteinatlas.org. Nothing is cached in the repo — the
data policy in `.gitignore` keeps fetched datasets out, and the eight targets'
values are small enough to live in the module as literals with provenance.

Source: https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip
        Human Protein Atlas, consensus immune-cell RNA, CC BY-SA 3.0.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import urllib.request
import zipfile

from bricks.harm_channel import EFFECTOR_SUBSETS, HPA_NTPM, MIN_EFFECTOR_NTPM

URL = "https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip"
KEEP_CELLS = tuple(EFFECTOR_SUBSETS) + (
    "T-reg", "naive B-cell", "memory B-cell", "NK-cell")


def fetch() -> dict[str, dict[str, float]]:
    with urllib.request.urlopen(URL, timeout=180) as resp:
        blob = resp.read()
    out: dict[str, dict[str, float]] = {g: {} for g in HPA_NTPM}
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        name = next(n for n in z.namelist() if n.endswith(".tsv"))
        with z.open(name) as fh:
            text = io.TextIOWrapper(fh, encoding="utf-8")
            for row in csv.DictReader(text, delimiter="\t"):
                gene, cell = row["Gene name"], row["Immune cell"]
                if gene in out and cell in KEEP_CELLS:
                    out[gene][cell] = float(row["nTPM"])
    return out


def ratio(v: dict[str, float]) -> float | None:
    eff = [v[c] for c in EFFECTOR_SUBSETS if c in v]
    m = sum(eff) / len(eff) if eff else 0.0
    return (v.get("T-reg", 0.0) / m) if m >= MIN_EFFECTOR_NTPM else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="compare the live download against the module's literals")
    args = ap.parse_args()

    live = fetch()

    if not args.check:
        print(f"{'target':<11} {'T-reg':>10} {'effT mean':>10} {'ratio':>8}")
        print("-" * 43)
        for gene, v in sorted(live.items(), key=lambda kv: -(ratio(kv[1]) or -1)):
            r = ratio(v)
            shown = f"{r:>8.2f}" if r is not None else "   below floor"
            print(f"{gene:<11} {v.get('T-reg', 0.0):>10.1f} "
                  f"{sum(v[c] for c in EFFECTOR_SUBSETS if c in v) / 4:>10.1f} {shown}")
        return 0

    drift = []
    for gene, stored in HPA_NTPM.items():
        for cell, value in stored.items():
            got = live.get(gene, {}).get(cell)
            if got is None:
                drift.append(f"{gene}/{cell}: missing from the live download")
            elif abs(got - value) > 0.05:
                drift.append(f"{gene}/{cell}: module {value}, live {got}")
    if drift:
        print("HPA values have drifted from the transcribed table:")
        for d in drift:
            print("  -", d)
        print("\nHPA versions its releases; update bricks/harm_channel.py and say "
              "which release, rather than silently accepting the new numbers.")
        return 1
    print(f"all {sum(len(v) for v in HPA_NTPM.values())} transcribed values match "
          "the live Human Protein Atlas download.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
