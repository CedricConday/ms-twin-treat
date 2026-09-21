"""Resolve every PMID in docs/TRIAL_ANCHORS.md and check it is about the subject.

WHY THIS EXISTS
---------------
On 2026-09-21 the parallel session reported that four of five PMIDs it recalled
from memory pointed at unrelated papers, and checked its own by live query
instead. Running the same check over the EXISTING table found one: cladribine's
CLARITY anchor was cited as PMID 20089950, which is

    "Signaling by the high-affinity HDL receptor scavenger receptor B type I."

a single digit away from 20089960, the real CLARITY paper. Cladribine is a
QUANTIFIED arm feeding the clinical gate, so a reader chasing the citation
would have found nothing that supports the number beside it.

WHAT IT CHECKS, AND WHAT IT CANNOT
------------------------------------
It resolves each identifier live and asks whether the title looks like it
belongs in this table at all — MS, EAE, or the myelin/immunology vocabulary the
mechanism citations use. **That is a smoke test, not verification.** It cannot
tell you the ARR in the row matches the paper; only reading the paper does
that. It catches the failure that actually happened: a well-formed identifier
pointing somewhere else entirely.

Network, ~25 requests, a few seconds each. Not part of the test suite: it needs
the internet, and a test that fails when the wifi drops teaches people to skip
tests.

Run:  python scripts/verify_anchors.py
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ANCHORS = Path(__file__).resolve().parent.parent / "docs" / "TRIAL_ANCHORS.md"

# Vocabulary any paper in this table should hit. Deliberately wide: the table
# holds trial reports AND mechanism citations (TNF biology, MBP peptides), and
# a narrow list would flag the latter as errors.
ON_TOPIC = (
    "sclerosis", "relapsing", "myelin", "encephalomyelitis", "demyelinat",
    "oligodendrocyte", "interferon", " ms ", "glatiramer", "natalizumab",
    "fingolimod", "teriflunomide", "fumarate", "ocrelizumab", "ofatumumab",
    "alemtuzumab", "cladribine", "daclizumab", "ponesimod", "ublituximab",
    "ozanimod", "laquinimod", "rituximab", "atacicept", "lymphocyte", "t cell",
    "t-cell", "tnf", "immune", "autoimmun", "brain", "neuro",
)


def resolve(pmid: str) -> tuple[str | None, str]:
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
           + urllib.parse.quote(f"EXT_ID:{pmid} AND SRC:MED")
           + "&format=json&resultType=core")
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            hits = json.load(r).get("resultList", {}).get("result", [])
    except Exception as exc:                      # network, not data
        return None, f"LOOKUP FAILED: {exc}"
    if not hits:
        return None, "NOT RESOLVED — no such record"
    return hits[0].get("pubYear"), hits[0].get("title", "")


def main() -> int:
    text = ANCHORS.read_text()
    pmids = sorted(set(re.findall(r"\b[0-9]{8}\b", text)))
    print(f"{len(pmids)} identifiers in {ANCHORS.name}\n")

    suspect = []
    for pmid in pmids:
        year, title = resolve(pmid)
        low = f" {title.lower()} "
        ok = year is not None and any(w in low for w in ON_TOPIC)
        flag = "    " if ok else " !! "
        print(f"{flag}{pmid}  {year or '----'}  {title[:78]}")
        if not ok:
            suspect.append((pmid, title))
        time.sleep(0.15)

    if suspect:
        print(f"\n{len(suspect)} identifier(s) do not look like they belong here:")
        for pmid, title in suspect:
            print(f"  {pmid}  {title[:70]}")
        print("\nA well-formed PMID pointing at the wrong paper is the failure this")
        print("catches. Check the row against the paper before trusting its number.")
        return 1
    print("\nEvery identifier resolves and every title is on topic.")
    print("This is a smoke test: it does NOT confirm the ARR in a row matches")
    print("the paper. Only reading the paper does that.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
