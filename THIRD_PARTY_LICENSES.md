# Third-party licenses & the model boundary

The **own code** in this repository (`spine/`, `bricks/`, `backtest/`, `data/`,
`results/`, `tests/`) is licensed **Apache-2.0** (see `LICENSE`). This document
records the licenses of everything else, and — the load-bearing part — states in
writing what this repository does and does not distribute.

## The boundary (read this first)
- **This repository bundles no model weights and no datasets.** Loaders fetch
  open data at runtime into a git-ignored cache; nothing is redistributed here.
- **No non-commercial or field-restricted model is a dependency.** The pipeline
  that runs uses only the permissive libraries listed below plus our own code.
- **Restricted models are optional and user-supplied.** Where the roadmap names a
  model under a non-permissive license, this repo ships only an *adapter
  interface* — never the model. If you choose to plug one in, you obtain it
  yourself and comply with *its* license. Our Apache-2.0 grant does not extend to
  it and cannot relicense it.

## Runtime dependencies (all permissive — verified)
| Package | License |
|---|---|
| numpy | BSD-3-Clause |
| scipy | BSD-3-Clause |
| pandas | BSD-3-Clause |
| scikit-learn | BSD-3-Clause |
| anndata | BSD-3-Clause |
| scanpy | BSD-3-Clause |
| mesa | Apache-2.0 |

## Development / optional
| Package | License | Use |
|---|---|---|
| pytest | MIT | test suite |
| matplotlib | Matplotlib/PSF (permissive, BSD-compatible) | `results/figures.py` deck assets |

## Data
| Dataset | Source | Terms | Handling |
|---|---|---|---|
| Kang 2018 (GSE96583, IFN-β PBMCs) | figshare (as used by `pertpy.data.kang_2018`) | per the figshare record (open, CC-BY family) | **downloaded at runtime, git-ignored, not redistributed** |

## Models named in the roadmap — NOT dependencies, NOT bundled
Documented in `../ms-twin/docs/RESEARCH_FINDINGS.md` as options; none is imported
or shipped by this repo. Each carries its own terms, which the *user* must honor
if they choose to integrate it:
| Model | License (as recorded) | Status here |
|---|---|---|
| scGPT (bowang-lab) | MIT (permissive) | candidate cell brick; **not yet integrated** |
| GenBio AIDO / GB.Cell / GB.Tissue … | `license:other` (custom, non-permissive) | **not used** |
| Arc State, CellFM | non-commercial | **not used** |
| PhysiCell / PK-Sim | BSD / GPLv2 respectively | referenced; **not bundled** |

> If any roadmap model's license contains a field-of-use restriction (e.g.
> "research only", "non-commercial", or a medical-use limitation), that
> restriction travels with the model and is **not** affected by this repository's
> Apache-2.0 license. Before integrating such a model into a commercial offering,
> read its exact license text and, where stakes warrant, obtain legal review.

*Last verified: 2026-08-19 against the installed environment.*
