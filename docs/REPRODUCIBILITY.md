# Reproducibility

Run everything from the repository root, in a virtual environment. The pins in
`requirements.txt` are a coherent set; a system Python with a newer numpy or
pandas will import and mostly run, which is worse than failing, because the
numbers it produces are not the ones in `results/RESULTS.md`.

## Environment

```bash
python -m venv .venv
. .venv/bin/activate                 # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt
```

That is enough for the test suite, the harness self-test, the toy pipeline and
the lint gate. It is deliberately **not** enough for the real-data path.

## The real-data path

`data/kang.py` needs scanpy, which pulls numba, umap-learn, seaborn and
statsmodels. It is kept out of `requirements.txt` so the common case stays
light, and imported lazily so a missing install fails with an instruction
rather than a traceback:

```bash
python -m pip install -r requirements-data.txt
```

scanpy is pinned to 1.12 — the newest release still compatible with
`pandas==2.2.3`. 1.12.1 and later require `pandas>=2.3`, which breaks the
pinned stack.

## Gates

```bash
python -m ruff check .               # lint
python -m pytest                     # invariant tests + input guards
python -m backtest.selftest          # the ruler validates itself (no data needed)
python spine/run_demo.py --quiet     # toy pipeline, end to end, self-labelling
```

With the real-data deps installed, and a first run that downloads the Kang
matrix into the git-ignored `data/cache/`:

```bash
python -m backtest.run_kang          # real IFN-β backtest, both nulls + ceiling
python -m bricks.cell_scgpt          # the cell model's LOCTO scorecard
python spine/run_demo.py --with-data --quiet
python -m results.experiment         # virtual cohort across real trial arms
python -m backtest.clinical          # the clinical gate
```

## For a comparable result

- use the documented seed arguments;
- record the commit SHA, the Python version, and which requirements file was
  installed;
- do not commit downloaded datasets or generated caches — `data/cache/` is
  git-ignored and stays that way;
- keep the complete command and configuration with any number you quote;
- treat `results/figures/` and cached outputs as derived artefacts, not sources.

A successful run means the software completed. It does not establish biological
or clinical validity; see [`QUALITY.md`](QUALITY.md).
