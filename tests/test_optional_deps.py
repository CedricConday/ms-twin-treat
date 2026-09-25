"""The core install must import without the real-data extras.

scanpy lives in requirements-data.txt, not requirements.txt, and the modules
that use it import it on demand. These tests block scanpy in a subprocess, so
they check that contract whether or not scanpy happens to be installed.
"""

import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# sys.modules[name] = None makes any later `import name` raise ImportError.
_BLOCK_SCANPY = "import sys; sys.modules['scanpy'] = None; "


def _run(code: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-c", _BLOCK_SCANPY + code],
                          cwd=REPO, capture_output=True, text=True, timeout=120)


@pytest.mark.parametrize("module", ["data.kang", "bricks.grn"])
def test_real_data_modules_import_without_scanpy(module):
    r = _run(f"import {module}")
    assert r.returncode == 0, r.stderr


def test_missing_scanpy_fails_with_the_install_instruction():
    r = _run("from data.kang import _scanpy; _scanpy()")
    assert r.returncode != 0
    assert "requirements-data.txt" in r.stderr
