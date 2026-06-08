"""Path + import bootstrap for Energetic Prometheus.

This module makes the two upstream codebases importable from a single place and
isolates their import-time side effects:

* **Prometheus** (``../../Prometheus``) — the radiative-transfer engine.  Its
  scripts use intra-package relative imports (``from . import constants``), so
  the package is imported as ``Prometheus.pythonScripts.gasProperties``.  We add
  the repo root *and* the ``Prometheus`` directory to ``sys.path`` so either the
  fully-qualified or the bare ``pythonScripts.*`` form resolves.

* **dishoom** (``../../dishoom``) — the energetics/escape physics.  Its modules
  print banners (``"dishoom dishoom"``, ASCII art) and pull a tangle of
  ``import *`` modules at import time.  :func:`silence_stdout` wraps those
  imports so nothing leaks onto the caller's stdout.

Nothing here computes physics; it only wires paths.  Import side effects run
exactly once thanks to the module cache.
"""

from __future__ import annotations

import contextlib
import io
import os
import sys

# ── Locate the sibling codebases ────────────────────────────────────────────
# This file lives at  <repo>/Energetic Prometheus/energetic_prometheus/_bootstrap.py
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_PKG_DIR, "..", ".."))
PROMETHEUS_DIR = os.path.join(REPO_ROOT, "Prometheus")
DISHOOM_DIR = os.path.join(REPO_ROOT, "dishoom")

for _p in (REPO_ROOT, PROMETHEUS_DIR, DISHOOM_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)


@contextlib.contextmanager
def silence_stdout():
    """Suppress stdout (dishoom prints banners at import time)."""
    saved = sys.stdout
    sys.stdout = io.StringIO()
    try:
        yield
    finally:
        sys.stdout = saved


def import_prometheus():
    """Import and return the Prometheus modules used by this package.

    Returns:
        tuple: ``(gasProperties, celestialBodies, geometryHandler, constants)``
        modules.
    """
    import Prometheus.pythonScripts.gasProperties as gasprop
    import Prometheus.pythonScripts.celestialBodies as bodies
    import Prometheus.pythonScripts.geometryHandler as geom
    import Prometheus.pythonScripts.constants as const
    return gasprop, bodies, geom, const
