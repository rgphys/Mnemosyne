"""Energetic Prometheus — in-memory coupling of dishoom escape physics and the
Prometheus radiative-transfer engine.

Instead of handing JSON setup files between the two codebases, this package lets
you compute a moon/planet mass-loss rate with **dishoom** and feed it straight
into a **Prometheus** transit, returning the transmission spectrum and lightcurve
as arrays.

Two layers:

* **Composable** — :mod:`escape` models (``ThermalSublimation``, ``TidalHeating``,
  ``SurfaceJeansEscape``, ``EnergyLimitedEscape``) expose dishoom quantities and
  build Prometheus scenarios; :mod:`scenarios` builders accept raw ``N``/``Mdot``;
  :mod:`grids` builds wavelength/spatial grids and moons.
* **Headline** — :class:`transit.EnergeticTransit` runs the whole pipeline and
  returns a :class:`transit.TransitResult` with ``.spectrum()`` / ``.lightcurve()``.

Quick start::

    from energetic_prometheus import EnergeticTransit, ThermalSublimation, grids

    planet = grids.find_planet('WASP-49b')
    moon   = grids.make_moon(planet)
    result = EnergeticTransit(ThermalSublimation(moon, T_surface=1500.0)).run()
    print(result.transit_depth() * 100, '% Na D2')
"""

from __future__ import annotations

from . import dishoom_adapter, escape, grids, scenarios, transit
from .dishoom_adapter import CONST
from .escape import (
    EnergyLimitedEscape,
    EscapeModel,
    SurfaceJeansEscape,
    ThermalSublimation,
    TidalHeating,
)
from .scenarios import (
    moon_exosphere_scenario,
    powerlaw_exosphere_scenario,
    radial_wind_scenario,
)
from .transit import EnergeticTransit, TransitResult

__all__ = [
    # headline
    "EnergeticTransit",
    "TransitResult",
    # escape models
    "EscapeModel",
    "ThermalSublimation",
    "TidalHeating",
    "SurfaceJeansEscape",
    "EnergyLimitedEscape",
    # scenario builders
    "moon_exosphere_scenario",
    "powerlaw_exosphere_scenario",
    "radial_wind_scenario",
    # constants + submodules
    "CONST",
    "grids",
    "scenarios",
    "escape",
    "transit",
    "dishoom_adapter",
]
