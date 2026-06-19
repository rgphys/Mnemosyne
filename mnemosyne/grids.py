"""Backward-compatible re-exports of grid/body/orbital-mechanics helpers.

These convenience builders — wavelength + spatial grids, planet/moon factories,
and the moon:planet orbital-mechanics relations — were **upstreamed into
Prometheus** because they are pure radiative-transfer / orbital-mechanics
helpers with no dishoom dependency.  They now live in:

* ``Prometheus.pythonScripts.gasProperties``  — :func:`na_d_grid`, :func:`line_grid`
* ``Prometheus.pythonScripts.geometryHandler`` — :func:`spatial_grid`,
  :func:`orbphase_window_from_hours`
* ``Prometheus.pythonScripts.celestialBodies`` — :func:`find_planet`,
  :func:`make_moon`, :func:`mean_motion_ratio`, :func:`optimal_midtransit_phase`,
  :func:`max_peak_shift_minutes`

This module re-exports them under their original ``mnemosyne.grids`` names so
existing callers keep working unchanged.  New code can import them straight from
Prometheus.
"""

from __future__ import annotations

from ._bootstrap import import_prometheus

_gasprop, _bodies, _geom, _const = import_prometheus()

# Wavelength grids (Prometheus gasProperties)
na_d_grid = _gasprop.na_d_grid
line_grid = _gasprop.line_grid

# Spatial / temporal grids (Prometheus geometryHandler)
spatial_grid = _geom.spatial_grid
orbphase_window_from_hours = _geom.orbphase_window_from_hours

# Bodies + moon orbital mechanics (Prometheus celestialBodies)
find_planet = _bodies.find_planet
make_moon = _bodies.make_moon
mean_motion_ratio = _bodies.mean_motion_ratio
optimal_midtransit_phase = _bodies.optimal_midtransit_phase
max_peak_shift_minutes = _bodies.max_peak_shift_minutes

__all__ = [
    "na_d_grid",
    "line_grid",
    "spatial_grid",
    "orbphase_window_from_hours",
    "find_planet",
    "make_moon",
    "mean_motion_ratio",
    "optimal_midtransit_phase",
    "max_peak_shift_minutes",
]
