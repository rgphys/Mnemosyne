"""Convenience builders for Prometheus wavelength/spatial grids and moons.

These wrap the Prometheus ``WavelengthGrid`` / ``geometryHandler.Grid`` /
``celestialBodies.Moon`` constructors with research-sensible defaults (Na-D
doublet window, stellar-radius integration limit) so a transit can be set up in
a couple of lines.  All lengths are CGS.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from ._bootstrap import import_prometheus

_gasprop, _bodies, _geom, _const = import_prometheus()


def find_planet(name: str) -> Any:
    """Look up a Prometheus ``Planet`` by name (from ``Resources/planets.csv``)."""
    planet = _bodies.AvailablePlanets().findPlanet(name)
    if planet is None:
        raise ValueError(f"Planet {name!r} not found in Prometheus resources.")
    return planet


def na_d_grid(lower_ang: float = 5880.0, upper_ang: float = 5910.0,
              widthHighRes: float = 4e-8, resolutionLow: float = 3e-9,
              resolutionHigh: float = 2e-10) -> Any:
    """A ``WavelengthGrid`` spanning the Na-D doublet (defaults 5880–5910 Å).

    Args:
        lower_ang, upper_ang: Wavelength bounds [Å].
        widthHighRes: High-res half-window around lines [cm].
        resolutionLow, resolutionHigh: Pixel sizes outside / inside the
            high-res window [cm].
    """
    return _gasprop.WavelengthGrid(
        lower_w=lower_ang * 1e-8, upper_w=upper_ang * 1e-8,
        widthHighRes=widthHighRes, resolutionLow=resolutionLow,
        resolutionHigh=resolutionHigh)


def line_grid(center_ang: float, half_window_ang: float = 15.0,
              **kwargs) -> Any:
    """A ``WavelengthGrid`` centred on an arbitrary line [Å]."""
    return na_d_grid(lower_ang=center_ang - half_window_ang,
                     upper_ang=center_ang + half_window_ang, **kwargs)


def spatial_grid(planet: Any, x_border_Rp: float = 12.0, x_steps: int = 25,
                 rho_steps: int = 60, phi_steps: int = 30,
                 orbphase_window: float = 0.0, orbphase_steps: int = 1,
                 rho_border: Optional[float] = None) -> Any:
    """A ``geometryHandler.Grid`` with defaults tuned for an extended exosphere.

    Args:
        planet: Prometheus ``Planet`` object.
        x_border_Rp: Half-length of the LOS chord, in planet radii.
        x_steps, rho_steps, phi_steps: Grid resolution.
        orbphase_window: Half-window of orbital phase [rad] (0 → single phase
            at mid-transit; >0 → a lightcurve).
        orbphase_steps: Number of orbital-phase samples.
        rho_border: Sky-plane integration radius [cm] (default: stellar radius;
            depth normalization is by the stellar disk — do not shrink it).
    """
    rho_border = planet.hostStar.R if rho_border is None else rho_border
    return _geom.Grid(
        x_midpoint=planet.a, x_border=x_border_Rp * planet.R, x_steps=x_steps,
        rho_border=rho_border, rho_steps=rho_steps, phi_steps=phi_steps,
        orbphase_border=orbphase_window, orbphase_steps=orbphase_steps)


def orbphase_window_from_hours(planet: Any, half_window_hours: float) -> float:
    """Convert a ±time half-window [hours] to an orbital-phase half-window [rad]."""
    period_hours = planet.orbitalPeriod * 24.0
    return (half_window_hours / period_hours) * 2.0 * np.pi


def make_moon(planet: Any, a_over_Rp: float = 1.7,
              R: Optional[float] = None,
              midTransitOrbphase: float = 0.375 * 2.0 * np.pi) -> Any:
    """A Prometheus ``Moon`` orbiting ``planet``.

    Args:
        planet: Host ``Planet`` object.
        a_over_Rp: Moon semi-major axis in planet radii (default 1.7).
        R: Moon radius [cm] (default Io radius).
        midTransitOrbphase: Moon orbital phase at mid-transit [rad].
    """
    R = _const.R_Io if R is None else R
    return _bodies.Moon(midTransitOrbphase=midTransitOrbphase, R=R,
                        a=a_over_Rp * planet.R, hostPlanet=planet)
