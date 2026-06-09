"""Builders that turn dishoom-derived quantities into Prometheus scenarios.

These are the low-level, composable seam: give them a particle number ``N`` or a
mass-loss rate ``Mdot`` and they return a fully-wired Prometheus scenario object
(constituent added, atomic line lookup attached to the wavelength grid) ready to
drop into a ``gasProperties.Atmosphere``.

Use them directly if you want to bypass the :mod:`.escape` objects, e.g.::

    scen = radial_wind_scenario(Mdot=2e7, planet=planet, species='NaI',
                                sigma_v=2e6, w_grid=w_grid)
    atmos = gasprop.Atmosphere([scen], hasOrbitalDopplerShift=True)
"""

from __future__ import annotations

from typing import Any, Optional

from ._bootstrap import import_prometheus

_gasprop, _bodies, _geom, _const = import_prometheus()


def _attach_constituent(scenario: Any, species: str, sigma_v: float,
                        w_grid: Any) -> Any:
    """Add an atomic constituent and its line lookup to a scenario."""
    scenario.addConstituent(species, sigma_v)
    scenario.constituents[-1].addLookupFunctionToConstituent(w_grid)
    return scenario


def moon_exosphere_scenario(N: float, q: float, moon: Any, species: str,
                            sigma_v: float, w_grid: Any) -> Any:
    """A ``MoonExosphere`` (power-law density, normalized by ``N``).

    Args:
        N: Total tracer particle number.
        q: Density power-law index (Io torus ≈ 3.34).
        moon: Prometheus ``Moon`` object.
        species: Species key (e.g. ``'NaI'``).
        sigma_v: Velocity dispersion [cm/s].
        w_grid: Prometheus ``WavelengthGrid``.
    """
    scen = _gasprop.MoonExosphere(N=N, q=q, moon=moon)
    return _attach_constituent(scen, species, sigma_v, w_grid)


def powerlaw_exosphere_scenario(N: float, q: float, planet: Any, species: str,
                                sigma_v: float, w_grid: Any) -> Any:
    """A planet-centred ``PowerLawExosphere`` normalized by ``N``."""
    scen = _gasprop.PowerLawExosphere(N=N, q=q, planet=planet)
    return _attach_constituent(scen, species, sigma_v, w_grid)


def radial_wind_scenario(Mdot: float, planet: Any, species: str,
                         sigma_v: float, w_grid: Any,
                         mu: Optional[float] = None,
                         wind_model: str = "parker", T: float = 1e4,
                         wind_mu: Optional[float] = None,
                         v_terminal: Optional[float] = None,
                         beta: float = 1.0,
                         v_base: Optional[float] = None,
                         r_inner: Optional[float] = None,
                         r_outer: Optional[float] = None) -> Any:
    """A ``RadialWindExosphere`` normalized by mass continuity from ``Mdot``.

    ``n(r) = Mdot / (4π r² v(r) μ)`` with either the exact isothermal Parker
    wind (``wind_model='parker'``) or a beta-law velocity profile.

    Args:
        Mdot: Tracer mass-loss rate [g/s].
        planet: Prometheus ``Planet`` object.
        species: Species key.
        sigma_v: Velocity dispersion [cm/s].
        w_grid: Prometheus ``WavelengthGrid``.
        mu: Tracer particle mass [g] (default Na atom mass).
        wind_model: ``'parker'`` or ``'beta'``.
        T: Isothermal wind temperature [K] (Parker dynamics).
        wind_mu: Bulk wind mean particle mass [g] for Parker dynamics.
        v_terminal, beta, v_base: beta-law parameters.
        r_inner, r_outer: Optional wind radial bounds [cm].
    """
    if mu is None:
        # Default the tracer particle mass to the species' own mass.
        mu = _const.AvailableSpecies().findSpecies(species).mass
    scen = _gasprop.RadialWindExosphere(
        Mdot=Mdot, mu=mu, v_terminal=v_terminal, beta=beta,
        r_inner=r_inner, r_outer=r_outer, v_base=v_base,
        wind_model=wind_model, T=T, planet=planet, wind_mu=wind_mu)
    return _attach_constituent(scen, species, sigma_v, w_grid)
