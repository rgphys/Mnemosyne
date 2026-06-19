"""Adapters over the scenario builders now living in Prometheus.

The builders that turn a particle number ``N`` or a mass-loss rate ``Mdot`` into
a fully-wired Prometheus density distribution (constituent + atomic line lookup
attached) were **upstreamed into Prometheus** — they construct Prometheus
objects and carry no dishoom dependency — and now live in
``Prometheus.pythonScripts.gasProperties``.

These thin wrappers re-expose them under the original ``mnemosyne.scenarios``
names and keyword (``w_grid``), so existing callers and :mod:`.escape` keep
working unchanged.  The only difference from calling Prometheus directly is the
``w_grid`` keyword, which Prometheus spells ``wavelengthGrid``.

    scen = radial_wind_scenario(Mdot=2e7, planet=planet, species='NaI',
                                sigma_v=2e6, w_grid=w_grid)
    atmos = gasprop.Atmosphere([scen], hasOrbitalDopplerShift=True)
"""

from __future__ import annotations

from typing import Any, Optional

from ._bootstrap import import_prometheus

_gasprop = import_prometheus()[0]


def moon_exosphere_scenario(N: float, q: float, moon: Any, species: str,
                            sigma_v: float, w_grid: Any) -> Any:
    """A ``MoonExosphere`` (power-law density, normalized by ``N``).

    See :func:`Prometheus.pythonScripts.gasProperties.moon_exosphere_scenario`.
    """
    return _gasprop.moon_exosphere_scenario(
        N=N, q=q, moon=moon, species=species, sigma_v=sigma_v,
        wavelengthGrid=w_grid)


def powerlaw_exosphere_scenario(N: float, q: float, planet: Any, species: str,
                                sigma_v: float, w_grid: Any) -> Any:
    """A planet-centred ``PowerLawExosphere`` normalized by ``N``.

    See
    :func:`Prometheus.pythonScripts.gasProperties.powerlaw_exosphere_scenario`.
    """
    return _gasprop.powerlaw_exosphere_scenario(
        N=N, q=q, planet=planet, species=species, sigma_v=sigma_v,
        wavelengthGrid=w_grid)


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

    See :func:`Prometheus.pythonScripts.gasProperties.radial_wind_scenario`.
    """
    return _gasprop.radial_wind_scenario(
        Mdot=Mdot, planet=planet, species=species, sigma_v=sigma_v,
        wavelengthGrid=w_grid, mu=mu, wind_model=wind_model, T=T,
        wind_mu=wind_mu, v_terminal=v_terminal, beta=beta, v_base=v_base,
        r_inner=r_inner, r_outer=r_outer)
