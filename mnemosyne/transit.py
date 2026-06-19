"""High-level Mnemosyne transit API.

:class:`EnergeticTransit` runs the full pipeline **in memory** — dishoom
energetics → Prometheus scenario → chord-summed transit — and returns a
:class:`TransitResult` carrying the transit-depth cube ``R(phase, λ)`` plus
convenience accessors for the transmission spectrum and the band-integrated
lightcurve.  No setup files, no output files.

The radiative-transfer orchestration (build the ``Atmosphere``/``Transit``, run
``sumOverChords``, package the cube) and the :class:`TransitResult` wrapper now
live in Prometheus (``gasProperties.run_transit`` / ``gasProperties.TransitResult``);
this class only contributes the dishoom-coupled scenario assembly and re-exports
``TransitResult`` for backward compatibility.

Example::

    from mnemosyne import EnergeticTransit, ThermalSublimation, grids

    planet = grids.find_planet('WASP-49b')
    moon   = grids.make_moon(planet)
    model  = ThermalSublimation(moon, T_surface=1500.0)
    result = EnergeticTransit(model).run()
    wav_ang, depth = result.wavelength_ang, result.transit_depth()
"""

from __future__ import annotations

from typing import Any, List, Optional

from . import grids as _grids
from ._bootstrap import import_prometheus

_gasprop, _bodies, _geom, _const = import_prometheus()

# TransitResult was upstreamed into Prometheus; re-export it so
# ``mnemosyne.TransitResult`` and ``mnemosyne.transit.TransitResult`` keep working.
TransitResult = _gasprop.TransitResult

# Na D2 / D1 rest wavelengths [Å], vacuum (matching the Prometheus LineList).
# Now sourced from Prometheus constants; kept here for backward compatibility.
NA_D2_ANG = _const.NA_D2_ANG
NA_D1_ANG = _const.NA_D1_ANG


class EnergeticTransit:
    """Run a dishoom→Prometheus transit in memory.

    Args:
        escape_model: An :class:`~mnemosyne.escape.EscapeModel`.
        w_grid: Prometheus ``WavelengthGrid`` (default: Na-D doublet window).
        s_grid: Prometheus ``geometryHandler.Grid`` (default: extended-exosphere
            grid at mid-transit, single orbital phase).
        sigma_v: Tracer velocity dispersion [cm/s] (default 2e6 = 20 km/s).
        hasOrbitalDopplerShift: Apply the orbital + wind Doppler shift.
        use_phoenix_star: Use the PHOENIX stellar spectrum (``addFstarFunction``).
            Set ``False`` for a flat star (much faster; fine for relative depths).
        extra_scenarios: Additional pre-built Prometheus scenarios to co-add.
    """

    def __init__(self, escape_model: Any,
                 w_grid: Optional[Any] = None, s_grid: Optional[Any] = None,
                 sigma_v: float = 2e6, hasOrbitalDopplerShift: bool = True,
                 use_phoenix_star: bool = True,
                 extra_scenarios: Optional[List[Any]] = None):
        self.escape_model = escape_model
        self.planet = escape_model.planet
        self.w_grid = w_grid if w_grid is not None else _grids.na_d_grid()
        self.s_grid = (s_grid if s_grid is not None
                       else _grids.spatial_grid(self.planet))
        self.sigma_v = sigma_v
        self.hasOrbitalDopplerShift = hasOrbitalDopplerShift
        self.use_phoenix_star = use_phoenix_star
        self.extra_scenarios = list(extra_scenarios or [])

    def build_scenarios(self) -> List[Any]:
        """Build the Prometheus scenario list (escape model + extras)."""
        scen = self.escape_model.build_scenario(self.w_grid, self.sigma_v)
        return [scen, *self.extra_scenarios]

    def run(self, max_memory_gb: float = 4.0) -> Any:
        """Execute the transit and return a :class:`TransitResult`.

        Delegates the radiative-transfer orchestration to
        ``Prometheus.pythonScripts.gasProperties.run_transit``.

        Args:
            max_memory_gb: Memory cap passed to ``Transit.sumOverChords``.
        """
        return _gasprop.run_transit(
            self.build_scenarios(), self.w_grid, self.s_grid,
            hasOrbitalDopplerShift=self.hasOrbitalDopplerShift,
            use_phoenix_star=self.use_phoenix_star,
            max_memory_gb=max_memory_gb)
