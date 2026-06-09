"""High-level Mnemosyne transit API.

:class:`EnergeticTransit` runs the full pipeline **in memory** — dishoom
energetics → Prometheus scenario → chord-summed transit — and returns a
:class:`TransitResult` carrying the transit-depth cube ``R(phase, λ)`` plus
convenience accessors for the transmission spectrum and the band-integrated
lightcurve.  No setup files, no output files.

Example::

    from mnemosyne import EnergeticTransit, ThermalSublimation, grids

    planet = grids.find_planet('WASP-49b')
    moon   = grids.make_moon(planet)
    model  = ThermalSublimation(moon, T_surface=1500.0)
    result = EnergeticTransit(model).run()
    wav_ang, depth = result.wavelength_ang, result.transit_depth()
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Optional

import numpy as np

from . import grids as _grids
from ._bootstrap import import_prometheus

_gasprop = import_prometheus()[0]

# Na D2 / D1 rest wavelengths [Å] for default bandpasses.
NA_D2_ANG = 5889.95
NA_D1_ANG = 5895.92


@dataclass
class TransitResult:
    """Result of a transit calculation.

    Attributes:
        wavelength_cm: Wavelength grid [cm], shape ``(n_wav,)``.
        R_2D: Transit depth cube ``R(phase, λ)``, shape ``(n_phase, n_wav)``.
        orbphase: Orbital-phase axis [rad], shape ``(n_phase,)``.
        planet: The Prometheus ``Planet`` used.
    """

    wavelength_cm: np.ndarray
    R_2D: np.ndarray
    orbphase: np.ndarray
    planet: Any

    # ── axis conversions ─────────────────────────────────────────────────────
    @property
    def wavelength_ang(self) -> np.ndarray:
        """Wavelength grid [Å]."""
        return self.wavelength_cm * 1e8

    @property
    def wavelength_um(self) -> np.ndarray:
        """Wavelength grid [µm]."""
        return self.wavelength_cm * 1e4

    # ── transmission spectrum ────────────────────────────────────────────────
    def spectrum(self) -> np.ndarray:
        """Phase-collapsed transit depth ``R(λ)`` (median over orbital phase)."""
        return np.median(self.R_2D, axis=0)

    def spectrum_normalized(self) -> np.ndarray:
        """Transmission spectrum normalized to its continuum (max → 1)."""
        spec = self.spectrum()
        return spec / spec.max()

    # ── derived line metrics ─────────────────────────────────────────────────
    def _masks(self, line_window_ang, continuum_exclude_ang):
        wav = self.wavelength_ang
        lo, hi = line_window_ang
        line_mask = (wav >= lo) & (wav <= hi)
        if continuum_exclude_ang is None:
            cont_mask = ~line_mask
        else:
            clo, chi = continuum_exclude_ang
            cont_mask = (wav < clo) | (wav > chi)
        return line_mask, cont_mask

    def transit_depth(self, line_window_ang=(NA_D2_ANG - 4.0, NA_D2_ANG + 4.0),
                      continuum_exclude_ang=(5884.0, 5902.0),
                      mode: str = "peak") -> float:
        """Excess absorption depth in a line window, vs a clean continuum.

        Args:
            line_window_ang: ``(lo, hi)`` line bandpass [Å].
            continuum_exclude_ang: ``(lo, hi)`` region to exclude from the
                continuum estimate [Å]; ``None`` uses everything outside the
                line window.
            mode: ``'peak'`` (deepest pixel, matches the Doppler-shifted moon
                cloud convention) or ``'mean'`` (band-averaged).

        Returns:
            Excess absorption as a fraction (multiply by 100 for percent).
        """
        spec = self.spectrum()
        line_mask, cont_mask = self._masks(line_window_ang,
                                           continuum_exclude_ang)
        R_cont = np.median(spec[cont_mask])
        if mode == "peak":
            return float(1.0 - spec[line_mask].min() / R_cont)
        elif mode == "mean":
            return float(1.0 - spec[line_mask].mean() / R_cont)
        raise ValueError(f"Unknown mode {mode!r}; use 'peak' or 'mean'.")

    def lightcurve(self, line_window_ang=(NA_D2_ANG - 0.375, NA_D2_ANG + 0.375),
                   continuum_exclude_ang=(5884.0, 5902.0), mode='mean') -> np.ndarray:
        """Band-integrated lightcurve ``L(phase)`` = line / continuum.

        Values < 1 mark net excess absorption at that orbital phase.  Requires a
        grid built with ``orbphase_steps > 1``.

        Args:
            line_window_ang: ``(lo, hi)`` line bandpass [Å].
            continuum_exclude_ang: continuum exclusion region [Å].
            mode: 'mean' for band-integrated flux, 'peak' for max line depth.
        """
        line_mask, cont_mask = self._masks(line_window_ang,
                                           continuum_exclude_ang)
        cont_per_phase = np.mean(self.R_2D[:, cont_mask], axis=1)
        
        if mode == 'mean':
            line_per_phase = np.mean(self.R_2D[:, line_mask], axis=1)
        elif mode == 'peak':
            line_per_phase = np.min(self.R_2D[:, line_mask], axis=1)
        else:
            raise ValueError(f"Unknown mode {mode!r}; use 'peak' or 'mean'.")
            
        return line_per_phase / cont_per_phase


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

    def run(self, max_memory_gb: float = 4.0) -> TransitResult:
        """Execute the transit and return a :class:`TransitResult`.

        Args:
            max_memory_gb: Memory cap passed to ``Transit.sumOverChords``.
        """
        atmos = _gasprop.Atmosphere(
            self.build_scenarios(),
            hasOrbitalDopplerShift=self.hasOrbitalDopplerShift)
        sim = _gasprop.Transit(atmos, self.w_grid, self.s_grid)
        sim.addWavelength()
        if self.use_phoenix_star:
            self.planet.hostStar.addFstarFunction(sim.wavelength)
        R_2D = sim.sumOverChords(max_memory_gb=max_memory_gb)

        orbphase = np.linspace(-self.s_grid.orbphase_border,
                               self.s_grid.orbphase_border, R_2D.shape[0])
        return TransitResult(wavelength_cm=np.asarray(sim.wavelength),
                             R_2D=np.asarray(R_2D), orbphase=orbphase,
                             planet=self.planet)
