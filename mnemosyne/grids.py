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


#  Optimal moon phase theory (see Test/midtransit_phase_proof.tex)
#
# The moon's sky-plane offset during transit is
#     y_m(theta) = a_p sin(theta) + a_m sin(theta_0 + N theta),
# with N the moon:planet mean-motion ratio.  Writing eps = a_m/a_p and
# expanding the lightcurve L = f(y_m/a_p) to first order in eps, the only
# time-antisymmetric (i.e. detectable against any symmetric bare-planet
# model) term is  eps sin(theta_0) f'(theta) cos(N theta), so every
# asymmetry observable is proportional to sin(theta_0) and maximised at
# quadrature.  The peak displacement admits an exact all-orders optimum at
# quadrature corrected by the moon's own motion during the displacement.


def mean_motion_ratio(planet: Any, a_over_Rp: float = 1.7) -> float:
    """Moon:planet mean-motion ratio N = sqrt(a_p^3 M_p / (a_m^3 M_star))."""
    a_m = a_over_Rp * planet.R
    return float(np.sqrt((planet.a**3 * planet.M) /
                         (a_m**3 * planet.hostStar.M)))


def optimal_midtransit_phase(planet: Any, a_over_Rp: float = 1.7,
                             branch: str = 'late') -> float:
    """Moon phase at mid-transit maximising the lightcurve peak shift [rad].

    Exact closed form (all orders in a_m/a_p, any monotone cloud profile):
    the moon must reach maximum sky-plane elongation exactly as its cloud
    crosses the stellar disk centre, which displaces the absorption peak by
    the maximum possible +/- arcsin(a_m/a_p) of planet phase.

        theta0* = 3*pi/2 - N*arcsin(a_m/a_p)   (branch='late',  peak after
                                                mid-transit, trailing moon)
        theta0* =   pi/2 + N*arcsin(a_m/a_p)   (branch='early', peak before
                                                mid-transit, leading moon)

    Both branches also sit on the flat |sin(theta0)| plateau of the
    detectability (antisymmetric-RMS) curve, within <1% of its maximum.

    Args:
        planet: Host ``Planet`` object.
        a_over_Rp: Moon semi-major axis in planet radii.
        branch: 'late' or 'early' peak displacement.
    """
    N = mean_motion_ratio(planet, a_over_Rp)
    delta = np.arcsin(a_over_Rp * planet.R / planet.a)
    if branch == 'late':
        return float(3 * np.pi / 2 - N * delta)
    if branch == 'early':
        return float(np.pi / 2 + N * delta)
    raise ValueError(f"branch must be 'late' or 'early', got {branch!r}")


def max_peak_shift_minutes(planet: Any, a_over_Rp: float = 1.7) -> float:
    """Maximum achievable lightcurve peak displacement [minutes].

    Delta_t = (T_p / 2 pi) arcsin(a_m/a_p): the time the planet takes to
    traverse one moon-orbit radius in the sky.  This bound is attained at
    ``optimal_midtransit_phase``.
    """
    delta = np.arcsin(a_over_Rp * planet.R / planet.a)
    return float(delta / (2 * np.pi) * planet.orbitalPeriod * 24.0 * 60.0)


def moon_phase_asymmetry_sweep(planet: Any, a_over_Rp: float = 1.7,
                               window_hours: float = 3.0,
                               cloud_sigma_factor: float = 3.0,
                               n_theta0: int = 3601) -> dict:
    """Geometric sweep of lightcurve-asymmetry metrics over theta0.

    Uses the Gaussian-kernel transit model L = exp(-y_m^2 / 2 sigma_y^2)
    with sigma_y = ``cloud_sigma_factor`` * a_m (no radiative transfer), so
    the full sweep is instantaneous.  Returns numeric metrics alongside the
    first-order closed-form predictions from the proof document.

    Returns dict with keys:
        theta0_deg, peak_shift_min, centroid_shift_min, antisym_rms,
        skewness  (numeric arrays over theta0), and
        pred_centroid_min, pred_antisym_rms, pred_skewness (closed-form
        arrays), theta0_opt_late_deg, theta0_opt_early_deg,
        max_peak_shift_min (scalars).
    """
    a_m = a_over_Rp * planet.R
    eps = a_m / planet.a
    N = mean_motion_ratio(planet, a_over_Rp)
    s = cloud_sigma_factor * eps
    h = N**2 * s**2 / 2
    P_hr = planet.orbitalPeriod * 24.0
    min_per_rad = P_hr * 60.0 / (2 * np.pi)
    window_ph = (window_hours / P_hr) * 2 * np.pi
    T = 2 * window_ph

    theta = np.linspace(-window_ph, window_ph, 8001)
    theta0 = np.radians(np.linspace(0, 360, n_theta0, endpoint=False))
    out = {k: np.zeros(n_theta0) for k in
           ('peak_shift_min', 'centroid_shift_min', 'antisym_rms', 'skewness')}
    for i, t0 in enumerate(theta0):
        y = np.sin(theta) + eps * np.sin(t0 + N * theta)
        lc = np.exp(-0.5 * (y / s)**2)
        out['peak_shift_min'][i] = theta[np.argmax(lc)] * min_per_rad
        W = lc.sum()
        mu = (theta * lc).sum() / W
        out['centroid_shift_min'][i] = mu * min_per_rad
        var = ((theta - mu)**2 * lc).sum() / W
        out['skewness'][i] = ((theta - mu)**3 * lc).sum() / W / var**1.5
        A = 0.5 * (lc - np.interp(-theta, theta, lc))
        out['antisym_rms'][i] = np.sqrt(np.mean(A**2))

    t0_deg = np.degrees(theta0)
    sin_t0 = np.sin(theta0)
    out.update(
        theta0_deg=t0_deg,
        pred_centroid_min=-eps * (1 - N**2 * s**2) * np.exp(-h)
                          * sin_t0 * min_per_rad,
        pred_antisym_rms=eps * np.abs(sin_t0) * np.sqrt(
            np.sqrt(np.pi) / (4 * s * T)
            * (1 + (1 - 2 * N**2 * s**2) * np.exp(-N**2 * s**2))),
        pred_skewness=(eps / s) * np.exp(-h) * N**2 * s**2
                      * (3 - N**2 * s**2) * sin_t0,
        theta0_opt_late_deg=np.degrees(optimal_midtransit_phase(
            planet, a_over_Rp, 'late')),
        theta0_opt_early_deg=np.degrees(optimal_midtransit_phase(
            planet, a_over_Rp, 'early')),
        max_peak_shift_min=max_peak_shift_minutes(planet, a_over_Rp),
    )
    return out
