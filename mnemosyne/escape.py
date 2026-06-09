"""Escape / source models — the composable physics layer.

Each model binds dishoom energetics to a Prometheus body (a ``Planet`` or a
``Moon``) and exposes:

* :meth:`~EscapeModel.mass_loss_rate` — the dishoom source rate [g/s];
* :meth:`~EscapeModel.particle_number` — the steady neutral inventory ``N``
  (mass-loss rate × ionization lifetime / particle mass), for number-normalized
  exospheres;
* :meth:`~EscapeModel.build_scenario` — a ready-to-run Prometheus scenario
  (constituent + line lookup attached), delegating to :mod:`.scenarios`.

Two families map onto two Prometheus scenario kinds:

* **Moon-outgassing** (:class:`ThermalSublimation`, :class:`TidalHeating`,
  :class:`SurfaceJeansEscape`) → ``MoonExosphere`` normalized by ``N``.
* **Planetary wind** (:class:`EnergyLimitedEscape`) → ``RadialWindExosphere``
  normalized by ``Mdot`` through mass continuity.

You can use a model purely for its numbers (just call ``mass_loss_rate``) or
hand it to :class:`~mnemosyne.transit.EnergeticTransit`.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Optional

from . import dishoom_adapter as dishoom
from . import scenarios

CONST = dishoom.CONST


class EscapeModel(ABC):
    """Abstract escape/source model bound to a Prometheus body.

    Attributes:
        planet: The host ``Planet`` (always present, sets the orbital distance).
        species: Absorbing species key (e.g. ``'NaI'``).
        particle_mass: Mass of one tracer particle [g].
    """

    scenario_kind: str = "exosphere"  # 'exosphere' or 'wind'

    def __init__(self, planet: Any, species: str = "NaI",
                 particle_mass: Optional[float] = None):
        self.planet = planet
        self.species = species
        self.particle_mass = (CONST.m_Na if particle_mass is None
                              else particle_mass)

    # ── physics ────────────────────────────────────────────────────────────
    @abstractmethod
    def mass_loss_rate(self) -> float:
        """Source / mass-loss rate of the tracer [g/s]."""

    @property
    def orbital_distance_AU(self) -> float:
        """Planet's orbital distance [AU]."""
        return self.planet.a / CONST.au2cm

    def ionization_lifetime(self) -> float:
        """Photoionization lifetime of the tracer [s] at the planet's orbit."""
        return dishoom.ionization_lifetime(self.orbital_distance_AU)

    def particle_number(self, lifetime_s: Optional[float] = None) -> float:
        """Steady-state neutral particle inventory ``N``.

        Args:
            lifetime_s: Loss lifetime [s]; defaults to the photoionization
                lifetime at the planet's orbital distance.
        """
        if lifetime_s is None:
            lifetime_s = self.ionization_lifetime()
        return dishoom.mdot_to_particle_number(
            self.mass_loss_rate(), lifetime_s, self.particle_mass)

    def summary(self) -> dict:
        """A dict of the key derived quantities (for logging/printing)."""
        return {
            "model": type(self).__name__,
            "species": self.species,
            "a_AU": self.orbital_distance_AU,
            "mass_loss_rate_gs": self.mass_loss_rate(),
            "tau_ion_s": self.ionization_lifetime(),
            "N_particles": self.particle_number(),
        }

    # ── coupling to Prometheus ───────────────────────────────────────────────
    @abstractmethod
    def build_scenario(self, w_grid: Any, sigma_v: float, **kwargs) -> Any:
        """Build the Prometheus scenario this model maps onto."""


# ── Moon-outgassing models → MoonExosphere ──────────────────────────────────
class _MoonOutgassingModel(EscapeModel):
    """Common machinery for moon-sourced, number-normalized exospheres."""

    scenario_kind = "exosphere"

    def __init__(self, moon: Any, q: float = 3.34, species: str = "NaI",
                 particle_mass: Optional[float] = None):
        super().__init__(moon.hostPlanet, species, particle_mass)
        self.moon = moon
        self.q = q

    def build_scenario(self, w_grid: Any, sigma_v: float,
                       lifetime_s: Optional[float] = None, **kwargs) -> Any:
        N = self.particle_number(lifetime_s)
        return scenarios.moon_exosphere_scenario(
            N=N, q=self.q, moon=self.moon, species=self.species,
            sigma_v=sigma_v, w_grid=w_grid)


class ThermalSublimation(_MoonOutgassingModel):
    """Thermal (Hertz-Knudsen) sublimation outgassing from a moon.

    ``Ṁ = X · P_vap(A,B,T₀) · √(m/2πkT₀) · 4πR²`` (van Lieshout vapor pressure),
    via :func:`dishoom_adapter.mdot_thermal_sublimation`.

    Args:
        moon: Prometheus ``Moon`` object.
        T_surface: Peak sub-stellar surface temperature [K].
        xi: Tracer (volcanic Na) mass fraction (default dishoom volcanic value).
        q: Density power-law index for the exosphere (default 3.34).
        A, B, mu_mineral: Mineral coefficients (default MgSiO3).
    """

    def __init__(self, moon: Any, T_surface: float,
                 xi: Optional[float] = None, q: float = 3.34,
                 species: str = "NaI", particle_mass: Optional[float] = None,
                 A: Optional[float] = None, B: Optional[float] = None,
                 mu_mineral: Optional[float] = None):
        super().__init__(moon, q, species, particle_mass)
        self.T_surface = T_surface
        self.xi = CONST.X_Na_volcanic if xi is None else xi
        self.A, self.B, self.mu_mineral = A, B, mu_mineral

    def mass_loss_rate(self) -> float:
        return dishoom.mdot_thermal_sublimation(
            self.T_surface, self.moon.R, self.xi,
            A=self.A, B=self.B, mu_mineral=self.mu_mineral)


class TidalHeating(_MoonOutgassingModel):
    """Tidal-heating driven outgassing (optionally plus thermal).

    Uses dishoom ``Mdot0_tidal`` (or ``Mdot0_tot`` when ``include_thermal``),
    where solid-body tidal power ``Ė(P)`` is converted to a source rate.

    Args:
        moon: Prometheus ``Moon`` object.
        P_days: Moon's orbital period [days] (sets the tidal heating).
        tidal_eff: Fraction of tidal power driving escape.
        xi: Tracer mass fraction.
        include_thermal: Add the thermal sublimation term (needs ``T_eq``).
        T_eq: Equilibrium temperature [K] (required if ``include_thermal``).
    """

    def __init__(self, moon: Any, P_days: float, tidal_eff: float = 0.4,
                 xi: Optional[float] = None, include_thermal: bool = False,
                 T_eq: Optional[float] = None, moon_mass: Optional[float] = None,
                 q: float = 3.34,
                 species: str = "NaI", particle_mass: Optional[float] = None,
                 A: Optional[float] = None, B: Optional[float] = None,
                 mu_mineral: Optional[float] = None):
        super().__init__(moon, q, species, particle_mass)
        self.P_days = P_days
        self.tidal_eff = tidal_eff
        self.xi = CONST.X_Na_volcanic if xi is None else xi
        self.include_thermal = include_thermal
        self.T_eq = T_eq
        # Prometheus Moon carries no mass; default to Io (dishoom's reference).
        self.moon_mass = CONST.m_Io if moon_mass is None else moon_mass
        self.A, self.B, self.mu_mineral = A, B, mu_mineral

    def mass_loss_rate(self) -> float:
        if self.include_thermal:
            if self.T_eq is None:
                raise ValueError("include_thermal=True requires T_eq [K].")
            return dishoom.mdot_thermal_plus_tidal(
                self.T_eq, self.moon.R, self.moon_mass, self.P_days, self.xi,
                self.tidal_eff, A=self.A, B=self.B, mu_mineral=self.mu_mineral)
        return dishoom.mdot_tidal(self.moon_mass, self.moon.R, self.P_days,
                                  self.tidal_eff, self.xi)


class SurfaceJeansEscape(_MoonOutgassingModel):
    """Surface-Jeans escape (Volkov+ 2011 fit) for Na from a moon.

    Wraps ``dMdt_SJeansMdot0_tot_Na``: the total source rate multiplied by the
    Jeans escape fraction ``R(λ,Kn)·(1+λ)e^{-λ}``.

    Args:
        moon: Prometheus ``Moon`` object.
        T_eq: Equilibrium temperature [K].
        P_days: Moon's orbital period [days].
        xi: Tracer (Na) mass fraction.
    """

    def __init__(self, moon: Any, T_eq: float, P_days: float,
                 xi: Optional[float] = None, moon_mass: Optional[float] = None,
                 q: float = 3.34,
                 species: str = "NaI", particle_mass: Optional[float] = None,
                 A: Optional[float] = None, B: Optional[float] = None,
                 mu_mineral: Optional[float] = None):
        super().__init__(moon, q, species, particle_mass)
        self.T_eq = T_eq
        self.P_days = P_days
        self.xi = CONST.X_Na_volcanic if xi is None else xi
        # Prometheus Moon carries no mass; default to Io (dishoom's reference).
        self.moon_mass = CONST.m_Io if moon_mass is None else moon_mass
        self.A, self.B, self.mu_mineral = A, B, mu_mineral

    def mass_loss_rate(self) -> float:
        return dishoom.mdot_surface_jeans_Na(
            self.T_eq, self.moon.R, self.moon_mass, self.P_days, self.xi,
            A=self.A, B=self.B, mu_mineral=self.mu_mineral,
            m_volatile=self.particle_mass)

    def exobase_radius(self, n0: float,
                       lifetime_s: Optional[float] = None) -> float:
        """Hydrostatic exobase radius [cm] (``Rx_hydrostatic``, Io-scaled)."""
        if lifetime_s is None:
            lifetime_s = self.ionization_lifetime()
        return dishoom.exobase_radius_hydrostatic(
            self.T_eq, self.P_days, self.xi, n0, lifetime_s,
            m_volatile=self.particle_mass)


# ── Planetary wind model → RadialWindExosphere ──────────────────────────────
class EnergyLimitedEscape(EscapeModel):
    """Energy-limited (XUV-driven) planetary wind.

    dishoom ``dMdt_ELescape`` sets the tracer mass-loss rate; Prometheus'
    ``RadialWindExosphere`` turns it into ``n(r) = Ṁ/(4π r² v μ)`` along the
    line of sight, using either the exact isothermal **Parker wind**
    (``wind_model='parker'``) or a beta-law.

    Provide the source either as an XUV **luminosity** ``L_xuv`` [erg/s] or a
    local XUV **flux** ``F_xuv`` [erg/s/cm²] (with ``R_abs``).

    Args:
        planet: Prometheus ``Planet`` object (the escaping body).
        L_xuv: Host XUV luminosity [erg/s] (mutually exclusive with ``F_xuv``).
        F_xuv: XUV flux at the planet [erg/s/cm²].
        R_abs: XUV absorption radius [cm] (defaults to ``planet.R``).
        eta: Heating efficiency (default 0.1).
        xi: Tracer mass fraction of the species in the wind (default 1.0).
        T_wind: Isothermal wind temperature [K] (Parker dynamics).
        wind_model: ``'parker'`` (default) or ``'beta'``.
        wind_mu: Mean particle mass of the *bulk* wind [g] that sets the Parker
            dynamics (default: H/He, 1.3 amu). The tracer's ``particle_mass``
            sets only the density normalization.
        v_terminal, beta, v_base: beta-law knobs (only for ``wind_model='beta'``).
    """

    scenario_kind = "wind"

    def __init__(self, planet: Any, L_xuv: Optional[float] = None,
                 F_xuv: Optional[float] = None, R_abs: Optional[float] = None,
                 eta: float = 0.1, xi: float = 1.0,
                 T_wind: float = 1e4, wind_model: str = "parker",
                 wind_mu: Optional[float] = None,
                 v_terminal: Optional[float] = None, beta: float = 1.0,
                 v_base: Optional[float] = None,
                 species: str = "NaI", particle_mass: Optional[float] = None):
        super().__init__(planet, species, particle_mass)
        if (L_xuv is None) == (F_xuv is None):
            raise ValueError("Provide exactly one of L_xuv or F_xuv.")
        self.L_xuv = L_xuv
        self.F_xuv = F_xuv
        self.R_abs = planet.R if R_abs is None else R_abs
        self.eta = eta
        self.xi = xi
        self.T_wind = T_wind
        self.wind_model = wind_model
        self.wind_mu = 1.3 * CONST.amu if wind_mu is None else wind_mu
        self.v_terminal = v_terminal
        self.beta = beta
        self.v_base = v_base

    def mass_loss_rate(self) -> float:
        if self.L_xuv is not None:
            mdot_bulk = dishoom.mdot_energy_limited(
                self.L_xuv, self.orbital_distance_AU, self.planet.R,
                self.planet.M)
            return self.xi * mdot_bulk
        return dishoom.mdot_energy_limited_xuv(
            self.F_xuv, self.planet.R, self.planet.M, self.R_abs,
            self.eta, self.xi)

    def build_scenario(self, w_grid: Any, sigma_v: float, **kwargs) -> Any:
        return scenarios.radial_wind_scenario(
            Mdot=self.mass_loss_rate(), planet=self.planet,
            species=self.species, sigma_v=sigma_v, w_grid=w_grid,
            mu=self.particle_mass, wind_model=self.wind_model,
            T=self.T_wind, wind_mu=self.wind_mu, v_terminal=self.v_terminal,
            beta=self.beta, v_base=self.v_base)
