"""Thin, typed adapter over dishoom's escape physics.

dishoom (A. V. Oza et al., built on A. Gebek's radiative-transfer code) remains
the **canonical** source of the energetics formulas — this module does not
re-derive them.  It imports the real dishoom functions once (with stdout
silenced), exposes them under documented names with explicit cgs units, and
collects the handful of mineral/atomic constants the escape models need.

All quantities are CGS.  ``Mdot`` values are mass-loss rates in **g/s**.

References (as cited inside dishoom):
    * Energy-limited escape — e.g. Watson+ 1981; Erkaev+ 2007; Oza+ 2019.
    * Thermal (Hertz-Knudsen) sublimation flux — van Lieshout+ 2014.
    * Surface-Jeans escape fit ``R(λ, Kn)`` — Volkov+ 2011.
    * Tidal heating ``Ė(P)`` — Cassidy+ 2009.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._bootstrap import silence_stdout
from . import constants as _c

#  Import dishoom once, quietly 
with silence_stdout():
    import moon_functions_5 as _moon
    import functions4 as _f4


#  Constants hub (CGS) 
@dataclass(frozen=True)
class DishoomConstants:
    """Frequently-needed constants, sourced from :mod:`.constants`.

    Provides a single ``CONST`` object so the rest of Mnemosyne
    never has to import from dishoom directly for constant values.
    """

    amu: float = _c.amu
    k_B: float = _c.k_B
    G: float = _c.G
    m_Na: float = _c.m_Na
    m_K: float = _c.m_K
    m_Fe: float = _c.m_Fe
    m_Mg: float = _c.m_Mg
    k_hv_Na: float = _c.k_hv_Na
    au2cm: float = _c.au2cm
    r_Io: float = _c.r_Io
    m_Io: float = _c.m_Io
    A_mgsio3: float = _c.MINERALS["MgSiO3"].A
    B_mgsio3: float = _c.MINERALS["MgSiO3"].B
    mu_mgsio3: float = _c.MINERALS["MgSiO3"].mu
    X_Na_volcanic: float = _c.X_Na_volcanic


CONST = DishoomConstants()

#  Mineral database — (A, B, mu) tuples for backward compatibility 
# The canonical database with references lives in :data:`constants.MINERALS`.
MINERAL_DB: dict[str, tuple[float, float, float]] = {
    name: (e.A, e.B, e.mu) for name, e in _c.MINERALS.items()
}


#  Energy-limited escape 
def mdot_energy_limited(L_xuv: float, a_AU: float, R_body: float,
                        M_body: float) -> float:
    """Energy-limited mass-loss rate [g/s] (efficiency η = 0.1).

    Wraps ``moon_functions_5.dMdt_ELescape``.

    Args:
        L_xuv: Host XUV luminosity [erg/s].
        a_AU: Body's orbital distance from the star [AU].
        R_body: Body radius [cm].
        M_body: Body mass [g].
    """
    return float(_moon.dMdt_ELescape(L_xuv, a_AU, R_body, M_body))


def mdot_energy_limited_xuv(F_xuv: float, R_body: float, M_body: float,
                            R_abs: float, eta: float = 0.1,
                            xi: float = 1.0) -> float:
    """Energy-limited mass-loss rate [g/s] from an explicit XUV flux.

    Wraps ``moon_functions_5.dMdt_ELescape_xuv``.

    Args:
        F_xuv: XUV flux at the body [erg/s/cm^2].
        R_body: Body radius [cm].
        M_body: Body mass [g].
        R_abs: Effective XUV absorption radius [cm].
        eta: Heating efficiency (default 0.1).
        xi: Tracer mass fraction of the species of interest (default 1.0).
    """
    return float(_moon.dMdt_ELescape_xuv(F_xuv, R_body, M_body, R_abs, eta, xi))


#  Thermal (Hertz-Knudsen) sublimation 
def mdot_thermal_sublimation(T_surface: float, R_body: float, xi: float,
                             A: float | None = None, B: float | None = None,
                             mu_mineral: float | None = None) -> float:
    """Thermal sublimation source rate [g/s] (no tidal term).

    Wraps ``moon_functions_5.Mdot0_thermal_notides``: a Hertz-Knudsen vapor flux
    over the full surface, scaled by the tracer fraction ``xi``.  Defaults to
    MgSiO3 (enstatite) coefficients.

    Args:
        T_surface: Sub-stellar surface temperature [K].
        R_body: Body radius [cm].
        xi: Tracer (e.g. volcanic Na) mass fraction.
        A, B: van Lieshout vapor-pressure coefficients (default MgSiO3).
        mu_mineral: Mineral molar mass [amu] (default MgSiO3).
    """
    A = CONST.A_mgsio3 if A is None else A
    B = CONST.B_mgsio3 if B is None else B
    mu_mineral = CONST.mu_mgsio3 if mu_mineral is None else mu_mineral
    return float(_moon.Mdot0_thermal_notides(A, B, T_surface, mu_mineral,
                                             R_body, xi))


def vapor_pressure(T: float, A: float | None = None,
                   B: float | None = None) -> float:
    """Saturation vapor pressure [dyne/cm^2] via ``functions4.Pvap_rocks``."""
    A = CONST.A_mgsio3 if A is None else A
    B = CONST.B_mgsio3 if B is None else B
    return float(_f4.Pvap_rocks(A, B, T))


#  Tidal-heating driven source 
def mdot_tidal(M_body: float, R_body: float, P_days: float,
               tidal_eff: float, xi: float) -> float:
    """Tidal-heating driven mass-loss rate [g/s].

    Wraps ``moon_functions_5.Mdot0_tidal`` (3-body solid-body tidal heating
    ``Ė(P)`` converted to a source rate via the gravitational binding energy).

    Args:
        M_body: Body mass [g].
        R_body: Body radius [cm].
        P_days: Orbital period of the body [days].
        tidal_eff: Fraction of tidal power going into escape.
        xi: Tracer mass fraction.
    """
    return float(_moon.Mdot0_tidal(M_body, R_body, P_days, tidal_eff, xi))


def mdot_thermal_plus_tidal(T_eq: float, R_body: float, M_body: float,
                            P_days: float, xi: float, tidal_eff: float,
                            A: float | None = None, B: float | None = None,
                            mu_mineral: float | None = None) -> float:
    """Combined thermal + tidal source rate [g/s] (``Mdot0_tot``)."""
    A = CONST.A_mgsio3 if A is None else A
    B = CONST.B_mgsio3 if B is None else B
    mu_mineral = CONST.mu_mgsio3 if mu_mineral is None else mu_mineral
    return float(_moon.Mdot0_tot(A, B, T_eq, mu_mineral, R_body, P_days,
                                 tidal_eff, M_body, xi))


#  Surface-Jeans escape + exobase radius 
def mdot_surface_jeans_Na(T_eq: float, R_body: float, M_body: float,
                          P_days: float, xi: float,
                          A: float | None = None, B: float | None = None,
                          mu_mineral: float | None = None,
                          m_volatile: float | None = None) -> float:
    """Surface-Jeans escape rate for Na [g/s].

    Wraps ``moon_functions_5.dMdt_SJeansMdot0_tot_Na`` — the total source rate
    multiplied by the Volkov+ 2011 Jeans-escape fraction ``R·(1+λ)e^{-λ}``.

    Args:
        T_eq: Equilibrium temperature [K].
        R_body: Body radius [cm].
        M_body: Body mass [g].
        P_days: Body orbital period [days].
        xi: Tracer (Na) mass fraction.
        A, B, mu_mineral: Mineral coefficients (default MgSiO3).
        m_volatile: Escaping-species mass [g] (default Na atom mass).
    """
    A = CONST.A_mgsio3 if A is None else A
    B = CONST.B_mgsio3 if B is None else B
    mu_mineral = CONST.mu_mgsio3 if mu_mineral is None else mu_mineral
    m_volatile = CONST.m_Na if m_volatile is None else m_volatile
    return float(_moon.dMdt_SJeansMdot0_tot_Na(A, B, T_eq, mu_mineral,
                                               m_volatile, R_body, M_body,
                                               P_days, xi))


def exobase_radius_hydrostatic(T_eq: float, P_days: float, xi: float,
                               n0: float, t_ion_s: float,
                               m_volatile: float | None = None,
                               frams: float = 1.0) -> float:
    """Hydrostatic exobase radius [cm] for an Io-like body.

    Wraps ``moon_functions_5.Rx_hydrostatic``.  Note this dishoom routine is
    hard-wired to Io's mass/radius internally; ``T_eq``/``n0``/``t_ion`` set the
    scale-height and column that locate the exobase.

    Args:
        T_eq: Equilibrium temperature [K].
        P_days: Body orbital period [days].
        xi: Na/volcanic gas fraction.
        n0: Base number density [cm^-3].
        t_ion_s: Ionization lifetime [s].
        m_volatile: Volatile mass [g] (default Na atom mass).
        frams: Ram-pressure scaling (default 1.0).
    """
    m_volatile = CONST.m_Na if m_volatile is None else m_volatile
    return float(_moon.Rx_hydrostatic(frams, T_eq, P_days, m_volatile,
                                      t_ion_s, xi, n0))


#  Helpers 
def ionization_lifetime(a_AU: float, k_hv_1AU: float | None = None) -> float:
    """Photoionization lifetime [s] at orbital distance ``a_AU``.

    ``τ_ion = 1 / (k_hv_1AU / a_AU^2)`` — the 1-AU rate diluted by the inverse
    square of the orbital distance.  Defaults to dishoom's Na rate.
    """
    k_hv_1AU = CONST.k_hv_Na if k_hv_1AU is None else k_hv_1AU
    return a_AU ** 2 / k_hv_1AU


def mdot_to_particle_number(mdot_gs: float, lifetime_s: float,
                            particle_mass: float | None = None) -> float:
    """Steady-state neutral particle inventory ``N = Ṁ·τ / m``.

    Args:
        mdot_gs: Mass source rate [g/s].
        lifetime_s: Loss (ionization) lifetime [s].
        particle_mass: Mass of one particle [g] (default Na atom mass).
    """
    particle_mass = CONST.m_Na if particle_mass is None else particle_mass
    return mdot_gs * lifetime_s / particle_mass
