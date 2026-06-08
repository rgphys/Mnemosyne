"""Energy-limited planetary wind → spectrum + lightcurve.

dishoom's energy-limited (XUV) escape sets the bulk mass-loss rate; a trace Na
fraction is advected in an isothermal **Parker wind** (Prometheus
``RadialWindExosphere``).  We pull both a phase-collapsed spectrum and a
band-integrated lightcurve from one cube.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import energetic_prometheus as ep
from energetic_prometheus import grids

planet = grids.find_planet("WASP-49b")

# Energy-limited escape: L_xuv ~ active-K-dwarf XUV; Na is a trace (xi=1e-3)
# species advected in the H/He Parker wind (wind_mu defaults to 1.3 amu).
model = ep.EnergyLimitedEscape(planet, L_xuv=1.4e28, T_wind=1e4, xi=1e-3,
                               wind_model="parker")
print(f"Na mass-loss rate: {model.mass_loss_rate():.3e} g/s")

# A multi-phase grid so we get a lightcurve (±2 h around mid-transit).
win = grids.orbphase_window_from_hours(planet, 2.0)
s_grid = grids.spatial_grid(planet, orbphase_window=win, orbphase_steps=11)

result = ep.EnergeticTransit(model, s_grid=s_grid,
                             use_phoenix_star=False).run()

# Phase-collapsed spectrum + depth.
depth_pct = result.transit_depth() * 100.0
print(f"Na D2 excess absorption (spectrum): {depth_pct:.3f} %")

# Lightcurve: line/continuum vs orbital phase (values < 1 = excess absorption).
lc = result.lightcurve()
phase_deg = np.degrees(result.orbphase)
print("\norbital phase [deg] -> lightcurve (line/continuum)")
for ph, val in zip(phase_deg, lc):
    print(f"  {ph:+6.2f}   {val:.5f}")
print(f"\ndeepest in-transit point: {lc.min():.5f}")
