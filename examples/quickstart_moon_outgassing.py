"""Quickstart: moon thermal outgassing → Na D transmission spectrum.

dishoom computes the Na source rate from a moon's thermal (Hertz-Knudsen)
sublimation; Prometheus turns the resulting exosphere into a Na-D transit depth.
Run from anywhere::

    python "examples/quickstart_moon_outgassing.py"
"""

import os
import sys

# Make the package importable without installation.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import energetic_prometheus as ep
from energetic_prometheus import grids

# 1. System: WASP-49b with an Io-like moon at 1.7 R_p.
planet = grids.find_planet("WASP-49b")
moon = grids.make_moon(planet, a_over_Rp=1.7)

# 2. dishoom energetics: thermal sublimation at a 1500 K sub-stellar point,
#    with the volcanic Na mass fraction (dishoom default).
model = ep.ThermalSublimation(moon, T_surface=1500.0)
s = model.summary()
print(f"Na source rate : {s['mass_loss_rate_gs']:.3e} g/s")
print(f"ionization life: {s['tau_ion_s']:.1f} s")
print(f"steady Na count: {s['N_particles']:.3e} atoms")

# 3. Prometheus transit (flat star → fast; drop use_phoenix_star for a real
#    PHOENIX stellar spectrum).
result = ep.EnergeticTransit(model, use_phoenix_star=False).run()

# 4. Outputs: spectrum array + a single excess-absorption number.
wav_ang = result.wavelength_ang
spec = result.spectrum_normalized()
depth_pct = result.transit_depth() * 100.0
print(f"\nNa D2 excess absorption: {depth_pct:.3f} %")
print(f"spectrum: {spec.size} points over "
      f"{wav_ang.min():.1f}-{wav_ang.max():.1f} A")
