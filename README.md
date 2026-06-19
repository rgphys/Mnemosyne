# Mnemosyne

In-memory coupling of **dishoom** (atmospheric escape / energetics) and
**Prometheus** (radiative transfer). Compute a moon or planet mass-loss rate
with dishoom and feed it *directly* into a Prometheus transit — no JSON setup
files, no output-file handoff. You get the transmission spectrum and lightcurve
back as NumPy arrays.

```
dishoom energetics  ►  Prometheus scenario  ►  chord-summed transit
 (Ṁ, N, exobase)         (Exosphere / Wind)        (R(phase, λ))
```

## Why

Previously, dishoom and Prometheus exchanged data through `setupFiles/*.txt` →
`output/*.txt`. Mnemosyne replaces that with a typed Python API so
escape physics and radiative transfer compose in a single process (sweeps,
retrievals, notebooks).

* **RT engine:** the actively-developed `Prometheus/` (not dishoom's bundled
  2022 fork). Its `RadialWindExosphere` already implements the exact isothermal
  Parker wind, so dishoom only needs to supply the *energetics* (Ṁ, N).
* **Energetics:** the real dishoom functions in `moon_functions_5.py` are
  imported through a thin adapter — dishoom stays canonical.

## Install / import

No packaging step. The package locates `../Prometheus` and `../dishoom`
automatically. Use the project venv (`source env/bin/activate`) and put the
folder on `sys.path`:

```python
import sys; sys.path.insert(0, "Mnemosyne")
import mnemosyne as mn
```

## Quick start

```python
import mnemosyne as mn
from mnemosyne import grids

planet = grids.find_planet("WASP-49b")
moon   = grids.make_moon(planet, a_over_Rp=1.7)

# dishoom: thermal sublimation source rate → Prometheus: Na-D transit
model  = mn.ThermalSublimation(moon, T_surface=1500.0)
result = mn.EnergeticTransit(model, use_phoenix_star=False).run()

print(result.transit_depth() * 100, "% Na D2")   # excess absorption
wav, spec = result.wavelength_ang, result.spectrum_normalized()
```

## Two API layers

### Headline — `EnergeticTransit`

Runs the whole pipeline and returns a `TransitResult` with:

* `.spectrum()` / `.spectrum_normalized()` — phase-collapsed `R(λ)`
* `.transit_depth(line_window_ang, continuum_exclude_ang, mode)` — excess
  absorption (peak or band-mean)
* `.lightcurve(line_window_ang, continuum_exclude_ang)` — line/continuum vs
  orbital phase (needs `orbphase_steps > 1`)
* `.wavelength_ang` / `.wavelength_um`, `.R_2D`, `.orbphase`

### Composable — escape models, scenario builders, grids

```python
from mnemosyne import grids, radial_wind_scenario

# build a Prometheus scenario from a raw mass-loss rate yourself:
w_grid = grids.na_d_grid()
scen   = radial_wind_scenario(Mdot=2e7, planet=planet, species="NaI",
                              sigma_v=2e6, w_grid=w_grid)
# ... drop `scen` into gasProperties.Atmosphere([...]) as usual.
```

> **These builders now live in Prometheus.** The grid/body/orbital helpers and
> the scenario builders were upstreamed into Prometheus because they carry no
> dishoom dependency — `gasProperties.na_d_grid` / `line_grid` /
> `moon_exosphere_scenario` / `powerlaw_exosphere_scenario` /
> `radial_wind_scenario` / `run_transit` / `TransitResult`,
> `geometryHandler.spatial_grid` / `orbphase_window_from_hours`, and
> `celestialBodies.find_planet` / `make_moon` / `mean_motion_ratio` /
> `optimal_midtransit_phase` / `max_peak_shift_minutes`. The `mnemosyne.grids`,
> `mnemosyne.scenarios`, and `mnemosyne.TransitResult` names are kept as
> backward-compatible re-exports (the only Mnemosyne-side difference is the
> `w_grid` keyword, which Prometheus spells `wavelengthGrid`). New code can
> import these straight from Prometheus.

## Escape models

| Model | dishoom physics | Prometheus scenario |
|-------|-----------------|---------------------|
| `ThermalSublimation`  | Hertz-Knudsen vapor flux (`Mdot0_thermal_notides`) | `MoonExosphere` (N) |
| `TidalHeating`        | tidal power → source (`Mdot0_tidal` / `Mdot0_tot`) | `MoonExosphere` (N) |
| `SurfaceJeansEscape`  | Volkov+ 2011 Jeans fraction (`dMdt_SJeansMdot0_tot_Na`); `.exobase_radius()` | `MoonExosphere` (N) |
| `EnergyLimitedEscape` | energy-limited XUV escape (`dMdt_ELescape[_xuv]`) | `RadialWindExosphere` (Ṁ, Parker/beta) |

Moon models normalize by particle number `N = Ṁ·τ_ion / m`; the wind model
normalizes by mass continuity `n(r) = Ṁ / (4π r² v(r) μ)`.

## Layout

```
mnemosyne/
  _bootstrap.py       path wiring + dishoom stdout silencing
  dishoom_adapter.py  typed wrappers over dishoom escape functions + constants
  constants.py        mineral vapor-pressure DB + escape-physics constants
  escape.py           EscapeModel classes (the dishoom-coupled physics layer)
  scenarios.py        thin adapter → Prometheus scenario builders (keeps w_grid)
  grids.py            re-exports of Prometheus grid/body/orbital helpers
  transit.py          EnergeticTransit (dishoom coupling) + TransitResult re-export
examples/             runnable quickstart + wind/lightcurve scripts
```

The pure radiative-transfer / orbital-mechanics pieces (grids, scenario
builders, `TransitResult`, the `run_transit` orchestration, and the moon
orbital-mechanics relations) live in **Prometheus**; what remains here is the
dishoom coupling — `dishoom_adapter`, the escape/source `constants`, the
`escape` models, and the `EnergeticTransit` driver.

## Notes & caveats

* **Units are CGS** throughout (lengths cm, Ṁ in g/s, σ_v in cm/s).
* **`use_phoenix_star=False`** uses a flat stellar continuum — much faster and
  fine for *relative* depths. Pass `True` for a real PHOENIX spectrum.
* **Trace species in winds:** for energy-limited escape, Na is a *trace*
  species — set `xi` (e.g. `1e-3`) so the Na mass-loss rate is a fraction of the
  bulk; `xi=1.0` treats all escaping mass as Na and saturates the line.
* **dishoom compatibility shim:** dishoom's `Mdot0_tot` calls a function
  (`Mdot0_thermal`) that was renamed to `Mdot0_thermal_notides` and never
  updated, so `Mdot0_tot` / `dMdt_SJeansMdot0_tot_Na` raise `NameError`
  upstream. `dishoom_adapter` restores the missing symbol **in memory only**
  (dishoom on disk is untouched) so the canonical routines run. See the comment
  in `dishoom_adapter.py`.
* **Moon mass:** Prometheus' `Moon` carries no mass; tidal/Jeans models default
  `moon_mass` to Io's mass (dishoom's reference body) — override as needed.
