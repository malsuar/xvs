# Exomoon Velocity Shadow Toolkit

Research-code implementation of the **Exomoon Velocity Shadow** (XVS), a
Doppler-resolved distortion in planetary transmission spectra produced when an
orbiting satellite removes absorption from selected projected-velocity regions
of a rotating atmospheric terminator.

The package name is currently `xds` for historical compatibility with the
early Exomoon Doppler Shadow notation used during development.

## What The Model Does

The forward model builds time--velocity absorption maps for a thin effective
planetary terminator. It includes:

- a rigidly rotating planetary terminator;
- local Gaussian absorption profiles;
- a Keplerian satellite orbit in the planet-centred frame;
- an illumination blocking function, where the moon blocks stellar rays before
  they reach the atmospheric terminator;
- an observer-line blocking function, where the moon occults transmitted light
  on its way to the observer;
- differential maps relative to the moonless planet;
- simple diagnostics, cadence estimates, and approximate detectability tools.

The dimensionless maps use a local line contrast through
`PlanetConfig.line_depth`. Physical amplitudes in stellar-flux units are applied
separately by the detectability module through an effective atmospheric annulus
depth.

## Installation

From the repository root:

```bash
python -m pip install -e .
```

For development and notebooks:

```bash
python -m pip install -e ".[dev,notebooks]"
```

## Quick Start

Generate a Fig. 2-style time--velocity map:

```bash
python scripts/reproduce_fig2.py --preset paper-showcase
```

Generate the appendix geometry figure:

```bash
python scripts/appendix_terminator_geometry.py
```

Run the refined detectability search used for the manuscript figures:

```bash
python scripts/two_stage_lhs_search.py
python scripts/paper_snr_scatter_maps.py
```

Outputs are written to `outputs/`.

## Minimal Python Example

```python
from xds import MoonConfig, PlanetConfig, SimulationConfig, run_xds_model
from xds.diagnostics import map_diagnostics

planet = PlanetConfig(v_eq=12.0, line_sigma=1.5, line_depth=1.0)
moon = MoonConfig(radius=0.12, semi_major_axis=3.0, period=21.24, phase=0.08)
sim = SimulationConfig(transit_duration=16.0)

result = run_xds_model(planet, moon, sim)
diagnostics = map_diagnostics(result.time, result.velocity, result.delta_absorption)
print(result.delta_absorption.shape)
print(diagnostics.centroid[:5])
```

## Reproducing Manuscript Assets

The main reproducibility entry points are:

- `scripts/reproduce_fig2.py`: XVS time--velocity showcase.
- `scripts/appendix_terminator_geometry.py`: discretised terminator and blocking
  geometry.
- `scripts/two_stage_lhs_search.py`: heuristic pre-selection and forward-model
  refinement.
- `scripts/paper_snr_scatter_maps.py`: detectability scatter maps.
- `notebooks/02_interactive_xds_fig2.ipynb`: interactive Fig. 2 exploration.
- `notebooks/03_interactive_exovel_sweet_spot_scatter.ipynb`: interactive
  detectability figure.

Large generated `.npz`, `.png`, and `.pdf` files in `outputs/` are intentionally
not required for installation. Regenerate them from the scripts above.

## Scientific Scope

### Molecular Injection Pilot

A separate CO line-by-line pilot uses traceable ExoMol data, radial atmospheric
opacity, instrumental resolution and exposure integration. It compares a
known-trajectory lunar signal with independent line-depth changes per exposure:

```bash
python -m pip install -e '.[molecular]'
python scripts/fetch_co_template_data.py
PYTHONPATH=src python scripts/molecular_injection_recovery.py
PYTHONPATH=src python scripts/audit_co_template.py
```

See [the experiment report](docs/co_template_injection_recovery.md) for assumptions,
convergence and injection controls. The native pilot signal is not detectable
under the adopted noise assumptions; amplified recovery controls are numerical
tests, not physical detection predictions. ExoMol data have a separate license.

The follow-up [atmospheric-degeneracy experiment](docs/co_atmospheric_degeneracy.md)
compares uniform and two-limb temperature, CO abundance and wind changes:

```bash
PYTHONPATH=src python scripts/atmospheric_nuisance_recovery.py
PYTHONPATH=src python scripts/audit_atmospheric_nuisance.py
```

It tests known-orbit recovery after atmospheric fitting, not a blind search or
a proof of uniqueness against arbitrary atmospheric structure.

The [bounded nonlinear follow-up](docs/co_nonlinear_coherent_recovery.md) fits
the actual CO family with smooth, bounded atmospheric histories and a coupled
transparent-to-opaque lunar mask:

```bash
PYTHONPATH=src python scripts/nonlinear_co_recovery.py --n-mc 16
```

Its reduced-noise control keeps the physical moon unamplified. Reported
`sqrt(delta chi2)` values are model separations, not calibrated discovery
significances. Cached spectral grids can be regenerated and are not tracked.

The next [phase/period-blind pilot](docs/co_blind_search_report.md) calibrates
a finite score search against independent atmospheric nulls. Its scope and
predeclared choices are in the [direction review and protocol](docs/research_direction_blind_protocol.md).

```bash
PYTHONPATH=src python scripts/build_blind_co_library.py
PYTHONPATH=src python scripts/validate_blind_co.py
PYTHONPATH=src python scripts/blind_co_recovery.py --tag pilot
PYTHONPATH=src python scripts/audit_blind_grid.py
PYTHONPATH=src python scripts/summarize_blind_comparisons.py --tag pilot
PYTHONPATH=src python scripts/plot_blind_co.py --tag pilot
```

The nonlinear CO emulator from the preceding step is required. The main pilot
performs 1408 searches on 704 simulated data sets and can take tens of minutes.
An identical command resumes completed trials; changed code/configuration
requires a new tag. Empirical probabilities apply to the calibrated score,
not to the subsequent nonlinear candidate-fit likelihood ratio.

### Limitations

This is a minimal forward model intended for reproducibility and exploratory
extensions. It is not yet a full retrieval or population-inference framework.
Current assumptions include a thin effective terminator, solid-body rotation,
Gaussian local line profiles, parallel stellar rays, geometric blocking masks,
and approximate detectability scaling.

## Development

Run the test suite with:

```bash
python -m pytest
```

## Citation

If you use this code, please cite the accompanying manuscript and the software
record described in `CITATION.cff`.
