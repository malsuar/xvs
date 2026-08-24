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
