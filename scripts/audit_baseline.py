"""Reproduce geometry and sampling checks without modifying the baseline model."""

import json
from dataclasses import replace

import numpy as np

from xds import MoonConfig, PlanetConfig, SimulationConfig, run_xds_model
from xds.config import AtmosphereConfig, InstrumentConfig, ObservationConfig
from xds.detectability import detectability_report
from xds.rotation import spin_axis
from xds.shadow import compute_shadow_masks
from xds.terminator import build_terminator


def main():
    grid = build_terminator(81, 55)
    projected_radius = np.linalg.norm(grid.positions[:, 1:], axis=1)
    summary = {
        "geometry": {
            "weight_sum": float(grid.weights.sum()),
            "weight_inside_projected_radius_0.9": float(grid.weights[projected_radius < 0.9].sum()),
            "negative_y_elements": int(np.sum(grid.positions[:, 1] < 0)),
            "spin_angle_to_observer_for_input_i0_deg": float(np.rad2deg(np.arccos(np.dot(spin_axis(0, 0), [-1, 0, 0])))),
            "spin_angle_to_observer_for_input_i90_deg": float(np.rad2deg(np.arccos(np.dot(spin_axis(90, 0), [-1, 0, 0])))),
        },
    }
    point = np.array([[0., 1., 0.]])
    upstream = compute_shadow_masks(point, np.array([2., 1., 0.]), .1, (-1, 0, 0), (-1, 0, 0))
    downstream = compute_shadow_masks(point, np.array([-2., 1., 0.]), .1, (-1, 0, 0), (-1, 0, 0))
    assert upstream.illumination[0] and not upstream.observer[0]
    assert downstream.observer[0] and not downstream.illumination[0]
    summary["ray_ordering_control_passed"] = True
    sim = SimulationConfig(n_times=61, n_velocities=121, n_longitudes=41, n_latitudes=31)
    rows = []
    for nt, nv in [(61, 121), (121, 121), (61, 241), (121, 241)]:
        result = run_xds_model(PlanetConfig(), MoonConfig(), replace(sim, n_times=nt, n_velocities=nv))
        report = detectability_report(result, AtmosphereConfig(), ObservationConfig(), InstrumentConfig())
        rows.append({"n_times": nt, "n_velocities": nv, "snr": report.matched_filter_snr})
    for row in rows:
        row["snr_ratio_to_base"] = row["snr"] / rows[0]["snr"]
    summary["sampling_at_fixed_physical_settings"] = rows
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
