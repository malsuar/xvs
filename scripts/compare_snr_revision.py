"""Controlled estimator comparison; not a population detectability forecast."""
import csv
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from xds import PlanetConfig, MoonConfig, SimulationConfig, run_xds_model
from xds.config import AtmosphereConfig, InstrumentConfig, ObservationConfig
from xds.detectability import detectability_report, atmospheric_line_depth
from xds.observing import observe_residual


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / "outputs"
    out.mkdir(exist_ok=True)
    atmosphere = AtmosphereConfig(planet_mass_jupiter=.7, planet_radius_jupiter=1.1,
                                 stellar_radius_solar=1., temperature=500.)
    planet = PlanetConfig(v_eq=12., line_sigma=1.5)
    moon = MoonConfig(radius=.0892, semi_major_axis=3., period=21.238292071543107, phase=.08)
    obs = ObservationConfig(continuum_snr_per_resolution_element=300.,
                            exposure_time=5/60, n_lines=200, systematic_floor=10e-6)
    inst = InstrumentConfig()
    base = SimulationConfig(transit_duration=16., geometry="projected_annulus",
                            spin_convention="observer", n_azimuths=240, n_radii=4,
                            n_times=1921, n_velocities=801, velocity_min=-40, velocity_max=40)
    cases = [("paper_script_radius_0.0892", planet, moon),
             ("paper_notebook_radius_0.12", planet, replace(moon, radius=.12)),
             ("slow_spin_3kms", replace(planet, v_eq=3.), moon),
             ("larger_moon_0.2", planet, replace(moon, radius=.2))]
    records = []
    first = None
    for name, p, m in cases:
        result = run_xds_model(p, m, base)
        old = detectability_report(result, atmosphere, obs, inst, method="legacy")
        new = detectability_report(result, atmosphere, obs, inst)
        records.append(dict(case=name, legacy_snr=old.matched_filter_snr,
                            observed_snr=new.matched_filter_snr,
                            ratio=new.matched_filter_snr/old.matched_filter_snr,
                            observed_peak_ppm=1e6*new.peak_xds_depth))
        if first is None:
            first = result
    with (root/"docs/snr_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    # Subsample one resolved map: the physical observing sequence stays fixed.
    convergence = []
    for stride in [4, 2, 1]:
        r = replace(first, time=first.time[::stride], velocity=first.velocity[::stride],
                    delta_absorption=first.delta_absorption[::stride, ::stride],
                    illumination_fraction=first.illumination_fraction[::stride],
                    observer_fraction=first.observer_fraction[::stride],
                    combined_fraction=first.combined_fraction[::stride])
        # Legacy drift diagnostic needs only these result arrays.
        old = detectability_report(r, atmosphere, obs, inst, method="legacy")
        new = detectability_report(r, atmosphere, obs, inst)
        convergence.append((len(r.time)*len(r.velocity), old.matched_filter_snr, new.matched_filter_snr))
    measured = observe_residual(first, atmospheric_line_depth(atmosphere), obs, inst)
    plt.rcParams.update({"font.size": 11})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    axes[0].plot(first.time, 1e6*atmospheric_line_depth(atmosphere)*np.max(abs(first.delta_absorption), axis=1), label="Intrinsic")
    axes[0].plot(measured.time, 1e6*np.max(abs(measured.signal), axis=1), label="LSF + 5 min exposures")
    axes[0].set(xlabel="Time [h]", ylabel="Peak deficit [ppm]")
    axes[0].legend()
    row = np.argmax(np.max(abs(measured.signal), axis=1))
    i = np.argmin(abs(first.time-measured.time[row]))
    axes[1].plot(first.velocity, -1e6*atmospheric_line_depth(atmosphere)*first.delta_absorption[i], label="Instantaneous")
    axes[1].step(measured.velocity, -1e6*measured.signal[row], where="mid", label="Observed bins")
    axes[1].set(xlabel="Velocity toward observer [km/s]", ylabel="Deficit [ppm]")
    axes[1].legend()
    a = np.array(convergence)
    axes[2].plot(a[:, 0], a[:, 1]/a[-1, 1], "o--", label="Legacy / finest")
    axes[2].plot(a[:, 0], a[:, 2]/a[-1, 2], "s-", label="Corrected / finest")
    axes[2].set(xlabel="Internal time x velocity samples", ylabel="Relative S/N")
    axes[2].legend()
    for ext in ["png", "pdf"]:
        fig.savefig(out/f"snr_revision_diagnostics.{ext}", dpi=180)
    print(records)
    print("CONVERGENCE", convergence)


if __name__ == "__main__":
    main()
