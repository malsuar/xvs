"""Print a physical detectability budget for the fiducial XDS model."""

from __future__ import annotations

from xds.config import (
    AtmosphereConfig,
    InstrumentConfig,
    MoonConfig,
    ObservationConfig,
    PlanetConfig,
    SimulationConfig,
)
from xds.detectability import detectability_report, detection_requirements
from xds.spectra import run_xds_model


def main() -> None:
    result = run_xds_model(
        PlanetConfig(v_eq=12.0, line_sigma=1.5),
        MoonConfig(radius=0.28, semi_major_axis=3.0, period=5.0, phase=0.08),
        SimulationConfig(
            transit_duration=5.0,
            n_times=241,
            n_longitudes=121,
            n_latitudes=81,
            n_velocities=301,
            use_illumination_shadow=True,
            use_observer_occultation=True,
        ),
    )
    atmosphere = AtmosphereConfig(
        temperature=500.0,
        mean_molecular_weight=2.3,
        effective_scale_heights=5.0,
        line_contrast=1.0,
        planet_mass_jupiter=1.0,
        planet_radius_jupiter=1.0,
        stellar_radius_solar=1.0,
    )
    observation = ObservationConfig(
        continuum_snr_per_resolution_element=200.0,
        exposure_time=1.0 / 60.0,
        n_transits=1,
        n_lines=30,
        systematic_floor=20.0e-6,
    )
    instrument = InstrumentConfig(spectral_resolution=100_000.0)

    report = detectability_report(result, atmosphere, observation, instrument)
    requirements = detection_requirements(report, observation, target_snr=5.0)

    print("Physical XDS detectability budget")
    print(f"scale_height_km = {report.scale_height_km:.2f}")
    print(f"atmospheric_line_depth_ppm = {1e6 * report.atmospheric_line_depth:.2f}")
    print(f"peak_xds_depth_ppm = {1e6 * report.peak_xds_depth:.2f}")
    print(f"peak_pixel_snr = {report.peak_pixel_snr:.3f}")
    print(f"matched_filter_snr = {report.matched_filter_snr:.3f}")
    print(f"n_effective_pixels = {report.n_effective_pixels}")
    print(f"smearing_penalty = {report.smearing_penalty:.3f}")
    print(f"cadence_limited_exposure_min = {60.0 * report.cadence_limited_exposure:.2f}")
    print(f"required_lines_one_transit_for_snr5 = {requirements.required_lines_one_transit:.1f}")
    print(f"required_transits_with_current_lines_for_snr5 = {requirements.required_transits:.1f}")
    if report.notes:
        print("notes = " + ", ".join(report.notes))


if __name__ == "__main__":
    main()
