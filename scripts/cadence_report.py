"""Print cadence and integration estimates for a fiducial XDS setup."""

from __future__ import annotations

from xds.config import InstrumentConfig, MoonConfig, PlanetConfig, SimulationConfig
from xds.observability import cadence_report
from xds.spectra import run_xds_model


def main() -> None:
    planet = PlanetConfig(v_eq=12.0, line_sigma=1.5)
    moon = MoonConfig(radius=0.28, semi_major_axis=3.0, period=5.0, phase=0.08)
    sim = SimulationConfig(
        transit_duration=5.0,
        n_times=401,
        n_longitudes=181,
        n_latitudes=121,
        use_illumination_shadow=True,
        use_observer_occultation=True,
    )
    instrument = InstrumentConfig(spectral_resolution=100_000.0, samples_per_contact=5)

    result = run_xds_model(planet, moon, sim)
    report = cadence_report(result, instrument)

    print("XDS cadence report")
    print(f"velocity_resolution_km_s = {report.velocity_resolution:.3f}")
    print(f"min_contact_duration_h = {report.min_contact_duration:.3f}")
    print(f"max_shadow_drift_rate_km_s_per_h = {report.max_shadow_drift_rate:.3f}")
    print(f"cadence_from_contacts_h = {report.cadence_from_contacts:.3f}")
    print(f"max_integration_from_smearing_h = {report.max_integration_from_smearing:.3f}")
    print(f"recommended_cadence_h = {report.recommended_cadence:.3f}")
    print(f"recommended_integration_time_h = {report.recommended_integration_time:.3f}")
    print("")
    for interval in report.intervals:
        print(
            f"{interval.kind:12s} start={interval.start:+.3f} h "
            f"end={interval.end:+.3f} h duration={interval.duration:.3f} h "
            f"peak_fraction={interval.peak_fraction:.4f}"
        )


if __name__ == "__main__":
    main()
