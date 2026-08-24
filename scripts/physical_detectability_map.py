"""Generate a first physical XDS detectability map.

The script keeps the dimensionless XDS morphology fixed and varies the
atmospheric scale height with an equilibrium-temperature proxy. This is not a
full exposure-time calculator, but it turns the conceptual efficiency plot into
a falsifiable SNR scaling with explicit assumptions.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from xds.config import (
    AtmosphereConfig,
    InstrumentConfig,
    MoonConfig,
    ObservationConfig,
    PlanetConfig,
    SimulationConfig,
)
from xds.detectability import detectability_report
from xds.efficiency import hill_radius_over_planet_radius
from xds.spectra import run_xds_model


def equilibrium_temperature(period_days: np.ndarray, stellar_mass_solar: float = 1.0) -> np.ndarray:
    """Zero-albedo, full-redistribution proxy for a Sun-like star."""

    period_years = np.asarray(period_days) / 365.25
    semi_major_axis_au = (period_years**2 * stellar_mass_solar) ** (1.0 / 3.0)
    return 278.0 / np.sqrt(semi_major_axis_au)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "outputs" / "physical_detectability_map.png"
    data_output = root / "outputs" / "physical_detectability_map.npz"
    output.parent.mkdir(parents=True, exist_ok=True)

    planet = PlanetConfig(v_eq=12.0, line_sigma=1.5)
    moon = MoonConfig(radius=0.28, semi_major_axis=3.0, period=5.0, phase=0.08)
    sim = SimulationConfig(
        transit_duration=5.0,
        n_times=241,
        n_longitudes=121,
        n_latitudes=81,
        n_velocities=301,
        use_illumination_shadow=True,
        use_observer_occultation=True,
    )
    result = run_xds_model(planet, moon, sim)

    instrument = InstrumentConfig(spectral_resolution=100_000.0, samples_per_contact=5)
    observation = ObservationConfig(
        continuum_snr_per_resolution_element=200.0,
        exposure_time=1.0 / 60.0,
        n_transits=1,
        n_lines=30,
        systematic_floor=20.0e-6,
    )

    period = np.logspace(0, 4, 160)
    temperature = equilibrium_temperature(period)
    hill = hill_radius_over_planet_radius(period)
    snr = np.empty_like(period)
    line_depth = np.empty_like(period)
    peak_depth = np.empty_like(period)
    scale_height = np.empty_like(period)

    for idx, (period_days, temp) in enumerate(zip(period, temperature)):
        atmosphere = AtmosphereConfig(
            temperature=float(temp),
            mean_molecular_weight=2.3,
            effective_scale_heights=5.0,
            line_contrast=1.0,
            planet_mass_jupiter=1.0,
            planet_radius_jupiter=1.0,
            stellar_radius_solar=1.0,
        )
        report = detectability_report(result, atmosphere, observation, instrument)
        snr[idx] = report.matched_filter_snr
        line_depth[idx] = report.atmospheric_line_depth
        peak_depth[idx] = report.peak_xds_depth
        scale_height[idx] = report.scale_height_km

    np.savez(
        data_output,
        period_days=period,
        temperature=temperature,
        hill_radius_rp=hill,
        snr=snr,
        atmospheric_line_depth=line_depth,
        peak_xds_depth=peak_depth,
        scale_height_km=scale_height,
    )

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.6), sharex=True, constrained_layout=True)
    axes[0].plot(period, snr, color="black", linewidth=2.0)
    axes[0].axhline(5.0, color="tab:red", linestyle="--", linewidth=1.2, label="SNR = 5")
    axes[0].set_ylabel("Matched-filter SNR")
    axes[0].legend()

    axes[1].plot(period, 1.0e6 * peak_depth, color="tab:blue", linewidth=2.0)
    axes[1].set_ylabel("Peak XDS depth [ppm]")

    axes[2].plot(period, temperature, color="tab:orange", linewidth=2.0, label="T_eq")
    axes[2].set_ylabel("T_eq [K]")
    axes[2].set_xlabel("Planet orbital period [days]")
    axes[2].set_xscale("log")
    ax2 = axes[2].twinx()
    ax2.plot(period, hill, color="tab:green", linewidth=1.5, alpha=0.8, label="Hill radius")
    ax2.set_ylabel("Hill radius [Rp]")

    for ax in axes:
        ax.axvspan(1.0, 10.0, color="tab:red", alpha=0.08)
        ax.axvspan(10.0, 300.0, color="tab:green", alpha=0.08)
        ax.axvspan(300.0, 10_000.0, color="tab:blue", alpha=0.08)

    fig.savefig(output, dpi=180)
    print(output)
    print(data_output)


if __name__ == "__main__":
    main()
