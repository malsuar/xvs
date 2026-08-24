"""Search for physically favourable XDS configurations.

This is an exploratory optimisation, not a population-synthesis model. It scans
planet mass/radius, stellar radius, planet orbital period, spin period, and moon
separation. Each candidate is filtered by Roche and Hill-stability limits, then
evaluated with the forward XDS model and the physical detectability estimator.
"""

from __future__ import annotations

import argparse
from pathlib import Path

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
from xds.dynamics import (
    breakup_rotation_period_hours,
    equatorial_velocity_km_s,
    fluid_roche_limit_rp,
    moon_period_hours,
    planet_density_g_cm3,
    stable_outer_moon_limit_rp,
    transit_duration_hours,
)
from xds.observability import cadence_report
from xds.spectra import run_xds_model


def equilibrium_temperature(
    planet_period_days: np.ndarray | float,
    stellar_radius_solar: np.ndarray | float,
    stellar_temperature: float = 5200.0,
    star_mass_solar: float = 0.8,
) -> np.ndarray:
    """Zero-albedo full-redistribution estimate."""

    period_years = np.asarray(planet_period_days, dtype=float) / 365.25
    a_au = (period_years**2 * star_mass_solar) ** (1.0 / 3.0)
    rstar_au = np.asarray(stellar_radius_solar, dtype=float) * 0.00465047
    return stellar_temperature * np.sqrt(rstar_au / (2.0 * a_au))


def sample_candidates(n: int = 220, seed: int = 7) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    return {
        "planet_mass_jupiter": 10 ** rng.uniform(np.log10(0.08), np.log10(2.0), n),
        "planet_radius_jupiter": rng.uniform(0.75, 1.75, n),
        "stellar_radius_solar": rng.uniform(0.45, 1.05, n),
        "planet_period_days": 10 ** rng.uniform(np.log10(3.0), np.log10(300.0), n),
        "rotation_period_hours": 10 ** rng.uniform(np.log10(3.0), np.log10(24.0), n),
        "moon_separation_rp": rng.uniform(1.8, 12.0, n),
        "moon_radius_rp": rng.uniform(0.12, 0.35, n),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-candidates", type=int, default=220)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument(
        "--tag",
        default=None,
        help="Output suffix. Defaults to n{n_candidates}_seed{seed}.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    output_dir = root / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    tag = args.tag or f"n{args.n_candidates}_seed{args.seed}"

    instrument = InstrumentConfig(spectral_resolution=100_000.0, samples_per_contact=5)
    observation = ObservationConfig(
        continuum_snr_per_resolution_element=300.0,
        exposure_time=1.0 / 60.0,
        n_transits=1,
        n_lines=200,
        systematic_floor=10.0e-6,
    )
    sim_base = dict(
        n_times=121,
        n_longitudes=81,
        n_latitudes=55,
        n_velocities=241,
        velocity_min=-45.0,
        velocity_max=45.0,
        use_illumination_shadow=True,
        use_observer_occultation=True,
    )

    samples = sample_candidates(n=args.n_candidates, seed=args.seed)
    n = len(samples["planet_mass_jupiter"])
    records: list[dict[str, float]] = []

    for idx in range(n):
        mp = float(samples["planet_mass_jupiter"][idx])
        rp = float(samples["planet_radius_jupiter"][idx])
        rstar = float(samples["stellar_radius_solar"][idx])
        pplanet = float(samples["planet_period_days"][idx])
        prot = float(samples["rotation_period_hours"][idx])
        am = float(samples["moon_separation_rp"][idx])
        rm = float(samples["moon_radius_rp"][idx])

        density = float(planet_density_g_cm3(mp, rp))
        roche = float(fluid_roche_limit_rp(density, moon_density_g_cm3=3.0))
        outer = float(stable_outer_moon_limit_rp(pplanet, mp, star_mass_solar=0.8, planet_radius_jupiter=rp))
        breakup = float(breakup_rotation_period_hours(mp, rp))
        if am <= max(roche, 1.0 + rm) or am >= outer:
            continue
        if prot <= 1.15 * breakup:
            continue

        transit_duration = float(transit_duration_hours(pplanet, rstar, star_mass_solar=0.8))
        if not np.isfinite(transit_duration) or transit_duration < 1.5:
            continue
        transit_duration = min(transit_duration, 16.0)

        teq = float(equilibrium_temperature(pplanet, rstar))
        if teq < 250.0:
            continue
        veq = float(equatorial_velocity_km_s(rp, prot))
        pm_hours = float(moon_period_hours(am, mp, rp))

        planet = PlanetConfig(v_eq=veq, line_sigma=1.5, line_depth=1.0)
        moon = MoonConfig(
            radius=rm,
            semi_major_axis=am,
            period=pm_hours,
            phase=0.08,
            inclination=90.0,
        )
        sim = SimulationConfig(transit_duration=transit_duration, **sim_base)
        result = run_xds_model(planet, moon, sim)
        atmosphere = AtmosphereConfig(
            temperature=teq,
            mean_molecular_weight=2.3,
            effective_scale_heights=5.0,
            line_contrast=1.0,
            planet_mass_jupiter=mp,
            planet_radius_jupiter=rp,
            stellar_radius_solar=rstar,
        )
        report = detectability_report(result, atmosphere, observation, instrument)
        cadence = cadence_report(result, instrument)

        records.append(
            {
                "snr": report.matched_filter_snr,
                "peak_xds_ppm": 1.0e6 * report.peak_xds_depth,
                "line_depth_ppm": 1.0e6 * report.atmospheric_line_depth,
                "scale_height_km": report.scale_height_km,
                "v_eq_km_s": veq,
                "planet_mass_jupiter": mp,
                "planet_radius_jupiter": rp,
                "stellar_radius_solar": rstar,
                "planet_period_days": pplanet,
                "rotation_period_hours": prot,
                "breakup_period_hours": breakup,
                "moon_separation_rp": am,
                "moon_period_hours": pm_hours,
                "moon_radius_rp": rm,
                "roche_limit_rp": roche,
                "outer_stability_limit_rp": outer,
                "transit_duration_hours": transit_duration,
                "temperature_k": teq,
                "min_contact_duration_hours": cadence.min_contact_duration,
                "integration_limit_minutes": 60.0 * report.cadence_limited_exposure,
                "smearing_penalty": report.smearing_penalty,
            }
        )

    if not records:
        raise RuntimeError("No valid candidate survived the dynamical filters.")

    names = list(records[0])
    arrays = {name: np.array([record[name] for record in records], dtype=float) for name in names}
    np.savez(output_dir / f"xds_sweet_spot_search_{tag}.npz", **arrays)

    order = np.argsort(arrays["snr"])[::-1]
    print("Top XDS candidates")
    for rank, row in enumerate(order[:10], start=1):
        print(
            f"{rank:02d} SNR={arrays['snr'][row]:.2f} "
            f"peak={arrays['peak_xds_ppm'][row]:.1f} ppm "
            f"Mp={arrays['planet_mass_jupiter'][row]:.2f} Mj "
            f"Rp={arrays['planet_radius_jupiter'][row]:.2f} Rj "
            f"Teq={arrays['temperature_k'][row]:.0f} K "
            f"Prot={arrays['rotation_period_hours'][row]:.1f} h "
            f"veq={arrays['v_eq_km_s'][row]:.1f} km/s "
            f"am={arrays['moon_separation_rp'][row]:.1f} Rp "
            f"Pm={arrays['moon_period_hours'][row]:.1f} h"
        )

    make_plots(arrays, output_dir, tag)
    make_moon_separation_scan(arrays, order[0], output_dir, observation, instrument, sim_base, tag)


def make_plots(arrays: dict[str, np.ndarray], output_dir: Path, tag: str) -> None:
    import matplotlib.pyplot as plt

    snr = arrays["snr"]

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.4), constrained_layout=True)
    ax = axes[0, 0]
    sc = ax.scatter(
        arrays["planet_mass_jupiter"],
        arrays["planet_radius_jupiter"],
        c=snr,
        s=45,
        cmap="viridis",
        edgecolor="none",
    )
    ax.set_xscale("log")
    ax.set_xlabel("Planet mass [Mj]")
    ax.set_ylabel("Planet radius [Rj]")
    ax.set_title("Mass-radius sweet spots")
    fig.colorbar(sc, ax=ax, label="SNR")

    ax = axes[0, 1]
    sc = ax.scatter(
        arrays["temperature_k"],
        arrays["v_eq_km_s"],
        c=snr,
        s=45,
        cmap="plasma",
        edgecolor="none",
    )
    ax.set_xlabel("Equilibrium temperature [K]")
    ax.set_ylabel("Equatorial velocity [km/s]")
    ax.set_title("Thermal column vs Doppler leverage")
    fig.colorbar(sc, ax=ax, label="SNR")

    ax = axes[1, 0]
    sc = ax.scatter(
        arrays["moon_separation_rp"],
        arrays["moon_period_hours"] / arrays["transit_duration_hours"],
        c=snr,
        s=45,
        cmap="magma",
        edgecolor="none",
    )
    ax.set_xlabel("Moon separation [Rp]")
    ax.set_ylabel("Moon period / transit duration")
    ax.set_title("Satellite timing optimisation")
    fig.colorbar(sc, ax=ax, label="SNR")

    ax = axes[1, 1]
    sc = ax.scatter(
        arrays["min_contact_duration_hours"],
        arrays["integration_limit_minutes"],
        c=snr,
        s=45,
        cmap="cividis",
        edgecolor="none",
    )
    ax.set_xlabel("Shortest contact duration [h]")
    ax.set_ylabel("Velocity-smearing exposure limit [min]")
    ax.set_title("Cadence feasibility")
    fig.colorbar(sc, ax=ax, label="SNR")

    fig.savefig(output_dir / f"xds_sweet_spot_diagnostics_{tag}.png", dpi=180)

    best = np.argsort(snr)[::-1][:40]
    labels = [
        "planet_mass_jupiter",
        "planet_radius_jupiter",
        "temperature_k",
        "v_eq_km_s",
        "moon_separation_rp",
        "moon_period_hours",
        "stellar_radius_solar",
        "snr",
    ]
    values = np.vstack([arrays[label][best] for label in labels])
    normed = (values - np.nanmin(values, axis=1)[:, None]) / (
        np.nanmax(values, axis=1)[:, None] - np.nanmin(values, axis=1)[:, None] + 1.0e-12
    )

    fig2, ax2 = plt.subplots(figsize=(10.0, 5.0), constrained_layout=True)
    image = ax2.imshow(normed, aspect="auto", cmap="viridis", origin="lower")
    ax2.set_yticks(np.arange(len(labels)), labels=[label.replace("_", " ") for label in labels])
    ax2.set_xlabel("Top candidates ranked by SNR")
    ax2.set_title("Normalised parameters among the top 40 candidates")
    fig2.colorbar(image, ax=ax2, label="normalised value")
    fig2.savefig(output_dir / f"xds_top_candidate_parameter_matrix_{tag}.png", dpi=180)


def make_moon_separation_scan(
    arrays: dict[str, np.ndarray],
    best_idx: int,
    output_dir: Path,
    observation: ObservationConfig,
    instrument: InstrumentConfig,
    sim_base: dict,
    tag: str,
) -> None:
    import matplotlib.pyplot as plt

    mp = float(arrays["planet_mass_jupiter"][best_idx])
    rp = float(arrays["planet_radius_jupiter"][best_idx])
    rstar = float(arrays["stellar_radius_solar"][best_idx])
    pplanet = float(arrays["planet_period_days"][best_idx])
    prot = float(arrays["rotation_period_hours"][best_idx])
    rm = float(arrays["moon_radius_rp"][best_idx])
    teq = float(arrays["temperature_k"][best_idx])
    transit_duration = float(arrays["transit_duration_hours"][best_idx])

    roche = float(arrays["roche_limit_rp"][best_idx])
    outer = float(arrays["outer_stability_limit_rp"][best_idx])
    lower = max(1.02 * roche, 1.05 * (1.0 + rm))
    upper = min(12.0, 0.95 * outer)
    if upper <= lower:
        return

    separations = np.linspace(lower, upper, 80)
    snr = np.empty_like(separations)
    contact = np.empty_like(separations)
    exposure_limit = np.empty_like(separations)
    peak_ppm = np.empty_like(separations)
    moon_period = moon_period_hours(separations, mp, rp)
    veq = float(equatorial_velocity_km_s(rp, prot))

    planet = PlanetConfig(v_eq=veq, line_sigma=1.5, line_depth=1.0)
    atmosphere = AtmosphereConfig(
        temperature=teq,
        mean_molecular_weight=2.3,
        effective_scale_heights=5.0,
        line_contrast=1.0,
        planet_mass_jupiter=mp,
        planet_radius_jupiter=rp,
        stellar_radius_solar=rstar,
    )
    sim = SimulationConfig(transit_duration=transit_duration, **sim_base)

    for idx, am in enumerate(separations):
        moon = MoonConfig(
            radius=rm,
            semi_major_axis=float(am),
            period=float(moon_period[idx]),
            phase=0.08,
            inclination=90.0,
        )
        result = run_xds_model(planet, moon, sim)
        report = detectability_report(result, atmosphere, observation, instrument)
        cadence = cadence_report(result, instrument)
        snr[idx] = report.matched_filter_snr
        contact[idx] = cadence.min_contact_duration
        exposure_limit[idx] = 60.0 * report.cadence_limited_exposure
        peak_ppm[idx] = 1.0e6 * report.peak_xds_depth

    np.savez(
        output_dir / f"xds_best_candidate_moon_separation_scan_{tag}.npz",
        moon_separation_rp=separations,
        moon_period_hours=moon_period,
        moon_period_over_transit=moon_period / transit_duration,
        snr=snr,
        min_contact_duration_hours=contact,
        exposure_limit_minutes=exposure_limit,
        peak_xds_ppm=peak_ppm,
        roche_limit_rp=roche,
        outer_stability_limit_rp=outer,
    )

    fig, axes = plt.subplots(4, 1, figsize=(8.8, 9.5), sharex=True, constrained_layout=True)
    axes[0].plot(separations, snr, color="black", linewidth=2.0)
    axes[0].axhline(5.0, color="tab:red", linestyle="--", linewidth=1.0)
    axes[0].set_ylabel("SNR")
    axes[0].set_title("Moon separation scan for the best sampled candidate")

    axes[1].plot(separations, peak_ppm, color="tab:blue", linewidth=2.0)
    axes[1].set_ylabel("Peak XDS [ppm]")

    axes[2].plot(separations, contact, color="tab:green", linewidth=2.0, label="contact")
    axes[2].set_ylabel("Contact duration [h]")
    ax2 = axes[2].twinx()
    ax2.plot(separations, moon_period / transit_duration, color="tab:orange", linewidth=1.6, label="Pm/Tdur")
    ax2.set_ylabel("Pm / Tdur")

    axes[3].plot(separations, exposure_limit, color="tab:purple", linewidth=2.0)
    axes[3].set_ylabel("Exposure limit [min]")
    axes[3].set_xlabel("Moon separation [Rp]")

    for ax in axes:
        ax.axvline(roche, color="0.25", linestyle=":", linewidth=1.0)
        ax.axvline(outer, color="0.25", linestyle=":", linewidth=1.0)
        ax.axvspan(0.0, roche, color="tab:red", alpha=0.08)
        ax.axvspan(outer, separations[-1], color="tab:red", alpha=0.08)

    fig.savefig(output_dir / f"xds_best_candidate_moon_separation_scan_{tag}.png", dpi=180)


if __name__ == "__main__":
    main()
