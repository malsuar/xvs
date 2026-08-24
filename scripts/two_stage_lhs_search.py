"""Two-stage XDS search: Latin hypercube screening plus forward-model refinement."""

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
from xds.observability import cadence_report
from xds.sampling import latin_hypercube_candidates
from xds.screening import fast_screen_candidates, select_top_indices
from xds.spectra import run_xds_model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-screen", type=int, default=100_000)
    parser.add_argument("--n-refine", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--tag", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tag = args.tag or f"lhs{args.n_screen}_ref{args.n_refine}_seed{args.seed}"
    root = Path(__file__).resolve().parents[1]
    output_dir = root / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

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

    samples = latin_hypercube_candidates(args.n_screen, seed=args.seed)
    screened = fast_screen_candidates(samples, observation)
    selected = select_top_indices(screened, args.n_refine)
    if selected.size == 0:
        raise RuntimeError("No valid candidates survived fast screening.")

    screen_payload = {
        key: value
        for key, value in screened.items()
        if key in samples or key in {"valid", "fast_score", "density_g_cm3", "temperature_k", "v_eq_km_s"}
    }
    screen_payload["selected_index"] = selected
    np.savez(output_dir / f"xds_lhs_screen_{tag}.npz", **screen_payload)

    records: list[dict[str, float]] = []
    for idx in selected:
        mp = float(screened["planet_mass_jupiter"][idx])
        rp = float(screened["planet_radius_jupiter"][idx])
        rstar = float(screened["stellar_radius_solar"][idx])
        prot = float(screened["rotation_period_hours"][idx])
        am = float(screened["moon_separation_rp"][idx])
        rm = float(screened["moon_radius_rp"][idx])
        teq = float(screened["temperature_k"][idx])
        veq = float(screened["v_eq_km_s"][idx])
        pm_hours = float(screened["moon_period_hours"][idx])
        transit_duration = float(screened["transit_duration_hours"][idx])

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
                "fast_score": float(screened["fast_score"][idx]),
                "peak_xds_ppm": 1.0e6 * report.peak_xds_depth,
                "line_depth_ppm": 1.0e6 * report.atmospheric_line_depth,
                "scale_height_km": report.scale_height_km,
                "density_g_cm3": float(screened["density_g_cm3"][idx]),
                "v_eq_km_s": veq,
                "planet_mass_jupiter": mp,
                "planet_radius_jupiter": rp,
                "stellar_radius_solar": rstar,
                "planet_period_days": float(screened["planet_period_days"][idx]),
                "rotation_period_hours": prot,
                "breakup_period_hours": float(screened["breakup_period_hours"][idx]),
                "moon_separation_rp": am,
                "moon_period_hours": pm_hours,
                "moon_radius_rp": rm,
                "roche_limit_rp": float(screened["roche_limit_rp"][idx]),
                "outer_stability_limit_rp": float(screened["outer_stability_limit_rp"][idx]),
                "transit_duration_hours": transit_duration,
                "temperature_k": teq,
                "min_contact_duration_hours": cadence.min_contact_duration,
                "integration_limit_minutes": 60.0 * report.cadence_limited_exposure,
                "smearing_penalty": report.smearing_penalty,
            }
        )

    names = list(records[0])
    arrays = {name: np.array([record[name] for record in records], dtype=float) for name in names}
    np.savez(output_dir / f"xds_lhs_refined_{tag}.npz", **arrays)

    order = np.argsort(arrays["snr"])[::-1]
    print(f"screened={args.n_screen} valid={int(np.sum(screened['valid']))} refined={len(records)}")
    print("Top refined XDS candidates")
    for rank, row in enumerate(order[:12], start=1):
        print(
            f"{rank:02d} SNR={arrays['snr'][row]:.2f} "
            f"peak={arrays['peak_xds_ppm'][row]:.1f} ppm "
            f"Mp={arrays['planet_mass_jupiter'][row]:.3f} Mj "
            f"Rp={arrays['planet_radius_jupiter'][row]:.2f} Rj "
            f"rho={arrays['density_g_cm3'][row]:.3f} "
            f"Teq={arrays['temperature_k'][row]:.0f} K "
            f"veq={arrays['v_eq_km_s'][row]:.1f} km/s "
            f"am={arrays['moon_separation_rp'][row]:.2f} Rp "
            f"Pm/Tdur={arrays['moon_period_hours'][row]/arrays['transit_duration_hours'][row]:.1f}"
        )


if __name__ == "__main__":
    main()
