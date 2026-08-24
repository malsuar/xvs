"""Generate a Fig. 2-style XDS time-velocity map."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from xds.config import MoonConfig, PlanetConfig, SimulationConfig
from xds.diagnostics import map_diagnostics
from xds.spectra import run_xds_model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=["fiducial", "lhs-best", "paper-showcase"], default="fiducial")
    parser.add_argument("--tag", default=None)
    parser.add_argument("--window-hours", type=float, default=None)
    parser.add_argument("--time-start", type=float, default=None)
    parser.add_argument("--time-end", type=float, default=None)
    parser.add_argument("--center-mode", choices=["midtransit", "peak"], default="midtransit")
    return parser.parse_args()


def build_configs(preset: str) -> tuple[PlanetConfig, MoonConfig, SimulationConfig, str]:
    if preset == "fiducial":
        planet = PlanetConfig(v_eq=12.0, line_sigma=1.5, line_depth=1.0)
        moon = MoonConfig(radius=0.28, semi_major_axis=3.0, period=5.0, phase=0.08)
        sim = SimulationConfig(
            transit_duration=5.0,
            n_times=241,
            n_longitudes=181,
            n_latitudes=121,
            n_velocities=501,
            velocity_min=-25.0,
            velocity_max=25.0,
            use_illumination_shadow=True,
            use_observer_occultation=True,
        )
        return planet, moon, sim, "fiducial"

    if preset == "lhs-best":
        planet = PlanetConfig(v_eq=10.052449786472337, line_sigma=1.5, line_depth=1.0)
        moon = MoonConfig(
            radius=0.34237930871226385,
            semi_major_axis=2.1680725360270703,
            period=58.950105370121975,
            phase=0.08,
        )
        sim = SimulationConfig(
            transit_duration=2.862057477656129,
            n_times=281,
            n_longitudes=181,
            n_latitudes=121,
            n_velocities=501,
            velocity_min=-25.0,
            velocity_max=25.0,
            use_illumination_shadow=True,
            use_observer_occultation=True,
        )
        return planet, moon, sim, "lhs_best"

    if preset == "paper-showcase":
        planet = PlanetConfig(v_eq=12.0, line_sigma=1.5, line_depth=1.0)
        moon = MoonConfig(
            radius=0.0892,
            semi_major_axis=3.0,
            period=21.238292071543107,
            phase=0.08,
        )
        sim = SimulationConfig(
            transit_duration=16.0,
            n_times=281,
            n_longitudes=181,
            n_latitudes=121,
            n_velocities=501,
            velocity_min=-25.0,
            velocity_max=25.0,
            use_illumination_shadow=True,
            use_observer_occultation=True,
        )
        return planet, moon, sim, "paper_showcase"

    raise ValueError(f"unknown preset: {preset}")


def with_sim_window(sim: SimulationConfig, window_hours: float | None) -> SimulationConfig:
    if window_hours is None:
        return sim
    return SimulationConfig(
        transit_duration=window_hours,
        n_times=sim.n_times,
        n_longitudes=sim.n_longitudes,
        n_latitudes=sim.n_latitudes,
        velocity_min=sim.velocity_min,
        velocity_max=sim.velocity_max,
        n_velocities=sim.n_velocities,
        star_direction=sim.star_direction,
        observer_direction=sim.observer_direction,
        use_illumination_shadow=sim.use_illumination_shadow,
        use_observer_occultation=sim.use_observer_occultation,
        include_mirrored_terminator=sim.include_mirrored_terminator,
    )


def recenter_moon_on_peak(planet: PlanetConfig, moon: MoonConfig, sim: SimulationConfig) -> MoonConfig:
    scan_sim = SimulationConfig(
        transit_duration=moon.period,
        n_times=721,
        n_longitudes=max(61, sim.n_longitudes // 2),
        n_latitudes=max(41, sim.n_latitudes // 2),
        velocity_min=sim.velocity_min,
        velocity_max=sim.velocity_max,
        n_velocities=101,
        star_direction=sim.star_direction,
        observer_direction=sim.observer_direction,
        use_illumination_shadow=sim.use_illumination_shadow,
        use_observer_occultation=sim.use_observer_occultation,
        include_mirrored_terminator=sim.include_mirrored_terminator,
    )
    scan = run_xds_model(planet, moon, scan_sim)
    idx = int(np.argmax(scan.combined_fraction))
    peak_time = float(scan.time[idx])
    phase_shift = peak_time / moon.period
    return MoonConfig(
        radius=moon.radius,
        semi_major_axis=moon.semi_major_axis,
        period=moon.period,
        phase=(moon.phase + phase_shift) % 1.0,
        inclination=moon.inclination,
        longitude_node=moon.longitude_node,
        argument_pericenter=moon.argument_pericenter,
        eccentricity=moon.eccentricity,
    )


def shift_moon_phase(moon: MoonConfig, time_offset_hours: float) -> MoonConfig:
    return MoonConfig(
        radius=moon.radius,
        semi_major_axis=moon.semi_major_axis,
        period=moon.period,
        phase=(moon.phase + time_offset_hours / moon.period) % 1.0,
        inclination=moon.inclination,
        longitude_node=moon.longitude_node,
        argument_pericenter=moon.argument_pericenter,
        eccentricity=moon.eccentricity,
    )


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    planet, moon, sim, default_tag = build_configs(args.preset)
    if (args.time_start is None) != (args.time_end is None):
        raise ValueError("--time-start and --time-end must be supplied together")
    if args.time_start is not None:
        if args.time_end <= args.time_start:
            raise ValueError("--time-end must be larger than --time-start")
        window_hours = args.time_end - args.time_start
        time_offset = 0.5 * (args.time_start + args.time_end)
    else:
        window_hours = args.window_hours
        time_offset = 0.0
    sim = with_sim_window(sim, window_hours)
    if args.center_mode == "peak":
        moon = recenter_moon_on_peak(planet, moon, sim)
    if time_offset != 0.0:
        moon = shift_moon_phase(moon, time_offset)
    tag = args.tag or default_tag
    output = root / "outputs" / f"fig2_xds_maps_{tag}.png"
    data_output = root / "outputs" / f"fig2_xds_maps_{tag}.npz"
    output.parent.mkdir(parents=True, exist_ok=True)

    result = run_xds_model(planet, moon, sim)
    plot_time = result.time + time_offset
    diagnostics = map_diagnostics(result.time, result.velocity, result.delta_absorption)
    np.savez(
        data_output,
        time=plot_time,
        model_time=result.time,
        velocity=result.velocity,
        reference_absorption=result.reference_absorption,
        absorption=result.absorption,
        delta_absorption=result.delta_absorption,
        moon_positions=result.moon_positions,
        illumination_fraction=result.illumination_fraction,
        observer_fraction=result.observer_fraction,
        combined_fraction=result.combined_fraction,
        centroid=diagnostics.centroid,
        skewness=diagnostics.skewness,
        blue_signal=diagnostics.blue_signal,
        red_signal=diagnostics.red_signal,
        shadow_velocity=diagnostics.shadow_velocity,
        shadow_strength=diagnostics.shadow_strength,
        planet_v_eq=planet.v_eq,
        moon_radius=moon.radius,
        moon_semi_major_axis=moon.semi_major_axis,
        moon_period=moon.period,
        moon_phase=moon.phase,
        transit_duration=sim.transit_duration,
        time_offset=time_offset,
    )

    extent = [
        result.velocity[0],
        result.velocity[-1],
        plot_time[0],
        plot_time[-1],
    ]
    vmax = np.max(result.reference_absorption)
    dmax = np.max(np.abs(result.delta_absorption))

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(15.5, 4.2),
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.0, 0.9]},
        constrained_layout=True,
    )

    panels = [
        (result.reference_absorption, "Moonless planet", "inferno", 0.0, vmax, r"Absorption $A_0(t,v)$"),
        (result.absorption, "Planet + moon", "inferno", 0.0, vmax, r"Absorption $A(t,v)$"),
        (result.delta_absorption, "Differential absorption", "RdBu_r", -dmax, dmax, r"Residual $\Delta A(t,v)$"),
    ]
    for ax, (data, title, cmap, vmin, vmax_panel, cbar_label) in zip(axes[:3], panels):
        im = ax.imshow(
            data,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax_panel,
        )
        ax.set_title(title)
        ax.set_xlabel("Velocity [km/s]")
        ax.set_ylabel("Time [h]")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046)
        cbar.set_label(cbar_label)

    ax = axes[3]
    ax.plot(plot_time, diagnostics.centroid, label="centroid")
    ax.plot(plot_time, diagnostics.shadow_velocity, label="shadow velocity")
    ax.plot(plot_time, 12.0 * result.illumination_fraction, label="illumination x12")
    ax.plot(plot_time, 12.0 * result.observer_fraction, label="observer x12")
    ax.axhline(0.0, color="0.4", linewidth=0.8)
    ax.set_title("Diagnostics")
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("Velocity [km/s] or scaled fraction")
    ax.legend(fontsize=8)

    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
