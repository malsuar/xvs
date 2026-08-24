"""Create an appendix figure for the EXOVEL terminator geometry."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from xds.config import MoonConfig, PlanetConfig, SimulationConfig
from xds.moon_orbit import moon_positions
from xds.rotation import line_of_sight_velocity
from xds.shadow import compute_shadow_masks
from xds.spectra import gaussian_profiles, weighted_absorption
from xds.terminator import build_terminator


def selected_mask_examples(
    planet: PlanetConfig,
    moon: MoonConfig,
    sim: SimulationConfig,
) -> tuple[tuple[float, np.ndarray, object], tuple[float, np.ndarray, object]]:
    times = np.linspace(-0.5 * moon.period, 0.5 * moon.period, 721)
    positions = moon_positions(times, moon)
    terminator = build_terminator(
        sim.n_longitudes,
        sim.n_latitudes,
        include_mirror=sim.include_mirrored_terminator,
    )
    best_illum = (-np.inf, 0, None)
    best_observer = (-np.inf, 0, None)
    for idx, moon_pos in enumerate(positions):
        masks = compute_shadow_masks(
            terminator.positions,
            moon_pos,
            moon.radius,
            sim.star_direction,
            sim.observer_direction,
        )
        illum = np.sum(terminator.weights[masks.illumination])
        obs = np.sum(terminator.weights[masks.observer])
        if illum > best_illum[0]:
            best_illum = (float(illum), idx, masks)
        if obs > best_observer[0]:
            best_observer = (float(obs), idx, masks)
    assert best_illum[2] is not None
    assert best_observer[2] is not None
    illum_idx = best_illum[1]
    obs_idx = best_observer[1]
    return (
        (float(times[illum_idx]), positions[illum_idx], best_illum[2]),
        (float(times[obs_idx]), positions[obs_idx], best_observer[2]),
    )


def draw_mask_panel(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    blocked: np.ndarray,
    point_sizes: np.ndarray,
    color: str,
    label: str,
    title: str,
    time: float,
    xlabel: str,
    ylabel: str,
) -> None:
    ax.scatter(
        x,
        y,
        c="0.86",
        s=point_sizes * 0.65,
        linewidths=0,
        alpha=0.75,
        label="unmasked",
    )
    ax.scatter(
        x[blocked],
        y[blocked],
        c=color,
        s=24,
        label=label,
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    ax.text(
        0.02,
        0.04,
        rf"$t={time:.2f}$ h",
        transform=ax.transAxes,
        fontsize=8,
        color="0.25",
    )


def make_figure(rotated_masks: bool = False) -> Path:
    root = Path(__file__).resolve().parents[1]
    output_dir = root / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_rotated_masks" if rotated_masks else ""
    output = output_dir / f"figA1_exovel_terminator_geometry{suffix}.png"

    planet = PlanetConfig(v_eq=12.0, line_sigma=1.5, line_depth=1.0)
    moon = MoonConfig(radius=0.16, semi_major_axis=2.6, period=8.0, phase=0.11)
    sim = SimulationConfig(
        transit_duration=8.0,
        n_times=241,
        n_longitudes=101,
        n_latitudes=71,
        n_velocities=501,
        velocity_min=-25.0,
        velocity_max=25.0,
        use_illumination_shadow=True,
        use_observer_occultation=True,
        include_mirrored_terminator=True,
    )

    terminator = build_terminator(
        sim.n_longitudes,
        sim.n_latitudes,
        include_mirror=sim.include_mirrored_terminator,
    )
    v_los = line_of_sight_velocity(
        terminator.positions,
        planet.v_eq,
        planet.spin_inclination,
        planet.spin_obliquity,
        sim.observer_direction,
    )
    (illum_time, _illum_pos, illum_masks), (obs_time, _obs_pos, obs_masks) = selected_mask_examples(
        planet, moon, sim
    )

    sky_y = terminator.positions[:, 1]
    sky_z = terminator.positions[:, 2]
    weights = terminator.weights / np.nanmax(terminator.weights)
    point_sizes = 5.0 + 28.0 * weights

    velocity = np.linspace(sim.velocity_min, sim.velocity_max, sim.n_velocities)
    kernels = gaussian_profiles(velocity, v_los, planet.line_sigma)
    reference_kernel = weighted_absorption(kernels, terminator.weights, planet.line_depth)
    illum_weights = np.where(illum_masks.illumination, 0.0, terminator.weights)
    obs_weights = np.where(obs_masks.observer, 0.0, terminator.weights)
    illum_kernel = weighted_absorption(kernels, illum_weights, planet.line_depth)
    obs_kernel = weighted_absorption(kernels, obs_weights, planet.line_depth)

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(15.6, 3.65),
        gridspec_kw={"width_ratios": [1.05, 0.85, 0.85, 0.95]},
        constrained_layout=True,
    )

    first_x = sky_z
    first_y = sky_y
    first_xlabel = r"Sky coordinate $z/R_{\rm p}$"
    first_ylabel = r"Sky coordinate $y/R_{\rm p}$"

    ax = axes[0]
    sc = ax.scatter(
        first_x,
        first_y,
        c=v_los,
        s=point_sizes,
        cmap="RdBu_r",
        vmin=-planet.v_eq,
        vmax=planet.v_eq,
        linewidths=0,
        alpha=0.92,
    )
    ax.set_title("Projected terminator")
    ax.set_xlabel(first_xlabel)
    ax.set_ylabel(first_ylabel)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal", adjustable="box")
    ax.text(
        0.02,
        0.04,
        r"point size $\propto w_i^0$",
        transform=ax.transAxes,
        fontsize=8,
        color="0.25",
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046)
    cbar.set_label(r"$v_{\rm los}$ [km s$^{-1}$]")

    if rotated_masks:
        mask_x = sky_z
        mask_y = sky_y
        mask_xlabel = r"Sky coordinate $z/R_{\rm p}$"
        mask_ylabel = r"Sky coordinate $y/R_{\rm p}$"
    else:
        mask_x = sky_y
        mask_y = sky_z
        mask_xlabel = r"Sky coordinate $y/R_{\rm p}$"
        mask_ylabel = r"Sky coordinate $z/R_{\rm p}$"

    draw_mask_panel(
        axes[1],
        mask_x,
        mask_y,
        illum_masks.illumination,
        point_sizes,
        "#1f78b4",
        r"$B_{\rm ill}=1$",
        "Illumination blocking",
        illum_time,
        mask_xlabel,
        mask_ylabel,
    )
    draw_mask_panel(
        axes[2],
        mask_x,
        mask_y,
        obs_masks.observer,
        point_sizes,
        "#e66101",
        r"$B_{\rm obs}=1$",
        "Observer-line blocking",
        obs_time,
        mask_xlabel,
        mask_ylabel,
    )

    ax = axes[3]
    ax.plot(velocity, reference_kernel, color="black", lw=2.0, label=r"$K_0(v)$")
    ax.plot(velocity, illum_kernel, color="#1f78b4", lw=2.0, label=r"$K_{\rm ill}(v,t)$")
    ax.plot(velocity, obs_kernel, color="#e66101", lw=2.0, label=r"$K_{\rm obs}(v,t)$")
    ax.fill_between(
        velocity,
        illum_kernel,
        reference_kernel,
        where=reference_kernel > illum_kernel,
        color="#1f78b4",
        alpha=0.22,
    )
    ax.fill_between(
        velocity,
        obs_kernel,
        reference_kernel,
        where=reference_kernel > obs_kernel,
        color="#e66101",
        alpha=0.22,
    )
    ax.set_title("Kernel truncation")
    ax.set_xlabel(r"Velocity [km s$^{-1}$]")
    ax.set_ylabel("Normalised absorption")
    ax.set_xlim(-16, 16)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8, frameon=True)

    for ax in axes:
        ax.tick_params(direction="in", top=True, right=True)

    fig.savefig(output, dpi=300)
    fig.savefig(output.with_suffix(".pdf"))
    return output


def main() -> None:
    print(make_figure(rotated_masks=False))
    print(make_figure(rotated_masks=True))


if __name__ == "__main__":
    main()
