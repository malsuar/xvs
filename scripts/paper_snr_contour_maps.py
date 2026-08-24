"""Create paper-style SNR maps from an XDS sweet-spot search."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="n5000_seed7")
    parser.add_argument("--bins", type=int, default=24)
    return parser.parse_args()


def binned_max(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> np.ndarray:
    grid = np.full((y_edges.size - 1, x_edges.size - 1), np.nan)
    x_idx = np.searchsorted(x_edges, x, side="right") - 1
    y_idx = np.searchsorted(y_edges, y, side="right") - 1
    valid = (x_idx >= 0) & (x_idx < grid.shape[1]) & (y_idx >= 0) & (y_idx < grid.shape[0])
    for ix, iy, value in zip(x_idx[valid], y_idx[valid], z[valid]):
        if np.isnan(grid[iy, ix]) or value > grid[iy, ix]:
            grid[iy, ix] = value
    return grid


def centers(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def plot_panel(ax, x, y, snr, x_edges, y_edges, xlabel, ylabel, xscale="linear", top=None):
    grid = binned_max(x, y, snr, x_edges, y_edges)
    masked_grid = np.ma.masked_invalid(grid)
    xx, yy = np.meshgrid(centers(x_edges), centers(y_edges))
    cmap = __import__("matplotlib").colormaps["inferno"].copy()
    cmap.set_bad("#eeeeee")
    image = ax.pcolormesh(x_edges, y_edges, masked_grid, shading="auto", cmap=cmap, vmin=0.0, vmax=25.0)
    levels = [5.0, 10.0, 20.0]
    if np.nanmax(grid) >= min(levels):
        valid_levels = [level for level in levels if np.nanmin(grid) <= level <= np.nanmax(grid)]
        if valid_levels:
            cs = ax.contour(xx, yy, masked_grid, levels=valid_levels, colors="white", linewidths=1.0)
            ax.clabel(cs, inline=True, fmt=lambda value: f"{value:g}", fontsize=8)
    if top is not None:
        ax.scatter(x[top], y[top], s=18, facecolors="none", edgecolors="#61dafb", linewidths=0.9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    return image


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    data_path = root / "outputs" / f"xds_sweet_spot_search_{args.tag}.npz"
    output_path = root / "outputs" / f"xds_paper_snr_contour_maps_{args.tag}.png"
    if not data_path.exists():
        raise FileNotFoundError(data_path)

    data = np.load(data_path)
    snr = data["snr"]
    top = np.argsort(snr)[::-1][:20]
    density = 1.326 * data["planet_mass_jupiter"] / data["planet_radius_jupiter"] ** 3
    pm_over_tdur = data["moon_period_hours"] / data["transit_duration_hours"]

    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.4), constrained_layout=True)

    image = plot_panel(
        axes[0, 0],
        data["planet_mass_jupiter"],
        data["planet_radius_jupiter"],
        snr,
        np.logspace(np.log10(0.08), np.log10(2.0), args.bins),
        np.linspace(0.75, 1.75, args.bins),
        r"$M_p$ [$M_J$]",
        r"$R_p$ [$R_J$]",
        xscale="log",
        top=top,
    )
    axes[0, 0].set_title("Inflated low-mass hosts")

    plot_panel(
        axes[0, 1],
        density,
        data["scale_height_km"],
        snr,
        np.logspace(np.log10(0.035), np.log10(2.5), args.bins),
        np.linspace(0.0, 1800.0, args.bins),
        r"$\rho_p$ [g cm$^{-3}$]",
        r"$H$ [km]",
        xscale="log",
        top=top,
    )
    axes[0, 1].set_title("Atmospheric column")

    plot_panel(
        axes[1, 0],
        data["temperature_k"],
        data["v_eq_km_s"],
        snr,
        np.linspace(250.0, 1050.0, args.bins),
        np.linspace(3.0, 42.0, args.bins),
        r"$T_{\rm eq}$ [K]",
        r"$v_{\rm eq}$ [km s$^{-1}$]",
        top=top,
    )
    axes[1, 0].set_title("Thermal-Doppler balance")

    plot_panel(
        axes[1, 1],
        data["moon_separation_rp"],
        pm_over_tdur,
        snr,
        np.linspace(1.7, 12.0, args.bins),
        np.linspace(0.0, 55.0, args.bins),
        r"$a_m/R_p$",
        r"$P_m/T_{\rm dur}$",
        top=top,
    )
    axes[1, 1].set_title("Satellite timing")

    cbar = fig.colorbar(image, ax=axes, shrink=0.92, pad=0.02)
    cbar.set_label("Maximum XDS SNR in bin")
    fig.savefig(output_path, dpi=220)
    print(output_path)


if __name__ == "__main__":
    main()
