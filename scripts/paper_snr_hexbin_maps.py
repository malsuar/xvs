"""Create cleaner paper-style SNR hexbin maps from an XDS search."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="n5000_seed7")
    parser.add_argument("--gridsize", type=int, default=26)
    return parser.parse_args()


def finite_limits(values: np.ndarray, pad: float = 0.03) -> tuple[float, float]:
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    span = hi - lo
    return lo - pad * span, hi + pad * span


def draw_hex_panel(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    snr: np.ndarray,
    top: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
    gridsize: int,
    xscale: str | None = None,
    extent: tuple[float, float, float, float] | None = None,
):
    if xscale == "log":
        x_plot = np.log10(x)
        top_x = np.log10(x[top])
        if extent is None:
            xlim = finite_limits(x_plot)
            ylim = finite_limits(y)
            extent = (xlim[0], xlim[1], ylim[0], ylim[1])
        hb = ax.hexbin(
            x_plot,
            y,
            C=snr,
            reduce_C_function=np.max,
            gridsize=gridsize,
            mincnt=1,
            extent=extent,
            cmap="inferno",
            vmin=0.0,
            vmax=25.0,
            linewidths=0.0,
        )
        ax.scatter(top_x, y[top], s=24, facecolors="none", edgecolors="#49c7f4", linewidths=1.0)
        ticks = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0])
        ticks = ticks[(ticks >= np.nanmin(x)) & (ticks <= np.nanmax(x))]
        ax.set_xticks(np.log10(ticks), [f"{tick:g}" for tick in ticks])
    else:
        if extent is None:
            xlim = finite_limits(x)
            ylim = finite_limits(y)
            extent = (xlim[0], xlim[1], ylim[0], ylim[1])
        hb = ax.hexbin(
            x,
            y,
            C=snr,
            reduce_C_function=np.max,
            gridsize=gridsize,
            mincnt=1,
            extent=extent,
            cmap="inferno",
            vmin=0.0,
            vmax=25.0,
            linewidths=0.0,
        )
        ax.scatter(x[top], y[top], s=24, facecolors="none", edgecolors="#49c7f4", linewidths=1.0)

    ax.contour(
        *hexbin_to_grid(hb),
        levels=[5.0, 10.0, 20.0],
        colors="white",
        linewidths=0.8,
        alpha=0.9,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_facecolor("#f3f3f3")
    return hb


def hexbin_to_grid(hb):
    offsets = hb.get_offsets()
    values = hb.get_array()
    if len(values) < 8:
        return np.array([[0.0]]), np.array([[0.0]]), np.array([[np.nan]])
    x_unique = np.unique(np.round(offsets[:, 0], 10))
    y_unique = np.unique(np.round(offsets[:, 1], 10))
    grid = np.full((len(y_unique), len(x_unique)), np.nan)
    x_index = {value: idx for idx, value in enumerate(x_unique)}
    y_index = {value: idx for idx, value in enumerate(y_unique)}
    for (xval, yval), zval in zip(np.round(offsets, 10), values):
        grid[y_index[yval], x_index[xval]] = zval
    xx, yy = np.meshgrid(x_unique, y_unique)
    return xx, yy, grid


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    data_path = root / "outputs" / f"xds_sweet_spot_search_{args.tag}.npz"
    output_path = root / "outputs" / f"xds_paper_snr_hexbin_maps_{args.tag}.png"
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
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.dpi": 140,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.1), constrained_layout=True)

    hb = draw_hex_panel(
        axes[0, 0],
        data["planet_mass_jupiter"],
        data["planet_radius_jupiter"],
        snr,
        top,
        r"$M_p$ [$M_J$]",
        r"$R_p$ [$R_J$]",
        "Inflated low-mass hosts",
        args.gridsize,
        xscale="log",
        extent=(np.log10(0.08), np.log10(2.0), 0.75, 1.75),
    )
    draw_hex_panel(
        axes[0, 1],
        density,
        data["scale_height_km"],
        snr,
        top,
        r"$\rho_p$ [g cm$^{-3}$]",
        r"$H$ [km]",
        "Atmospheric column",
        args.gridsize,
        xscale="log",
        extent=(np.log10(0.035), np.log10(2.5), 0.0, 1800.0),
    )
    draw_hex_panel(
        axes[1, 0],
        data["temperature_k"],
        data["v_eq_km_s"],
        snr,
        top,
        r"$T_{\rm eq}$ [K]",
        r"$v_{\rm eq}$ [km s$^{-1}$]",
        "Thermal-Doppler balance",
        args.gridsize,
        extent=(250.0, 1050.0, 3.0, 42.0),
    )
    draw_hex_panel(
        axes[1, 1],
        data["moon_separation_rp"],
        pm_over_tdur,
        snr,
        top,
        r"$a_m/R_p$",
        r"$P_m/T_{\rm dur}$",
        "Satellite timing",
        args.gridsize,
        extent=(1.7, 12.0, 0.0, 55.0),
    )

    cbar = fig.colorbar(hb, ax=axes, shrink=0.92, pad=0.02)
    cbar.set_label("Maximum XDS SNR")
    fig.savefig(output_path, dpi=240)
    print(output_path)


if __name__ == "__main__":
    main()
