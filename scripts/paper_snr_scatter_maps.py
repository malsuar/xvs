"""Create clean paper-style scatter maps from an XDS search."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="n5000_seed7")
    parser.add_argument("--input", default=None, help="Optional .npz file to plot.")
    parser.add_argument("--output", default=None, help="Optional output PNG path.")
    return parser.parse_args()


def draw_panel(ax, x, y, snr, top, xlabel, ylabel, title, xscale=None, xlim=None, ylim=None):
    low = snr < 5.0
    high = snr >= 5.0
    ax.scatter(x[low], y[low], s=10, c="#c9c9c9", alpha=0.35, linewidths=0.0, rasterized=True)
    sc = ax.scatter(
        x[high],
        y[high],
        s=34,
        c=snr[high],
        cmap="inferno",
        vmin=5.0,
        vmax=27.0,
        edgecolors="black",
        linewidths=0.25,
    )
    ax.scatter(
        x[top],
        y[top],
        s=64,
        facecolors="none",
        edgecolors="#19bff2",
        linewidths=1.2,
    )
    if xscale:
        ax.set_xscale(xscale)
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.axhline(0, color="none")
    return sc


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    data_path = Path(args.input) if args.input else root / "outputs" / f"xds_sweet_spot_search_{args.tag}.npz"
    output_path = Path(args.output) if args.output else root / "outputs" / f"xds_paper_snr_scatter_maps_{args.tag}.png"
    if not data_path.exists():
        raise FileNotFoundError(data_path)

    data = np.load(data_path)
    snr = data["snr"]
    detectable = np.flatnonzero(snr >= 5.0)
    top = detectable[np.argsort(snr[detectable])[::-1]][:20]
    if "density_g_cm3" in data.files:
        density = data["density_g_cm3"]
    else:
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
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)
    sc = draw_panel(
        axes[0, 0],
        data["planet_mass_jupiter"],
        data["planet_radius_jupiter"],
        snr,
        top,
        r"$M_p$ [$M_J$]",
        r"$R_p$ [$R_J$]",
        "Inflated low-mass hosts",
        xscale="log",
        xlim=(0.07, 2.2),
        ylim=(0.72, 1.78),
    )
    draw_panel(
        axes[0, 1],
        density,
        data["scale_height_km"],
        snr,
        top,
        r"$\rho_p$ [g cm$^{-3}$]",
        r"$H$ [km]",
        "Atmospheric column",
        xscale="log",
        xlim=(0.03, 2.7),
        ylim=(0, 1800),
    )
    draw_panel(
        axes[1, 0],
        data["temperature_k"],
        data["v_eq_km_s"],
        snr,
        top,
        r"$T_{\rm eq}$ [K]",
        r"$v_{\rm eq}$ [km s$^{-1}$]",
        "Thermal-Doppler balance",
        xlim=(240, 1060),
        ylim=(3, 42),
    )
    draw_panel(
        axes[1, 1],
        data["moon_separation_rp"],
        pm_over_tdur,
        snr,
        top,
        r"$a_m/R_p$",
        r"$P_m/T_{\rm dur}$",
        "Satellite timing",
        xlim=(1.6, 12.2),
        ylim=(0, 55),
    )

    cbar = fig.colorbar(sc, ax=axes, shrink=0.92, pad=0.02)
    cbar.set_label("XDS SNR for detectable candidates")
    fig.savefig(output_path, dpi=240)
    print(output_path)


if __name__ == "__main__":
    main()
