"""Generate a conceptual XDS efficiency curve."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from xds.efficiency import hill_radius_over_planet_radius, xds_efficiency_metric


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "outputs" / "fig3_efficiency.png"
    output.parent.mkdir(parents=True, exist_ok=True)

    period = np.logspace(0, 4, 600)
    metric = xds_efficiency_metric(period)
    hill = hill_radius_over_planet_radius(period)

    fig, ax = plt.subplots(figsize=(8.5, 4.2), constrained_layout=True)
    ax.plot(period, metric, color="black", linewidth=2.0)
    ax.fill_between(period, 0.0, metric, where=period < 10.0, alpha=0.18, color="tab:red", label="hot")
    ax.fill_between(
        period,
        0.0,
        metric,
        where=(period >= 10.0) & (period < 300.0),
        alpha=0.18,
        color="tab:green",
        label="warm",
    )
    ax.fill_between(period, 0.0, metric, where=period >= 300.0, alpha=0.18, color="tab:blue", label="cold")
    ax.set_xscale("log")
    ax.set_xlabel("Planet orbital period [days]")
    ax.set_ylabel("Normalised heuristic XDS efficiency")
    ax.legend()

    top = ax.secondary_xaxis(
        "top",
        functions=(
            lambda p: hill_radius_over_planet_radius(p),
            lambda h: np.interp(h, hill, period),
        ),
    )
    top.set_xlabel("Approximate Hill radius [Rp]")

    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
