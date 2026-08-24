"""Fast physical screening for XDS sweet-spot searches."""

from __future__ import annotations

import numpy as np

from .config import AtmosphereConfig, ObservationConfig
from .detectability import atmospheric_line_depth, scale_height_km
from .dynamics import (
    breakup_rotation_period_hours,
    equatorial_velocity_km_s,
    fluid_roche_limit_rp,
    moon_period_hours,
    planet_density_g_cm3,
    stable_outer_moon_limit_rp,
    transit_duration_hours,
)


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


def fast_screen_candidates(
    samples: dict[str, np.ndarray],
    observation: ObservationConfig,
    star_mass_solar: float = 0.8,
    moon_density_g_cm3: float = 3.0,
) -> dict[str, np.ndarray]:
    """Filter candidates and compute an approximate ranking score."""

    mp = samples["planet_mass_jupiter"]
    rp = samples["planet_radius_jupiter"]
    rstar = samples["stellar_radius_solar"]
    pplanet = samples["planet_period_days"]
    prot = samples["rotation_period_hours"]
    am = samples["moon_separation_rp"]
    rm = samples["moon_radius_rp"]

    density = planet_density_g_cm3(mp, rp)
    roche = fluid_roche_limit_rp(density, moon_density_g_cm3=moon_density_g_cm3)
    outer = stable_outer_moon_limit_rp(pplanet, mp, star_mass_solar=star_mass_solar, planet_radius_jupiter=rp)
    breakup = breakup_rotation_period_hours(mp, rp)
    tdur = transit_duration_hours(pplanet, rstar, star_mass_solar=star_mass_solar)
    teq = equilibrium_temperature(pplanet, rstar, star_mass_solar=star_mass_solar)
    veq = equatorial_velocity_km_s(rp, prot)
    pm = moon_period_hours(am, mp, rp)

    valid = (
        (am > np.maximum(roche, 1.0 + rm))
        & (am < outer)
        & (prot > 1.15 * breakup)
        & np.isfinite(tdur)
        & (tdur >= 1.5)
        & (teq >= 250.0)
    )

    height = np.empty_like(teq)
    line_depth = np.empty_like(teq)
    for idx in range(teq.size):
        atmosphere = AtmosphereConfig(
            temperature=float(teq[idx]),
            mean_molecular_weight=2.3,
            effective_scale_heights=5.0,
            line_contrast=1.0,
            planet_mass_jupiter=float(mp[idx]),
            planet_radius_jupiter=float(rp[idx]),
            stellar_radius_solar=float(rstar[idx]),
        )
        height[idx] = scale_height_km(atmosphere)
        line_depth[idx] = atmospheric_line_depth(atmosphere)

    moon_area = rm**2
    doppler = veq / np.sqrt(veq**2 + 1.5**2)
    timing = np.exp(-0.5 * ((pm / np.clip(tdur, 1e-6, None) - 8.0) / 9.0) ** 2)
    proximity = np.exp(-0.35 * np.maximum(am - 2.0, 0.0))
    score = line_depth * moon_area * doppler * timing * proximity
    score *= np.sqrt(observation.n_lines * observation.n_transits)
    score = np.where(valid, score, -np.inf)

    out = dict(samples)
    out.update(
        {
            "valid": valid,
            "fast_score": score,
            "density_g_cm3": density,
            "roche_limit_rp": roche,
            "outer_stability_limit_rp": outer,
            "breakup_period_hours": breakup,
            "transit_duration_hours": np.minimum(tdur, 16.0),
            "temperature_k": teq,
            "v_eq_km_s": veq,
            "moon_period_hours": pm,
            "scale_height_km": height,
            "line_depth_ppm": 1.0e6 * line_depth,
        }
    )
    return out


def select_top_indices(screened: dict[str, np.ndarray], n_refine: int) -> np.ndarray:
    score = screened["fast_score"]
    valid = np.flatnonzero(np.isfinite(score))
    if valid.size == 0:
        return valid
    order = valid[np.argsort(score[valid])[::-1]]
    return order[: min(n_refine, order.size)]
