from xds.dynamics import (
    breakup_rotation_period_hours,
    equatorial_velocity_km_s,
    moon_period_hours,
    planet_density_g_cm3,
)


def test_moon_period_increases_with_separation():
    assert moon_period_hours(5.0, 1.0, 1.0) > moon_period_hours(2.0, 1.0, 1.0)


def test_equatorial_velocity_increases_for_faster_spin():
    assert equatorial_velocity_km_s(1.0, 5.0) > equatorial_velocity_km_s(1.0, 10.0)


def test_breakup_period_is_positive():
    assert breakup_rotation_period_hours(1.0, 1.0) > 0.0


def test_density_decreases_for_inflated_planet():
    assert planet_density_g_cm3(1.0, 1.5) < planet_density_g_cm3(1.0, 1.0)
