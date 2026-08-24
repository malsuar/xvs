from xds.config import AtmosphereConfig
from xds.detectability import atmospheric_line_depth, scale_height_km


def test_scale_height_increases_with_temperature():
    cold = AtmosphereConfig(temperature=300.0)
    hot = AtmosphereConfig(temperature=900.0)
    assert scale_height_km(hot) > scale_height_km(cold)


def test_line_depth_increases_with_scale_heights():
    shallow = AtmosphereConfig(effective_scale_heights=2.0)
    deep = AtmosphereConfig(effective_scale_heights=8.0)
    assert atmospheric_line_depth(deep) > atmospheric_line_depth(shallow)
