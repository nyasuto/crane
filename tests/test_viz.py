import numpy as np

from crane.viz import link_points


def test_link_points_heelstrike_symmetry():
    """heel-strike 瞬間 (φ=2θ) には swing 足が接地高さ (y=0) にある。"""
    theta = -0.2
    hip, swing = link_points(theta, 2 * theta, foot_x=0.0)
    assert np.isclose(hip[1], np.cos(theta))
    assert np.isclose(swing[1], 0.0, atol=1e-12)


def test_link_points_upright():
    """θ=0, φ=0 で両脚とも鉛直、swing 足は stance 足と重なる。"""
    hip, swing = link_points(0.0, 0.0, foot_x=1.5)
    assert np.allclose(hip, [1.5, 1.0])
    assert np.allclose(swing, [1.5, 0.0])


def test_link_points_abs_matches_relative_convention():
    """絶対角版と (θ, φ) 版の整合。"""
    from crane.viz import link_points, link_points_abs

    hip1, swing1 = link_points(0.2, 0.5, foot_x=1.0)
    hip2, swing2 = link_points_abs(0.2, 0.2 - 0.5, foot_x=1.0)
    assert np.allclose(hip1, hip2)
    assert np.allclose(swing1, swing2)
