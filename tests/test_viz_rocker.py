import numpy as np

from crane.viz import rocker_joints


def test_rocker_contact_at_ground():
    """stance 円弧足の接触点が地面 (y=0) にある。"""
    pts = rocker_joints(np.array([0.2, -0.2, 0, 0]), L=1.0, R=0.3, contact_x=0.0)
    # pts = (stance_contact, stance_C, hip, swing_C, swing_lowest)
    assert np.isclose(pts[0][1], 0.0, atol=1e-12)


def test_rocker_curvature_center_height_R():
    """stance 曲率中心 C は高さ R。"""
    pts = rocker_joints(np.array([0.2, -0.2, 0, 0]), L=1.0, R=0.3, contact_x=0.0)
    assert np.isclose(pts[1][1], 0.3, atol=1e-12)
