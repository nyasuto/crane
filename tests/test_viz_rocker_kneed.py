import numpy as np

from crane.viz import rocker_kneed_joints


def test_stance_contact_at_ground_and_C_height_R():
    # pts = (stance_contact, stance_C, knee_st, hip, knee_sw, swing_C, swing_lowest)
    pts = rocker_kneed_joints(
        np.array([0.2, -0.2, -0.2, 0, 0, 0]), l_t=0.5, l_s=0.5, R=0.3, contact_x=0.0
    )
    assert np.isclose(pts[0][1], 0.0, atol=1e-12)  # 接触点 y=0
    assert np.isclose(pts[1][1], 0.3, atol=1e-12)  # 曲率中心 C 高さ R
