import numpy as np

from crane.viz import kneed_joints


def test_kneed_joints_heelstrike_symmetry():
    """heel-strike 配置 (θ_th=θ_sh=−θ_st) で swing 足が接地高さにある。"""
    th = -0.18
    x = np.array([th, -th, -th, 0, 0, 0])
    pts = kneed_joints(x, l_t=0.5, l_s=0.5, foot_x=0.0)
    # pts = (stance_foot, stance_knee, hip, swing_knee, swing_foot)
    assert np.isclose(pts[4][1], 0.0, atol=1e-12)


def test_kneed_joints_straight_legs_match_two_segments():
    """両膝伸展時、knee は foot–hip 線分上にある。"""
    x = np.array([0.2, -0.3, -0.3, 0, 0, 0])
    pts = kneed_joints(x, l_t=0.5, l_s=0.5, foot_x=1.0)
    foot, knee_st, hip = pts[0], pts[1], pts[2]
    seg = hip - foot
    rel = knee_st - foot
    cross = seg[0] * rel[1] - seg[1] * rel[0]
    assert abs(cross) < 1e-12
