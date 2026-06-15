# tests/test_rocker_kneed_stride.py
import numpy as np

from crane import references_kneed as kref
from crane.models.rocker_kneed import RockerKneedParams, make_rocker_kneed
from crane.stride import stride

P = RockerKneedParams(
    m_h=kref.M_HIP,
    m_t=kref.M_THIGH,
    m_s=kref.M_SHANK,
    l_t=kref.L_THIGH,
    l_s=kref.L_SHANK,
    b_t=kref.B_THIGH,
    b_s=kref.B_SHANK,
    R=0.2,
    gamma=kref.GAMMA_GAIT,
    g=kref.G,
)
MODEL = make_rocker_kneed(P)


def test_stride_passes_both_phases_and_returns_to_section():
    """seed から unlocked→knee-strike→locked→heel-strike を経て断面に戻る。"""
    result = stride(MODEL, MODEL.lift(np.array(kref.SECTION_GUESS)))
    x_end = result.x_end
    assert np.isclose(x_end[1], -x_end[0], atol=1e-8)  # θ_th = −θ_st
    assert np.isclose(x_end[2], -x_end[0], atol=1e-8)  # θ_sh = −θ_st
    assert result.t_step > 0.1
    flexion = result.x[2] - result.x[1]  # θ_sh − θ_th
    assert np.max(np.abs(flexion)) > 1e-3  # 膝屈曲が現れる
