import numpy as np

from crane import references_kneed as ref
from crane.models.kneed import KneedParams, make_kneed
from crane.stride import stride

P = KneedParams(
    m_h=ref.M_HIP,
    m_t=ref.M_THIGH,
    m_s=ref.M_SHANK,
    l_t=ref.L_THIGH,
    l_s=ref.L_SHANK,
    b_t=ref.B_THIGH,
    b_s=ref.B_SHANK,
    gamma=ref.GAMMA_GAIT,
    g=ref.G,
)
MODEL = make_kneed(P)


def test_stride_passes_both_phases_and_returns_to_section():
    """文献 seed から 1 stride: knee-strike → heel-strike を経て断面に戻る。"""
    result = stride(MODEL, MODEL.lift(np.array(ref.SECTION_GUESS)))
    x_end = result.x_end
    assert np.isclose(x_end[1], -x_end[0], atol=1e-8)  # θ_th = −θ_st
    assert np.isclose(x_end[2], -x_end[0], atol=1e-8)  # θ_sh = −θ_st
    assert result.t_step > 0.1
    # trajectory 中に膝の屈曲が現れている（unlocked 相の実在確認）
    flexion = result.x[2] - result.x[1]  # θ_sh − θ_th
    assert np.max(np.abs(flexion)) > 1e-3


def test_stride_reports_knee_flexion_sign():
    """屈曲方向が references_kneed の KNEE_FLEXION_SIGN と一致（規約照合）。"""
    result = stride(MODEL, MODEL.lift(np.array(ref.SECTION_GUESS)))
    flexion = result.x[2] - result.x[1]
    peak = flexion[np.argmax(np.abs(flexion))]
    assert np.sign(peak) == ref.KNEE_FLEXION_SIGN
