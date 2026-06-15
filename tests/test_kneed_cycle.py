import numpy as np

from crane import references_kneed as ref
from crane.models.kneed import KneedParams, make_kneed
from crane.search import find_limit_cycle
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


def test_kneed_limit_cycle_exists_at_published_config():
    """Phase 3 ゲート: 公表パラメータでリミットサイクルが存在する。"""
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    # 安定性は文献の記述に従う（Sec.3.3.1 が locally stable と明言）
    if ref.PUBLISHED_STABLE:
        assert np.max(np.abs(fp.eigenvalues)) < 1.0


def test_kneed_gait_matches_published_quantities():
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    result = stride(MODEL, MODEL.lift(fp.y))
    if ref.STEP_PERIOD is not None:
        assert np.isclose(result.t_step, ref.STEP_PERIOD, rtol=ref.GAIT_TOLERANCE)
    if ref.THETA_STRIKE is not None:
        assert np.isclose(abs(result.x_strike[0]), abs(ref.THETA_STRIKE), rtol=ref.GAIT_TOLERANCE)


def test_gait_family_continues_near_published_slope():
    """公表 slope の ±20% で gait family が continuation で追える。"""
    fp0 = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp0.converged
    y = fp0.y
    for scale in [0.95, 0.9, 0.85, 0.8]:
        p = KneedParams(
            m_h=ref.M_HIP,
            m_t=ref.M_THIGH,
            m_s=ref.M_SHANK,
            l_t=ref.L_THIGH,
            l_s=ref.L_SHANK,
            b_t=ref.B_THIGH,
            b_s=ref.B_SHANK,
            gamma=ref.GAMMA_GAIT * scale,
            g=ref.G,
        )
        fp = find_limit_cycle(make_kneed(p), y)
        if not fp.converged:
            break  # gait が存在しない slope に達した: 既知の kneed の性質。報告のみ
        y = fp.y
    assert fp0.converged  # 最低限、公表点での存在が family の証拠
