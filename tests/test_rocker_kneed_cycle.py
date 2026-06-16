# tests/test_rocker_kneed_cycle.py
"""Phase 3.6 文献ゲート: McGeer 1990b round-foot kneed の公表 config でリミットサイクルを
発見し、論文値と照合する。一次ゲートは R→0 退化（test_rocker_kneed_pointfoot_limit.py）。
本ファイルは AUXILIARY ゲート: McGeer は分布質量脚を使い本実装は点質量近似なので、固有値の
厳密一致は期待しない。よって hard-assert するのは「cycle 存在 + 安定 + step period（許容内）」、
固有値は dominant |λ| の参考照合（次元差 3 vs 4 でクラッシュさせない）。
"""

import numpy as np

from crane import references_mcgeer_knees as ref
from crane.models.rocker_kneed import RockerKneedParams, make_rocker_kneed
from crane.search import find_limit_cycle
from crane.stride import stride

P = RockerKneedParams(
    m_h=ref.M_HIP,
    m_t=ref.M_THIGH,
    m_s=ref.M_SHANK,
    l_t=ref.L_THIGH,
    l_s=ref.L_SHANK,
    b_t=ref.B_THIGH,
    b_s=ref.B_SHANK,
    R=ref.R_FOOT,
    gamma=ref.GAMMA_GAIT,
    g=ref.G,
)
MODEL = make_rocker_kneed(P)


def test_rocker_kneed_limit_cycle_exists_and_stable():
    """McGeer config でサイクルが存在し、安定（max|λ|<1）。一次主張: 存在＋安定。"""
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    assert fp.eigenvalues is not None
    if ref.PUBLISHED_STABLE:
        assert np.max(np.abs(fp.eigenvalues)) < 1.0


def test_rocker_kneed_step_period_matches_published():
    """step period を秒で論文値と照合。STEP_PERIOD は √(l/g) 単位なので秒に変換する（罠 1）。"""
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    result = stride(MODEL, MODEL.lift(fp.y))
    if ref.STEP_PERIOD is not None:
        # √(l/g) → 秒（本モデルは G=9.81, l=l_t+l_s で次元を持つ）
        expected = ref.STEP_PERIOD * np.sqrt((ref.L_THIGH + ref.L_SHANK) / ref.G)
        assert np.isclose(result.t_step, expected, rtol=ref.GAIT_TOLERANCE)


def test_rocker_kneed_dominant_eigenvalue_reference():
    """固有値は参考比較のみ（点質量 vs 分布質量 + 断面次元差 3 vs 4）。

    本実装の first-return map は 3D なので fp.eigenvalues は 3 個、McGeer の
    POINCARE_EIGENVALUES は 4 個。配列の直接 allclose は shape 不一致でクラッシュするので
    行わない。dominant |λ| を報告し、安定だけを hard-assert する。
    """
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    ours_max = float(np.max(np.abs(fp.eigenvalues)))
    # 安定は必須（dominant |λ| < 1）
    assert ours_max < 1.0
    if ref.POINCARE_EIGENVALUES is not None:
        theirs_max = float(np.max(np.abs(np.array(ref.POINCARE_EIGENVALUES))))
        # 参考: 両者とも安定マージン（|λ|<1）の同オーダー。点質量近似のため厳密一致は非期待。
        # 同一性は assert しない（分布質量差 + 3 vs 4 次元差）。報告のみ。
        print(f"dominant |lambda|: ours={ours_max:.4f} mcgeer={theirs_max:.4f}")
        assert theirs_max < 1.0  # McGeer も安定（Table 1 全固有値 |z|<1）
