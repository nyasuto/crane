"""Poincaré shooting と固有値解析のテスト。"""

import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle


MODEL = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))


def test_long_period_fixed_point_matches_garcia1998():
    """Phase 1 ゲート: 不動点が文献値と一致。

    許容誤差は絶対 1.5e-3: 参照値は論文印刷の O(γ) 展開係数からの導出値で、
    打ち切り誤差 O(γ^(5/3)) ≈ 4e-4〜8e-4 を持つ（references.py の docstring 参照）。
    真の数値不動点はこの範囲で参照値からずれるのが正しい挙動。
    符号・規約・転記ミスは ≥1e-2 のずれを生むので、このゲートで十分検出できる。
    """
    # ±5% は stride 領域外（即 StrideError）、±2% は basin 境界近くで J の固有値が +1 を跨ぎ Newton が発散。
    # ±1% は健全に収束（数値確認済み）
    guess = np.array([ref.LONG_PERIOD_THETA * 1.01, ref.LONG_PERIOD_THETA_DOT * 0.99])
    fp = find_limit_cycle(MODEL, guess)
    assert fp.converged
    assert np.isclose(fp.y[0], ref.LONG_PERIOD_THETA, rtol=0, atol=1.5e-3)
    assert np.isclose(fp.y[1], ref.LONG_PERIOD_THETA_DOT, rtol=0, atol=1.5e-3)


def test_short_period_fixed_point_matches_garcia1998():
    """short-period gait の不動点も文献値と一致（同じ許容誤差の根拠）。"""
    guess = np.array([ref.SHORT_PERIOD_THETA, ref.SHORT_PERIOD_THETA_DOT])
    fp = find_limit_cycle(MODEL, guess)
    assert fp.converged
    assert np.isclose(fp.y[0], ref.SHORT_PERIOD_THETA, rtol=0, atol=1.5e-3)
    assert np.isclose(fp.y[1], ref.SHORT_PERIOD_THETA_DOT, rtol=0, atol=1.5e-3)


def test_long_period_gait_is_stable():
    """γ=0.009 の long-period gait は漸近安定 (max|λ| < 1)。"""
    fp = find_limit_cycle(MODEL, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    assert fp.converged
    assert np.max(np.abs(fp.eigenvalues)) < 1.0
    if ref.LONG_PERIOD_EIGENVALUE_ABS is not None:
        ours = np.sort(np.abs(fp.eigenvalues))
        theirs = np.sort(np.array(ref.LONG_PERIOD_EIGENVALUE_ABS))
        assert np.allclose(ours, theirs, rtol=0.05)


def test_convergence_history_is_logged():
    """収束履歴が残る（観察可能性）。残差が初期より終端で小さい。"""
    # ±5% は stride 領域外（即 StrideError）、±2% は basin 境界近くで J の固有値が +1 を跨ぎ Newton が発散。
    # ±1% は健全に収束（数値確認済み）
    guess = np.array([ref.LONG_PERIOD_THETA * 1.01, ref.LONG_PERIOD_THETA_DOT * 0.99])
    fp = find_limit_cycle(MODEL, guess)
    residuals = [r for _, r in fp.history]
    assert len(residuals) >= 2
    assert residuals[-1] < residuals[0]


def test_backtracking_never_increases_residual():
    """Armijo: 残差が増加するステップは受理しない（履歴は単調非増加）。

    seed = fp*0.6 は basin 内だが、素の Newton では初手で大きくオーバーシュート
    して残差が 2.6e-2 → 1.18 へ膨張し 11 反復かかる（issue #1）。Armijo
    十分減少条件を入れると、各ステップで残差が非増加になり収束も速くなる。
    """
    fp_true = np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    fp = find_limit_cycle(MODEL, fp_true * 0.6)
    assert fp.converged
    residuals = [r for _, r in fp.history]
    assert len(residuals) >= 2
    for prev, nxt in zip(residuals, residuals[1:]):
        assert nxt <= prev, f"residual increased: {prev:.2e} -> {nxt:.2e}"
