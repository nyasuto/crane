import numpy as np
import pytest

from crane import references as ref
from crane.models.simplest import SimplestParams, lift
from crane.stride import StrideError, stride


P = SimplestParams(gamma=ref.GAMMA_REF)


def test_stride_from_reference_guess_returns_to_section():
    """文献不動点近傍から1歩進むと、終状態が断面 (φ=2θ) 上にある。"""
    x0 = lift(ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT)
    result = stride(P, x0)
    assert result.t_step > 1.0
    assert np.isclose(result.x_end[1], 2.0 * result.x_end[0], atol=1e-8)


def test_stride_records_trajectory():
    """軌跡が時系列で記録される（観察可能性）。"""
    x0 = lift(ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT)
    result = stride(P, x0)
    assert result.t.shape[0] == result.x.shape[1]
    assert result.t[0] == 0.0
    assert np.all(np.diff(result.t) > 0)


def test_stride_raises_when_no_heelstrike():
    """heel-strike に至らない初期条件では StrideError。"""
    x0 = lift(0.001, 0.0)  # ほぼ直立静止 → 歩かない
    with pytest.raises(StrideError):
        stride(P, x0)


def test_stride_short_period_reaches_true_strike():
    """short-period 参照点からの1歩は τ₀=π 近傍で下降接地する。

    θ<0 だが足が上昇中の偽交差 (t≈2.13) を弾けることのゲート。
    許容幅は O(γ) 近似のずれ込み（references.py 参照）。
    """
    x0 = lift(ref.SHORT_PERIOD_THETA, ref.SHORT_PERIOD_THETA_DOT)
    result = stride(P, x0)
    assert np.isclose(result.t_step, np.pi, atol=0.1)
    assert result.x_strike[0] < -0.1  # 真の strike (θ≈-0.195)。偽交差なら -0.048
