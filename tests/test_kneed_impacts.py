import numpy as np

from crane import references_kneed as ref
from crane.models.kneed import (
    KneedParams,
    heelstrike_map,
    kinetic_energy,
    kneestrike_map,
)

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


def test_kneestrike_locks_and_dissipates():
    """knee-strike: 位置不変・速度等値化・KE 非増加。"""
    x_pre = np.array([-0.05, 0.25, 0.25, -1.3, -0.6, 1.8])  # θ_th=θ_sh（lock 面）
    x_post = kneestrike_map(x_pre, P)
    assert np.allclose(x_post[:3], x_pre[:3])  # 位置不変
    assert x_post[4] == x_post[5]  # 速度ロック
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_kneestrike_zero_velocity_is_fixed():
    x_pre = np.array([-0.05, 0.25, 0.25, 0.0, 0.0, 0.0])
    x_post = kneestrike_map(x_pre, P)
    assert np.allclose(x_post[3:], np.zeros(3), atol=1e-12)


def test_kneestrike_already_locked_is_identity():
    """既に速度が揃っている状態の knee-strike は恒等（撃力ゼロ）。"""
    x_pre = np.array([-0.05, 0.25, 0.25, -1.3, 0.7, 0.7])
    x_post = kneestrike_map(x_pre, P)
    assert np.allclose(x_post, x_pre, atol=1e-10)


def test_heelstrike_swaps_and_dissipates():
    """heel-strike: 脚交換幾何 + KE 非増加 + 新 swing 整列（剛体のまま）。"""
    th = -0.18
    x_pre = np.array([th, -th, -th, -1.4, -0.9, -0.9])  # locked, strike 面上
    x_post = heelstrike_map(x_pre, P)
    assert x_post[0] == -th  # 新 stance = 旧 swing 角
    assert x_post[1] == th  # 新 swing 大腿 = 旧 stance 角
    assert x_post[2] == th  # 新 swing 脛も整列
    assert x_post[4] == x_post[5]  # 衝突中は剛体: θ̇_th⁺ = θ̇_sh⁺
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_heelstrike_zero_velocity_is_fixed():
    th = -0.18
    x_pre = np.array([th, -th, -th, 0.0, 0.0, 0.0])
    x_post = heelstrike_map(x_pre, P)
    assert np.allclose(x_post[3:], np.zeros(3), atol=1e-12)
