import numpy as np

from crane.models.simplest import SimplestParams, dynamics, heelstrike_map, lift


P = SimplestParams(gamma=0.009)


def test_equilibrium_state_has_zero_acceleration():
    """θ=γ, φ=0, 静止状態では加速度ゼロ（直立平衡）。"""
    x = [P.gamma, 0.0, 0.0, 0.0]
    xdot = dynamics(0.0, x, P)
    assert np.allclose(xdot, [0.0, 0.0, 0.0, 0.0], atol=1e-15)


def test_heelstrike_geometry_is_exact():
    """θ⁺=−θ⁻, φ⁺=−2θ⁻ は厳密（脚交換の鏡映）。"""
    x_minus = np.array([-0.2, -0.4, -0.25, -0.01])
    x_plus = heelstrike_map(x_minus)
    assert x_plus[0] == 0.2
    assert x_plus[1] == 0.4


def test_heelstrike_velocity_map():
    """θ̇⁺ = cos(2θ⁻)θ̇⁻, φ̇⁺ = cos(2θ⁻)(1−cos(2θ⁻))θ̇⁻。"""
    theta, theta_dot = -0.2, -0.25
    x_plus = heelstrike_map(np.array([theta, 2 * theta, theta_dot, -0.01]))
    c = np.cos(2 * theta)
    assert np.isclose(x_plus[2], c * theta_dot)
    assert np.isclose(x_plus[3], c * (1 - c) * theta_dot)


def test_post_impact_state_is_on_section():
    """衝突直後は φ=2θ かつ φ̇=(1−cos2θ)θ̇ が成立 = lift と整合。"""
    x_plus = heelstrike_map(np.array([-0.2, -0.4, -0.25, -0.01]))
    expected = lift(x_plus[0], x_plus[2])
    assert np.allclose(x_plus, expected)
