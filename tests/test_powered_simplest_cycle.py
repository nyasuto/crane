import numpy as np

from crane.models.powered_simplest import PoweredSimplestParams, make_powered_simplest
from crane.search import find_limit_cycle
from crane.stride import stride


def _ke(theta_dot):
    # simplest walker の KE = (1/2)θ̇²（hip 点質量 M=1, 脚 massless, l=1）
    return 0.5 * theta_dot**2


def _find_level_cycle():
    """平地 γ=0 の動力サイクルを連続接続で探す。

    γ=0 では θ サブ系（θ̈=sinθ）が P・φ から独立なので、受動不動点近傍を
    seed にした単純な P スキャンは自明な静止不動点 (y→0) へ落ちる。
    そこで受動 γ=0.009 サイクルから出発し、γ を 0 まで下げながら push_off を
    増やして毎ステップ前解で seed する継続法を使う（plan Task 4 step 2）。
    非自明（|θ|>0.05）で収束した最終 γ=0 解を返す。
    """
    gammas = [0.009, 0.007, 0.005, 0.003, 0.0015, 0.0]
    push_offs = [0.0, 0.05, 0.07, 0.085, 0.10, 0.115]
    y = np.array([0.2003109, -0.1998325])  # 受動 γ=0.009 長周期不動点
    push_off = 0.0
    model = None
    for gamma, p in zip(gammas, push_offs):
        model = make_powered_simplest(PoweredSimplestParams(gamma=gamma, push_off=p))
        fp = find_limit_cycle(model, y)
        if fp.converged and abs(fp.y[0]) > 0.05:  # 自明な静止解を排除
            y, push_off = fp.y, p
    if gamma == 0.0 and fp.converged and abs(fp.y[0]) > 0.05:
        return push_off, model, fp
    raise AssertionError("平地動力サイクルが見つからない")


def test_level_ground_limit_cycle_exists_and_stable():
    push_off, model, fp = _find_level_cycle()
    assert fp.converged
    assert np.max(np.abs(fp.eigenvalues)) < 1.0
    assert fp.y[1] < 0.0  # 前進方向


def test_energy_balance_pushoff_work_equals_collision_loss():
    push_off, model, fp = _find_level_cycle()
    result = stride(model, model.lift(fp.y))
    td_strike = result.x_strike[2]
    td_end = result.x_end[2]
    ke_after_pushoff = 0.5 * (td_strike**2 + push_off**2)
    collision_loss = ke_after_pushoff - _ke(td_end)
    pushoff_work = 0.5 * push_off**2
    assert np.isclose(pushoff_work, collision_loss, atol=1e-3)
