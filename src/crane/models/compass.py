# src/crane/models/compass.py
"""Goswami 点足 compass gait。derive レイヤーで記号導出し lambdify。

座標: x = [θ_st, θ_sw, θ̇_st, θ̇_sw]（slope 法線基準の絶対角、進行 +x）。
stance 足原点の pinned パラメータ化。strike 面: g(x) = θ_st + θ_sw = 0。
受理: θ_st < 0 かつ swing 足下降（ẏ_sw = −l sinθ_st (θ̇_st+θ̇_sw) < 0、
θ_st<0 では θ̇_st+θ̇_sw < 0 と等価）。

衝突（plastic）: pre は pre-pinned 系で、post は post-pinned 系（脚ラベル
交換済みの同形の式）で評価した角運動量を等置する:
  1. 系全体: post の pivot 回り = pre の swing 足先回り（同一物理点 = 新接地点）
  2. trailing 脚: post の swing 脚 hip 回り = pre の stance 脚 hip 回り（同一物理脚）
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum
from crane.derive.lagrange import derive_qdd
from crane.model import HybridModel, PhaseSpec


@dataclass(frozen=True)
class CompassParams:
    m: float  # 脚質量 [kg]
    m_h: float  # hip 質量 [kg]
    a: float  # 足 → 脚質量点 距離 [m]
    b: float  # 脚質量点 → hip 距離 [m]
    gamma: float  # 斜面角 [rad]
    g: float = 9.81

    @property
    def l(self) -> float:  # noqa: E743
        return self.a + self.b


@lru_cache(maxsize=1)
def _build():
    """記号導出して lambdify した関数群を返す（モジュール内キャッシュ）。"""
    th_st, th_sw, w_st, w_sw = sp.symbols("th_st th_sw w_st w_sw")
    m, m_h, a, b, gamma, g = sp.symbols("m m_h a b gamma g", positive=True)
    length = a + b
    q, qd = [th_st, th_sw], [w_st, w_sw]

    def down(theta):
        """絶対角 theta の脚に沿って hip から下向きの単位ベクトル。"""
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    hip = sp.Matrix([-length * sp.sin(th_st), length * sp.cos(th_st)])
    p_st = hip + b * down(th_st)  # stance 脚質量点（足から a）
    p_sw = hip + b * down(th_sw)  # swing 脚質量点（hip から b）
    foot_sw = hip + length * down(th_sw)

    def vel(r):
        return r.jacobian(q) * sp.Matrix(qd)

    T = (
        m_h * vel(hip).dot(vel(hip)) + m * vel(p_st).dot(vel(p_st)) + m * vel(p_sw).dot(vel(p_sw))
    ) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])  # slope frame の重力方向
    V = -g * (m_h * hip.dot(g_dir) + m * p_st.dot(g_dir) + m * p_sw.dot(g_dir))

    qdd = derive_qdd(q, qd, T, V, simplify=False)
    args = (th_st, th_sw, w_st, w_sw, m, m_h, a, b, gamma, g)
    f_qdd = sp.lambdify(args, [qdd[0], qdd[1]], "numpy")
    f_energy = sp.lambdify(args, T + V, "numpy")
    f_kinetic = sp.lambdify(args, T, "numpy")

    # --- 衝突写像（docstring の定式化） ---
    bodies = [(m_h, hip), (m, p_st), (m, p_sw)]
    origin = sp.Matrix([0, 0])
    L_sys_pivot = angular_momentum(bodies, q, qd, origin)
    L_sys_swfoot = angular_momentum(bodies, q, qd, foot_sw)
    L_st_hip = angular_momentum([(m, p_st)], q, qd, hip)
    L_sw_hip = angular_momentum([(m, p_sw)], q, qd, hip)

    wp_st, wp_sw = sp.symbols("wp_st wp_sw")  # post 速度（post ラベル）
    swap = {th_st: th_sw, th_sw: th_st, w_st: wp_st, w_sw: wp_sw}
    eqs = [
        L_sys_pivot.subs(swap, simultaneous=True) - L_sys_swfoot,
        L_sw_hip.subs(swap, simultaneous=True) - L_st_hip,
    ]
    A, rhs = sp.linear_eq_to_matrix(eqs, [wp_st, wp_sw])
    qd_post = A.LUsolve(rhs)
    f_impact = sp.lambdify(args, [qd_post[0], qd_post[1]], "numpy")

    return f_qdd, f_energy, f_kinetic, f_impact


def _args(x, p: CompassParams):
    return (*x, p.m, p.m_h, p.a, p.b, p.gamma, p.g)


def dynamics(t: float, x, p: CompassParams):
    f_qdd, _, _, _ = _build()
    qdd = f_qdd(*_args(x, p))
    return [x[2], x[3], qdd[0], qdd[1]]


def energy(x, p: CompassParams) -> float:
    _, f_energy, _, _ = _build()
    return float(f_energy(*_args(x, p)))


def kinetic_energy(x, p: CompassParams) -> float:
    _, _, f_kinetic, _ = _build()
    return float(f_kinetic(*_args(x, p)))


def heelstrike_map(x: np.ndarray, p: CompassParams) -> np.ndarray:
    """衝突写像 + 脚ラベル交換。post 速度は post ラベルで返る。"""
    _, _, _, f_impact = _build()
    wp = f_impact(*_args(x, p))
    return np.array([x[1], x[0], float(wp[0]), float(wp[1])])


def make_compass(p: CompassParams) -> HybridModel:
    return HybridModel(
        phases=(
            PhaseSpec(
                dynamics=lambda t, x: dynamics(t, x, p),
                event_value=lambda x: x[0] + x[1],
                event_accept=lambda x: x[0] < 0.0 and (x[2] + x[3]) < 0.0,
                impact=lambda x: heelstrike_map(x, p),
            ),
        ),
        lift=lambda y: np.array([y[0], -y[0], y[1], y[2]]),
        project=lambda x: np.array([x[0], x[2], x[3]]),
    )
