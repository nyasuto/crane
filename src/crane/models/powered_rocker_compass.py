"""動力付き rocker_compass。rocker_compass に pre-emptive 撃力 push-off を足す。

後脚（trailing = old stance）軸方向の撃力 Î_po = push_off·(P_st-hip)/|hip-P_st| を
heelstrike 衝突の直前に注入。撃力は接触点 P_st と hip を結ぶ脚線上にあり hip を通るため、
hip まわりのモーメントはゼロ。よって角運動量保存の eq1（系の新接触点まわり）にのみ
push-off モーメント M_po = cross2d(P_st - foot_sw_contact, Î_po) を加え、eq2（swing 脚の
hip まわり）は不変。push_off=0 で rocker_compass の衝突写像に厳密一致。

push-off 公式は実装前 de-risk で検証済み（push_off→0 で diff 0.0、push_off>0 でエネルギー注入）。
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum, cross2d
from crane.derive.lagrange import derive_qdd
from crane.model import HybridModel, PhaseSpec


@dataclass(frozen=True)
class PoweredRockerCompassParams:
    m: float
    m_h: float
    c: float
    rho: float
    L: float
    R: float
    gamma: float
    push_off: float  # pre-emptive 撃力 push-off の大きさ（0 で受動）
    g: float = 9.81


@lru_cache(maxsize=1)
def _build():
    th_st, th_sw, w_st, w_sw = sp.symbols("th_st th_sw w_st w_sw")
    m, m_h, c, rho, L, R, gamma, g, po = sp.symbols("m m_h c rho L R gamma g po", positive=True)
    q, qd = [th_st, th_sw], [w_st, w_sw]

    def down(theta):
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    C_st = sp.Matrix([-R * th_st, R])
    hip = C_st - (L - R) * down(th_st)
    p_st1 = hip + (c - rho) * down(th_st)
    p_st2 = hip + (c + rho) * down(th_st)
    p_sw1 = hip + (c - rho) * down(th_sw)
    p_sw2 = hip + (c + rho) * down(th_sw)
    C_sw = hip + (L - R) * down(th_sw)
    foot_sw_contact = sp.Matrix([C_sw[0], 0])
    P_st = sp.Matrix([-R * th_st, 0])

    bodies = [
        (m_h, hip),
        (m / 2, p_st1),
        (m / 2, p_st2),
        (m / 2, p_sw1),
        (m / 2, p_sw2),
    ]

    def vel(r):
        return r.jacobian(q) * sp.Matrix(qd)

    T = sum(mass * vel(r).dot(vel(r)) for mass, r in bodies) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])
    V = -g * sum(mass * r.dot(g_dir) for mass, r in bodies)

    qdd = derive_qdd(q, qd, T, V, simplify=False)
    params = (m, m_h, c, rho, L, R, gamma, g, po)
    args = (th_st, th_sw, w_st, w_sw, *params)
    f_qdd = sp.lambdify(args, [qdd[0], qdd[1]], "numpy")
    f_kinetic = sp.lambdify(args, T, "numpy")

    st_masses = [(m / 2, p_st1), (m / 2, p_st2)]
    sw_masses = [(m / 2, p_sw1), (m / 2, p_sw2)]
    L_sys_stance = angular_momentum(bodies, q, qd, P_st)
    L_sys_swfoot = angular_momentum(bodies, q, qd, foot_sw_contact)
    L_st_hip = angular_momentum(st_masses, q, qd, hip)
    L_sw_hip = angular_momentum(sw_masses, q, qd, hip)

    # pre-emptive push-off: 後脚脚線方向の撃力（hip を通るため eq2 は不変、eq1 のみ）
    leg = hip - P_st
    legn = sp.sqrt(leg.dot(leg))
    impulse = po * (P_st - hip) / legn  # エネルギー注入符号（de-risk 検証済み）
    M_po = cross2d(P_st - foot_sw_contact, impulse)

    wp_st, wp_sw = sp.symbols("wp_st wp_sw")
    swap = {th_st: th_sw, th_sw: th_st, w_st: wp_st, w_sw: wp_sw}
    eqs = [
        L_sys_stance.subs(swap, simultaneous=True) - L_sys_swfoot - M_po,
        L_sw_hip.subs(swap, simultaneous=True) - L_st_hip,
    ]
    A, rhs = sp.linear_eq_to_matrix(eqs, [wp_st, wp_sw])
    qd_post = A.LUsolve(rhs)
    f_impact = sp.lambdify(args, [qd_post[0], qd_post[1]], "numpy")

    return f_qdd, f_kinetic, f_impact


def _args(x, p: PoweredRockerCompassParams):
    return (*x, p.m, p.m_h, p.c, p.rho, p.L, p.R, p.gamma, p.g, p.push_off)


def dynamics(t: float, x, p: PoweredRockerCompassParams):
    f_qdd = _build()[0]
    qdd = f_qdd(*_args(x, p))
    return [x[2], x[3], qdd[0], qdd[1]]


def powered_heelstrike_map(x: np.ndarray, p: PoweredRockerCompassParams) -> np.ndarray:
    """衝突写像（push-off モーメント付き）＋脚交換。push_off=0 で受動に一致。"""
    wp = _build()[2](*_args(x, p))
    return np.array([x[1], x[0], float(wp[0]), float(wp[1])])


def make_powered_rocker_compass(p: PoweredRockerCompassParams) -> HybridModel:
    return HybridModel(
        phases=(
            PhaseSpec(
                dynamics=lambda t, x: dynamics(t, x, p),
                event_value=lambda x: x[0] + x[1],
                event_accept=lambda x: x[0] < 0.0 and (x[2] + x[3]) < 0.0,
                impact=lambda x: powered_heelstrike_map(x, p),
            ),
        ),
        lift=lambda y: np.array([y[0], -y[0], y[1], y[2]]),
        project=lambda x: np.array([x[0], x[2], x[3]]),
    )
