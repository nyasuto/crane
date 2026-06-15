# src/crane/models/kneed.py
"""点足 kneed walker。derive レイヤーで記号導出し lambdify。

状態 x = [θ_st, θ_th, θ_sh, θ̇_st, θ̇_th, θ̇_sh]（slope 法線基準の絶対角）。
θ_st: stance 脚（膝ロック・直線）、θ_th: swing 大腿、θ_sh: swing 脛。
locked 相は θ_th=θ_sh を 6D 埋め込み（q̈_th=q̈_sh）で維持する。

相機械: unlocked → knee-strike → locked → heel-strike（早見表参照）。
簡略化（docstring 明記必須）:
- knee-strike イベント g=θ_sh−θ_th=0 の transversal 交差は常に受理
  （ハイパーエクステンション抑止はロック係合そのもの）
- unlocked 相中の swing 足の地面通過（scuff）はイベントにしない（文献慣行）

heel-strike 定式化（Hsu Chen 2007 / McGeer 機械式ロックと整合）:
- 衝突の瞬間は両脚とも剛体（膝ロック維持）。撃力通過後に新 swing 膝がアンロック。
- 保存則: (i) 系全体を新接地点回りで保存、(ii) 新 swing 脚（旧 stance 脚、剛体のまま）
  を hip 回りで保存。post は θ̇_th⁺ = θ̇_sh⁺（衝突中は剛体）。
- 【当初案との差異】Hsu Chen (2007) Sec. 3.2.2 は q̇3⁺ = q̇2⁺ を運動学的仮定として
  採用（保存則 2 本のみ）。当初案の「新 swing 脛 knee 回り保存」3 本目は不採用。
  differences are small (m_s/m_t = 0.1) but systematic; recorded in references_kneed.py.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum
from crane.derive.lagrange import derive_qdd
from crane.model import HybridModel, PhaseSpec


@dataclass(frozen=True)
class KneedParams:
    m_h: float  # hip 質量 [kg]
    m_t: float  # 大腿質量 [kg]
    m_s: float  # 脛質量 [kg]
    l_t: float  # 大腿長 hip→knee [m]
    l_s: float  # 脛長 knee→foot [m]
    b_t: float  # hip→大腿質量点 [m]
    b_s: float  # knee→脛質量点 [m]
    gamma: float  # 斜面角 [rad]
    g: float = 9.81

    @property
    def l(self) -> float:  # noqa: E743
        return self.l_t + self.l_s


@lru_cache(maxsize=1)
def _build():
    """記号導出して lambdify した関数群を返す（モジュール内キャッシュ）。"""
    th_st, th_th, th_sh = sp.symbols("th_st th_th th_sh")
    w_st, w_th, w_sh = sp.symbols("w_st w_th w_sh")
    m_h, m_t, m_s = sp.symbols("m_h m_t m_s", positive=True)
    l_t, l_s, b_t, b_s = sp.symbols("l_t l_s b_t b_s", positive=True)
    gamma, g = sp.symbols("gamma g", positive=True)
    length = l_t + l_s

    def down(theta):
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    # stance 脚（直線、足原点）: hip は足から脚全長ぶん上
    hip = sp.Matrix([-length * sp.sin(th_st), length * sp.cos(th_st)])
    knee_st = hip + l_t * down(th_st)
    p_th_st = hip + b_t * down(th_st)  # stance 大腿質量点
    p_sh_st = knee_st + b_s * down(th_st)  # stance 脛質量点
    # swing 脚
    p_th_sw = hip + b_t * down(th_th)
    knee_sw = hip + l_t * down(th_th)
    p_sh_sw = knee_sw + b_s * down(th_sh)
    foot_sw = knee_sw + l_s * down(th_sh)

    q3, qd3 = [th_st, th_th, th_sh], [w_st, w_th, w_sh]
    bodies = [(m_h, hip), (m_t, p_th_st), (m_s, p_sh_st), (m_t, p_th_sw), (m_s, p_sh_sw)]

    def vel(r, q, qd):
        return r.jacobian(q) * sp.Matrix(qd)

    T3 = sum(mass * vel(r, q3, qd3).dot(vel(r, q3, qd3)) for mass, r in bodies) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])
    V = -g * sum(mass * r.dot(g_dir) for mass, r in bodies)

    params = (m_h, m_t, m_s, l_t, l_s, b_t, b_s, gamma, g)
    args6 = (th_st, th_th, th_sh, w_st, w_th, w_sh, *params)

    # unlocked: 3 DOF
    qdd3 = derive_qdd(q3, qd3, T3, V, simplify=False)
    f_qdd3 = sp.lambdify(args6, [qdd3[0], qdd3[1], qdd3[2]], "numpy")
    f_energy = sp.lambdify(args6, T3 + V, "numpy")
    f_kinetic = sp.lambdify(args6, T3, "numpy")

    # locked: θ_th=θ_sh=θ_sw の 2 DOF を導出して 6D に埋め込む
    th_sw, w_sw = sp.symbols("th_sw w_sw")
    lock = {th_th: th_sw, th_sh: th_sw, w_th: w_sw, w_sh: w_sw}
    T2 = T3.subs(lock, simultaneous=True)
    V2 = V.subs(lock, simultaneous=True)
    qdd2 = derive_qdd([th_st, th_sw], [w_st, w_sw], T2, V2, simplify=False)
    f_qdd2 = sp.lambdify((th_st, th_sw, w_st, w_sw, *params), [qdd2[0], qdd2[1]], "numpy")

    # --- knee-strike 衝突（脚交換なし、3 速度 → 2 速度）---
    L_sys_pivot3 = angular_momentum(bodies, q3, qd3, sp.Matrix([0, 0]))
    L_swleg_hip3 = angular_momentum([(m_t, p_th_sw), (m_s, p_sh_sw)], q3, qd3, hip)
    wp_st, wp_sw = sp.symbols("wp_st wp_sw")
    post_ks = {w_st: wp_st, w_th: wp_sw, w_sh: wp_sw}  # 位置不変（θ_th=θ_sh 上）
    eqs_ks = [
        L_sys_pivot3.subs(post_ks, simultaneous=True) - L_sys_pivot3,
        L_swleg_hip3.subs(post_ks, simultaneous=True) - L_swleg_hip3,
    ]
    A_ks, rhs_ks = sp.linear_eq_to_matrix(eqs_ks, [wp_st, wp_sw])
    qd_post_ks = A_ks.LUsolve(rhs_ks)
    f_impact_ks = sp.lambdify(args6, [qd_post_ks[0], qd_post_ks[1]], "numpy")

    # --- heel-strike 衝突（脚交換、衝突中は両脚剛体 = 膝ロック維持。2 速度 → 2 速度）---
    # 【Task 2 文献調査に基づく定式化: Hsu Chen 2007 / McGeer の機械式ロックと整合】
    # pre は locked 配置（θ_th=θ_sh）の pre-pinned 系で評価。
    # post は post-pinned 系（同形の式）に swap 代入:
    #   位置: th_st→th_th（新 stance = 旧 swing 角）、th_th→th_st, th_sh→th_st
    #   速度: w_st→wq_st, w_th→wq_leg, w_sh→wq_leg（新 swing 脚は衝突中剛体）
    wq_st, wq_leg = sp.symbols("wq_st wq_leg")
    swap_hs = {
        th_st: th_th,
        th_th: th_st,
        th_sh: th_st,
        w_st: wq_st,
        w_th: wq_leg,
        w_sh: wq_leg,
    }
    L_stleg_hip3 = angular_momentum([(m_t, p_th_st), (m_s, p_sh_st)], q3, qd3, hip)
    L_sys_swfoot3 = angular_momentum(bodies, q3, qd3, foot_sw)
    eqs_hs = [
        L_sys_pivot3.subs(swap_hs, simultaneous=True) - L_sys_swfoot3,
        L_swleg_hip3.subs(swap_hs, simultaneous=True) - L_stleg_hip3,
    ]
    A_hs, rhs_hs = sp.linear_eq_to_matrix(eqs_hs, [wq_st, wq_leg])
    qd_post_hs = A_hs.LUsolve(rhs_hs)
    f_impact_hs = sp.lambdify(args6, [qd_post_hs[0], qd_post_hs[1]], "numpy")

    return f_qdd3, f_qdd2, f_energy, f_kinetic, f_impact_ks, f_impact_hs


def _args(x, p: KneedParams):
    return (*x, p.m_h, p.m_t, p.m_s, p.l_t, p.l_s, p.b_t, p.b_s, p.gamma, p.g)


def dynamics_unlocked(t: float, x, p: KneedParams):
    """unlocked 相（3 DOF）の時間微分を返す。"""
    f_qdd3 = _build()[0]
    qdd = f_qdd3(*_args(x, p))
    return [x[3], x[4], x[5], qdd[0], qdd[1], qdd[2]]


def dynamics_locked(t: float, x, p: KneedParams):
    """locked 相: 2 DOF を 6D に埋め込み（θ_th スロットを θ_sw として使う）。"""
    f_qdd2 = _build()[1]
    qdd = f_qdd2(
        x[0], x[1], x[3], x[4], p.m_h, p.m_t, p.m_s, p.l_t, p.l_s, p.b_t, p.b_s, p.gamma, p.g
    )
    return [x[3], x[4], x[5], qdd[0], qdd[1], qdd[1]]


def energy(x, p: KneedParams) -> float:
    """全力学的エネルギー（運動 + ポテンシャル）。"""
    return float(_build()[2](*_args(x, p)))


def kinetic_energy(x, p: KneedParams) -> float:
    """運動エネルギー。"""
    return float(_build()[3](*_args(x, p)))


def kneestrike_map(x: np.ndarray, p: KneedParams) -> np.ndarray:
    """knee-strike: 脚交換なし。swing 膝がロックし 3 速度 → 2 速度。"""
    wp = _build()[4](*_args(x, p))
    return np.array([x[0], x[1], x[2], float(wp[0]), float(wp[1]), float(wp[1])])


def heelstrike_map(x: np.ndarray, p: KneedParams) -> np.ndarray:
    """heel-strike: 脚交換。衝突中は両脚剛体（膝ロック維持、Hsu Chen 2007 /
    McGeer の機械式ロックと整合）。アンロックは次相開始時 = post で
    θ̇_th⁺ = θ̇_sh⁺。post 速度は post ラベル。"""
    wq = _build()[5](*_args(x, p))
    return np.array([x[1], x[0], x[0], float(wq[0]), float(wq[1]), float(wq[1])])


def make_kneed(p: KneedParams) -> HybridModel:
    """相機械: unlocked → knee-strike → locked → heel-strike。

    断面 = heel-strike 直後。y = (θ_st, θ̇_st, θ̇_sw) の 3 次元
    （位置 θ_th = θ_sh = −θ_st、速度 θ̇_th = θ̇_sh = θ̇_sw が成立）。
    heel-strike 受理は compass と同じ閉形式（脚全長 l で
    ẏ_sw = −l sin(θ_st)(θ̇_st + θ̇_sw) < 0、locked 相では θ̇_sw = x[4]）。
    """
    unlocked = PhaseSpec(
        dynamics=lambda t, x: dynamics_unlocked(t, x, p),
        event_value=lambda x: x[2] - x[1],  # θ_sh − θ_th = 0 で knee-strike
        event_accept=lambda x: True,  # transversal 交差は常にロック係合
        impact=lambda x: kneestrike_map(x, p),
    )
    locked = PhaseSpec(
        dynamics=lambda t, x: dynamics_locked(t, x, p),
        event_value=lambda x: x[0] + x[1],  # θ_st + θ_sw = 0 で heel-strike
        event_accept=lambda x: x[0] < 0.0 and (x[3] + x[4]) < 0.0,
        impact=lambda x: heelstrike_map(x, p),
    )
    return HybridModel(
        phases=(unlocked, locked),
        lift=lambda y: np.array([y[0], -y[0], -y[0], y[1], y[2], y[2]]),
        project=lambda x: np.array([x[0], x[3], x[4]]),
    )
