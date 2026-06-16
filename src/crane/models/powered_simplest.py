"""Kuo 2002 の動力付き simplest walker。

simplest walker（Garcia 1998）に pre-emptive 撃力 push-off を足し、平地 γ=0 を歩く。
連続相は simplest.dynamics を γ 経由で再利用し、heel-strike だけを
「後脚軸方向の撃力 P → 角運動量保存衝突 → 脚交換」の合成写像に差し替える。

push-off 写像（reduced 座標, M=l=g=1）: c=cos(2θ), s=sin(2θ),
    θ̇⁺ = c·θ̇ + s·P
P=0 で受動 heelstrike_map に厳密退化。push-off 仕事 = P²/2。
"""

from dataclasses import dataclass

import numpy as np

from crane.model import HybridModel, PhaseSpec
from crane.models.simplest import SimplestParams, dynamics, lift


@dataclass(frozen=True)
class PoweredSimplestParams:
    gamma: float  # slope angle [rad]（平地は 0）
    push_off: float  # pre-emptive 撃力 push-off の大きさ P（P=0 で受動）


def powered_heelstrike_map(x: np.ndarray, push_off: float) -> np.ndarray:
    """合成衝突写像: 後脚軸撃力 push_off → 角運動量保存衝突 → 脚交換。

    P=0 で simplest.heelstrike_map に厳密一致。
    """
    theta, _phi, theta_dot, _phi_dot = x
    c = np.cos(2.0 * theta)
    s = np.sin(2.0 * theta)
    td_plus = c * theta_dot + s * push_off
    return np.array([-theta, -2.0 * theta, td_plus, (1.0 - c) * td_plus])


def make_powered_simplest(p: PoweredSimplestParams) -> HybridModel:
    """パラメータ束縛済みの HybridModel。連続相は simplest と同一（γ 経由）。"""
    sp = SimplestParams(gamma=p.gamma)
    return HybridModel(
        phases=(
            PhaseSpec(
                dynamics=lambda t, x: dynamics(t, x, sp),
                event_value=lambda x: x[1] - 2.0 * x[0],
                event_accept=lambda x: x[0] < 0.0 and (x[3] - 2.0 * x[2]) > 0.0,
                impact=lambda x: powered_heelstrike_map(x, p.push_off),
            ),
        ),
        lift=lambda y: lift(y[0], y[1]),
        project=lambda x: np.array([x[0], x[2]]),
    )
