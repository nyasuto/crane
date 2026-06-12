"""stride 写像: 衝突直後状態から積分 → heel-strike → 衝突写像 → 次の衝突直後状態。

heel-strike 判定の幾何学:
    面 g = φ−2θ = 0 は1歩の中で複数回交差することがある。
    swing 足高さ変化率 ẏ_sw = sin(θ)·ġ（ġ = φ̇−2θ̇）。
    θ<0 かつ ġ>0 のとき ẏ_sw<0 → swing 足が下降接地 → 本物の heel-strike。
    ġ≤0 の交差は scuff-dip 出口（足が上昇中）であり無視する。
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from crane.models.simplest import SimplestParams, dynamics, heelstrike_map

T_BURN = 1e-3  # event 無効の冒頭区間（衝突直後 g=0 と scuffing 直後の再発火対策）


class StrideError(RuntimeError):
    """heel-strike に到達しなかった（転倒・停止など）。"""


@dataclass
class StrideResult:
    x_end: np.ndarray  # 衝突写像適用後（次の断面上の状態）
    x_strike: np.ndarray  # 衝突写像適用前（heel-strike 瞬間の状態）
    t_step: float
    t: np.ndarray  # 軌跡時刻列
    x: np.ndarray  # 軌跡状態列 shape (4, N)


def _heelstrike_event(t, x, p):
    return x[1] - 2.0 * x[0]


_heelstrike_event.terminal = True
_heelstrike_event.direction = 0


def stride(
    p: SimplestParams,
    x0: np.ndarray,
    *,
    t_max: float = 10.0,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> StrideResult:
    """1歩分の写像。副作用なし。"""
    t0 = 0.0
    x = np.asarray(x0, dtype=float)
    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []

    def _append(t_seg: np.ndarray, x_seg: np.ndarray) -> None:
        if ts:  # 2セグメント目以降は境界点の重複を落とす
            t_seg, x_seg = t_seg[1:], x_seg[:, 1:]
        ts.append(t_seg)
        xs.append(x_seg)

    while t0 < t_max:
        # burn-in: event なしで微小区間を進める
        burn = solve_ivp(dynamics, (t0, t0 + T_BURN), x, args=(p,), rtol=rtol, atol=atol)
        _append(burn.t, burn.y)
        t0, x = burn.t[-1], burn.y[:, -1]

        if t0 >= t_max:
            break  # burn-in が t_max を踏み越えた場合はループを抜けて StrideError へ

        sol = solve_ivp(
            dynamics,
            (t0, t_max),
            x,
            args=(p,),
            events=_heelstrike_event,
            rtol=rtol,
            atol=atol,
        )
        _append(sol.t, sol.y)

        if sol.t_events[0].size == 0:
            break  # t_max まで heel-strike なし

        t_e = sol.t_events[0][0]
        x_e = sol.y_events[0][0]
        g_dot = x_e[3] - 2.0 * x_e[2]
        if (
            x_e[0] < 0.0 and g_dot > 0.0
        ):  # 本物の heel-strike: stance が鉛直越え かつ swing 足が下降接地
            return StrideResult(
                x_end=heelstrike_map(x_e),
                x_strike=x_e,
                t_step=t_e,
                t=np.concatenate(ts),
                x=np.concatenate(xs, axis=1),
            )
        t0, x = t_e, x_e  # scuffing: 無視して続行

    raise StrideError(f"no heelstrike before t_max={t_max}")
