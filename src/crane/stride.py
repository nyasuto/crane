# src/crane/stride.py
"""stride 写像: 衝突直後状態から積分 → heel-strike → 衝突写像 → 次の衝突直後状態。

strike 面 g(x)=0 は1歩の間に複数回交差しうる（foot scuffing）。本物の
heel-strike の判定は model.strike_accept に委譲する（Phase 1 の教訓:
simplest walker では「θ<0 かつ swing 足が下降」が必要だった）。
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from crane.model import HybridModel

T_BURN = 1e-3  # event 無効の冒頭区間（衝突直後 g=0 と scuffing 直後の再発火対策）


class StrideError(RuntimeError):
    """heel-strike に到達しなかった（転倒・停止など）。"""


@dataclass
class StrideResult:
    x_end: np.ndarray  # 衝突写像適用後（次の断面上の状態）
    x_strike: np.ndarray  # 衝突写像適用前（heel-strike 瞬間の状態）
    t_step: float
    t: np.ndarray  # 軌跡時刻列（厳密単調増加）
    x: np.ndarray  # 軌跡状態列 shape (n_dof*2, N)


def stride(
    model: HybridModel,
    x0: np.ndarray,
    *,
    t_max: float = 10.0,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> StrideResult:
    """1歩分の写像。副作用なし。"""

    def event(t, x):
        return model.strike_value(x)

    event.terminal = True
    event.direction = 0

    t0 = 0.0
    x = np.asarray(x0, dtype=float)
    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []

    def append(t_seg, x_seg):
        if ts:  # 2セグメント目以降は境界点の重複を落とす
            t_seg, x_seg = t_seg[1:], x_seg[:, 1:]
        ts.append(t_seg)
        xs.append(x_seg)

    while t0 < t_max:
        burn = solve_ivp(model.dynamics, (t0, t0 + T_BURN), x, rtol=rtol, atol=atol)
        append(burn.t, burn.y)
        t0, x = burn.t[-1], burn.y[:, -1]
        if t0 >= t_max:
            break

        sol = solve_ivp(model.dynamics, (t0, t_max), x, events=event, rtol=rtol, atol=atol)
        append(sol.t, sol.y)

        if sol.t_events[0].size == 0:
            break

        t_e = sol.t_events[0][0]
        x_e = sol.y_events[0][0]
        if model.strike_accept(x_e):
            return StrideResult(
                x_end=model.impact(x_e),
                x_strike=x_e,
                t_step=t_e,
                t=np.concatenate(ts),
                x=np.concatenate(xs, axis=1),
            )
        t0, x = t_e, x_e  # scuffing: 無視して続行

    raise StrideError(f"no heelstrike before t_max={t_max}")
