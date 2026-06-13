# src/crane/stride.py
"""stride 写像: 衝突直後状態から各相を順に積分し、最終 heel-strike 後の状態を返す。

各相の終端面 g(x)=0 は複数回交差しうる（foot scuffing 等）。本物の相終端の
判定は PhaseSpec.event_accept に委譲する（Phase 1 の教訓:
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
    x_strike: np.ndarray  # 最終相の衝突写像適用前（heel-strike 瞬間の状態）
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
    """1歩分の写像（多相）。副作用なし。t_max は stride 全体の予算。"""
    t0 = 0.0
    x = np.asarray(x0, dtype=float)
    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []

    def append(t_seg, x_seg):
        if ts:  # 2セグメント目以降は境界点の重複を落とす
            t_seg, x_seg = t_seg[1:], x_seg[:, 1:]
        ts.append(t_seg)
        xs.append(x_seg)

    for i_phase, phase in enumerate(model.phases):

        def event(t, x_, _p=phase):
            return _p.event_value(x_)

        event.terminal = True
        event.direction = 0

        while True:
            if t0 >= t_max:
                raise StrideError(
                    f"phase {i_phase}: no terminal event before t_max={t_max} (t={t0:.3f})"
                )
            burn = solve_ivp(phase.dynamics, (t0, t0 + T_BURN), x, rtol=rtol, atol=atol)
            append(burn.t, burn.y)
            t0, x = burn.t[-1], burn.y[:, -1]
            if t0 >= t_max:
                continue  # 次周の冒頭ガードで StrideError

            sol = solve_ivp(phase.dynamics, (t0, t_max), x, events=event, rtol=rtol, atol=atol)
            append(sol.t, sol.y)
            if sol.t_events[0].size == 0:
                raise StrideError(f"phase {i_phase}: no terminal event before t_max={t_max}")
            t_e = sol.t_events[0][0]
            x_e = sol.y_events[0][0]
            if phase.event_accept(x_e):
                break
            t0, x = t_e, x_e  # 偽交差（scuffing 等）: 無視して続行

        x_strike, t0, x = x_e, t_e, phase.impact(x_e)

    return StrideResult(
        x_end=x,
        x_strike=x_strike,  # 最終相の衝突直前状態（heel-strike 瞬間）
        t_step=t0,
        t=np.concatenate(ts),
        x=np.concatenate(xs, axis=1),
    )
