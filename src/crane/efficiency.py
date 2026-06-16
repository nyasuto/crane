"""受動歩行の効率指標（model 非依存）。

主指標は一歩あたり相対損失 δ = (KE_pre - KE_post)/KE_pre（車輪極限で →0）。
機械的 COT はリミットサイクル上で勾配に固定される（≈ sin γ）ので副次の内部チェック用。
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def relative_loss(ke_pre: float, ke_post: float) -> float:
    """一歩あたり相対損失 δ = (KE_pre - KE_post)/KE_pre。車輪極限で →0。"""
    return (ke_pre - ke_post) / ke_pre


def step_collision_loss(
    x_strike: np.ndarray,
    x_end: np.ndarray,
    kinetic_energy: Callable[[np.ndarray], float],
) -> tuple[float, float, float]:
    """heel-strike の衝突損失。(loss, ke_pre, ke_post) を返す。

    x_strike: 衝突直前の全状態、x_end: 衝突写像適用後（次断面）の全状態。
    KE は脚ラベル交換で不変なので ke_post は衝突直後 KE に等しい。
    """
    ke_pre = float(kinetic_energy(x_strike))
    ke_post = float(kinetic_energy(x_end))
    return ke_pre - ke_post, ke_pre, ke_post


def mechanical_cot(loss: float, *, m: float, g: float, step_length: float) -> float:
    """機械的 cost of transport = 衝突損失 / (m·g·一歩進行距離)。

    step_length は斜面（接触フレーム）方向の進行距離（鉛直落下 = step_length·sin γ）。
    リミットサイクル上ではエネルギー収支より ≈ sin γ になるはず（内部チェック）。
    """
    return loss / (m * g * step_length)
