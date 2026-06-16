# src/crane/basin.py
"""Basin of attraction の分類とスライス計算。

各 Poincaré 断面点から stride 写像を前進反復し「収束 / 転倒 / 未決」を分類する。
stride 写像は副作用なしの純関数なので、グリッド各点は独立に並列計算できる。
"""

from __future__ import annotations

import numpy as np

from crane.model import HybridModel
from crane.stride import StrideError, stride

CONVERGED = 0
FELL = 1
UNDECIDED = 2


def classify_ic(
    model: HybridModel,
    y0: np.ndarray,
    fixed_point: np.ndarray,
    *,
    max_strides: int = 20,
    converge_tol: float = 1e-3,
) -> int:
    """断面点 y0 から stride を最大 max_strides 反復し basin 分類を返す。

    deviation<converge_tol で CONVERGED / StrideError で FELL /
    上限まで未収束で UNDECIDED。

    収束判定は一発ヒューリスティック: 軌道が一時的に converge_tol 以内を
    通過しただけでも CONVERGED と判定する。局所安定な固定点に対する標準的
    かつ安全な basin ヒューリスティックである。
    """
    fp = np.asarray(fixed_point, dtype=float)
    x = model.lift(np.asarray(y0, dtype=float))
    for _ in range(max_strides):
        try:
            x = stride(model, x).x_end
        except StrideError:
            return FELL
        if np.linalg.norm(model.project(x) - fp) < converge_tol:
            return CONVERGED
    return UNDECIDED
