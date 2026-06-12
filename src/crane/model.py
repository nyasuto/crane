# src/crane/model.py
"""Hybrid dynamics モデルの共通インターフェース。"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HybridModel:
    """1歩 = 連続相 + 衝突写像 の hybrid モデル（パラメータ束縛済み）。

    - dynamics(t, x) -> xdot
    - strike_value(x) -> float   (g(x)=0 が strike 面)
    - strike_accept(x) -> bool   (g=0 交差が本物の heel-strike か)
    - impact(x) -> x'            (衝突写像 + 脚ラベル交換)
    - lift(y) -> x               (Poincaré 断面座標 -> 全状態)
    - project(x) -> y            (全状態 -> 断面座標)
    """

    dynamics: Callable[[float, np.ndarray], list | np.ndarray]
    strike_value: Callable[[np.ndarray], float]
    strike_accept: Callable[[np.ndarray], bool]
    impact: Callable[[np.ndarray], np.ndarray]
    lift: Callable[[np.ndarray], np.ndarray]
    project: Callable[[np.ndarray], np.ndarray]
