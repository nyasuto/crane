# src/crane/model.py
"""Hybrid dynamics モデルの共通インターフェース（多相対応）。"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PhaseSpec:
    """連続相 1 つ: 力学 + 終端イベント + 相末衝突写像。

    - dynamics(t, x) -> xdot
    - event_value(x) -> float   (g(x)=0 が相終端面)
    - event_accept(x) -> bool   (g=0 交差を相終端として受理するか)
    - impact(x) -> x'           (相終端の衝突写像。最終相では脚交換を含む)
    """

    dynamics: Callable[[float, np.ndarray], list | np.ndarray]
    event_value: Callable[[np.ndarray], float]
    event_accept: Callable[[np.ndarray], bool]
    impact: Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class HybridModel:
    """1 stride = phases を順に通過する相機械（パラメータ束縛済み）。

    - phases: 連続相の列。各相は event 受理時に impact を適用して次相へ
    - lift(y) -> x / project(x) -> y: Poincaré 断面（最終相の衝突直後）座標変換
    """

    phases: tuple[PhaseSpec, ...]
    lift: Callable[[np.ndarray], np.ndarray]
    project: Callable[[np.ndarray], np.ndarray]
