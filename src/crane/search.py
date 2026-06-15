# src/crane/search.py
"""Poincaré 写像の不動点探索（Newton shooting）と安定性解析。"""

from dataclasses import dataclass, field

import numpy as np

from crane.model import HybridModel
from crane.stride import StrideError, stride


@dataclass
class FixedPoint:
    y: np.ndarray  # 断面座標
    eigenvalues: np.ndarray | None
    converged: bool
    history: list[tuple[np.ndarray, float]] = field(default_factory=list)


def poincare_map(model: HybridModel, y: np.ndarray, *, n_strides: int = 1) -> np.ndarray:
    """断面 (heel-strike 直後) 上の reduced 写像（n_strides 歩合成）。"""
    x = model.lift(y)
    for _ in range(n_strides):
        x = stride(model, x).x_end
    return model.project(x)


def _try_poincare(
    model: HybridModel, y: np.ndarray, *, n_strides: int = 1
) -> tuple[np.ndarray | None, float | None]:
    try:
        img = poincare_map(model, y, n_strides=n_strides)
        return img, float(np.linalg.norm(img - y))
    except StrideError:
        return None, None


def _jacobian(
    model: HybridModel, y: np.ndarray, *, n_strides: int = 1, h: float = 1e-7
) -> np.ndarray:
    n = y.size
    J = np.empty((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = h
        J[:, j] = (
            poincare_map(model, y + e, n_strides=n_strides)
            - poincare_map(model, y - e, n_strides=n_strides)
        ) / (2.0 * h)
    return J


def find_limit_cycle(
    model: HybridModel,
    y_guess: np.ndarray,
    *,
    tol: float = 1e-12,
    max_iter: int = 30,
    n_strides: int = 1,
) -> FixedPoint:
    """Newton 法で S(y) = y を解く。

    初期推測またはニュートンステップがバスン外に出た場合はバックトラッキングで縮小する。
    収束履歴を必ず残す。
    """
    y = np.asarray(y_guess, dtype=float).copy()
    history: list[tuple[np.ndarray, float]] = []

    # 初期推測がバスン外なら即座に失敗を返す（呼び出し側がバスン内の推測を供給する責任）
    img, norm = _try_poincare(model, y, n_strides=n_strides)
    if img is None:
        return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)

    assert img is not None and norm is not None
    history.append((y.copy(), norm))

    if norm < tol:
        J = _jacobian(model, y, n_strides=n_strides)
        return FixedPoint(y=y, eigenvalues=np.linalg.eigvals(J), converged=True, history=history)

    c_armijo = 1e-4  # Armijo 十分減少係数
    for _ in range(max_iter - 1):
        cur_norm = float(np.linalg.norm(img - y))
        J = _jacobian(model, y, n_strides=n_strides)
        # Newton 更新: F(y) = S(y) - y = 0 → (J_S - I) * Δy = S(y) - y → y_new = y - Δy
        step = np.linalg.solve(J - np.eye(y.size), img - y)

        # バックトラッキングライン探索: バスン外、または残差が十分減少しない
        # ステップ（Armijo: ‖F(y_new)‖ ≤ (1 − c·α)‖F(y)‖ を満たさない）を縮小する。
        # これで残差増加ステップを排除し、収束をロバストにする (issue #1)。
        alpha = 1.0
        y_new = y - alpha * step
        new_img, new_norm = _try_poincare(model, y_new, n_strides=n_strides)
        for _ in range(10):
            if new_img is not None and new_norm <= (1.0 - c_armijo * alpha) * cur_norm:
                break
            alpha *= 0.5
            y_new = y - alpha * step
            new_img, new_norm = _try_poincare(model, y_new, n_strides=n_strides)
        else:
            return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)

        assert new_img is not None and new_norm is not None
        y = y_new
        img = new_img
        history.append((y.copy(), new_norm))

        if new_norm < tol:
            J = _jacobian(model, y, n_strides=n_strides)
            return FixedPoint(
                y=y, eigenvalues=np.linalg.eigvals(J), converged=True, history=history
            )

    return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)
