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


def poincare_map(model: HybridModel, y: np.ndarray) -> np.ndarray:
    """断面 (heel-strike 直後) 上の reduced 写像。"""
    return model.project(stride(model, model.lift(y)).x_end)


def _try_poincare(model: HybridModel, y: np.ndarray) -> tuple[np.ndarray | None, float | None]:
    try:
        img = poincare_map(model, y)
        return img, float(np.linalg.norm(img - y))
    except StrideError:
        return None, None


def _jacobian(model: HybridModel, y: np.ndarray, h: float = 1e-7) -> np.ndarray:
    n = y.size
    J = np.empty((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = h
        J[:, j] = (poincare_map(model, y + e) - poincare_map(model, y - e)) / (2.0 * h)
    return J


def find_limit_cycle(
    model: HybridModel,
    y_guess: np.ndarray,
    *,
    tol: float = 1e-12,
    max_iter: int = 30,
) -> FixedPoint:
    """Newton 法で S(y) = y を解く。

    初期推測またはニュートンステップがバスン外に出た場合はバックトラッキングで縮小する。
    収束履歴を必ず残す。
    """
    y = np.asarray(y_guess, dtype=float).copy()
    history: list[tuple[np.ndarray, float]] = []

    # 初期推測がバスン外なら即座に失敗を返す（呼び出し側がバスン内の推測を供給する責任）
    img, norm = _try_poincare(model, y)
    if img is None:
        return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)

    assert img is not None and norm is not None
    history.append((y.copy(), norm))

    if norm < tol:
        J = _jacobian(model, y)
        return FixedPoint(y=y, eigenvalues=np.linalg.eigvals(J), converged=True, history=history)

    for _ in range(max_iter - 1):
        J = _jacobian(model, y)
        # Newton 更新: F(y) = S(y) - y = 0 → (J_S - I) * Δy = S(y) - y → y_new = y - Δy
        step = np.linalg.solve(J - np.eye(y.size), img - y)

        # バックトラッキングライン探索: バスン外のステップを縮小する
        alpha = 1.0
        y_new = y - alpha * step
        new_img, new_norm = _try_poincare(model, y_new)
        for _ in range(10):
            if new_img is not None:
                break
            alpha *= 0.5
            y_new = y - alpha * step
            new_img, new_norm = _try_poincare(model, y_new)
        else:
            return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)

        assert new_img is not None and new_norm is not None
        y = y_new
        img = new_img
        history.append((y.copy(), new_norm))

        if new_norm < tol:
            J = _jacobian(model, y)
            return FixedPoint(
                y=y, eigenvalues=np.linalg.eigvals(J), converged=True, history=history
            )

    return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)
