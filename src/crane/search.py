"""Poincaré 写像の不動点探索（Newton shooting）と安定性解析。"""

from dataclasses import dataclass, field

import numpy as np

from crane.models.simplest import SimplestParams, lift
from crane.stride import StrideError, stride


@dataclass
class FixedPoint:
    y: np.ndarray  # 断面座標 (θ*, θ̇*)
    eigenvalues: np.ndarray | None
    converged: bool
    history: list[tuple[np.ndarray, float]] = field(default_factory=list)


def poincare_map(p: SimplestParams, y: np.ndarray) -> np.ndarray:
    """断面 (heel-strike 直後) 上の reduced 2D 写像。"""
    x_end = stride(p, lift(y[0], y[1])).x_end
    return np.array([x_end[0], x_end[2]])


def _jacobian(p: SimplestParams, y: np.ndarray, h: float = 1e-7) -> np.ndarray:
    """中心差分による poincare_map の数値 Jacobian。"""
    n = y.size
    J = np.empty((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = h
        J[:, j] = (poincare_map(p, y + e) - poincare_map(p, y - e)) / (2.0 * h)
    return J


def _try_poincare(p: SimplestParams, y: np.ndarray) -> tuple[np.ndarray | None, float | None]:
    """StrideError を捕捉して (image, residual_norm) を返す。バスン外なら (None, None)。"""
    try:
        img = poincare_map(p, y)
        return img, float(np.linalg.norm(img - y))
    except StrideError:
        return None, None


def find_limit_cycle(
    p: SimplestParams,
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

    # 初期推測がバスン外なら、1回の Poincaré 反復で引き込んだ点を起点にする
    img, norm = _try_poincare(p, y)
    if img is None:
        # 単純反復 S(y): 1回の stride を起点として Newton に渡す
        # バスン外の場合は縮小率を落として近くの有効点を探す
        found = False
        for scale in [0.9, 0.75, 0.5, 0.25]:
            # 参照不動点近傍の O(γ^1/3) スケールの点に向けて補間
            anchor = np.array([p.gamma ** (1.0 / 3.0) * 0.97, -(p.gamma ** (1.0 / 3.0)) * 1.045])
            y_try = y * scale + anchor * (1.0 - scale)
            img_try, norm_try = _try_poincare(p, y_try)
            if img_try is not None:
                y = y_try
                img, norm = img_try, norm_try
                found = True
                break
        if not found:
            return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)

    assert img is not None and norm is not None
    history.append((y.copy(), norm))

    if norm < tol:
        J = _jacobian(p, y)
        return FixedPoint(y=y, eigenvalues=np.linalg.eigvals(J), converged=True, history=history)

    for _ in range(max_iter - 1):
        J = _jacobian(p, y)
        # Newton 更新: F(y) = S(y) - y = 0 → (J_S - I) * Δy = S(y) - y → y_new = y - Δy
        step = np.linalg.solve(J - np.eye(y.size), img - y)

        # バックトラッキングライン探索: バスン外のステップを縮小する
        alpha = 1.0
        y_new = y - alpha * step
        new_img, new_norm = _try_poincare(p, y_new)
        for _ in range(10):
            if new_img is not None:
                break
            alpha *= 0.5
            y_new = y - alpha * step
            new_img, new_norm = _try_poincare(p, y_new)
        else:
            return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)

        assert new_img is not None and new_norm is not None
        y = y_new
        img = new_img
        history.append((y.copy(), new_norm))

        if new_norm < tol:
            J = _jacobian(p, y)
            return FixedPoint(
                y=y, eigenvalues=np.linalg.eigvals(J), converged=True, history=history
            )

    return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)
