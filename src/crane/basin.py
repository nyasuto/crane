# src/crane/basin.py
"""Basin of attraction の分類とスライス計算。

各 Poincaré 断面点から stride 写像を前進反復し「収束 / 転倒 / 未決」を分類する。
stride 写像は副作用なしの純関数なので、グリッド各点は独立に並列計算できる。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from multiprocessing import Pool

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


@dataclass(frozen=True)
class BasinResult:
    grid: np.ndarray  # (resolution, resolution) 分類コード。row=ax1, col=ax0
    ax0_vals: np.ndarray  # 横軸（axes[0]）の値
    ax1_vals: np.ndarray  # 縦軸（axes[1]）の値
    axes: tuple[int, int]
    fixed_point: np.ndarray
    basin_fraction: float  # 窓内 CONVERGED 率（比較スカラー）
    model_name: str


# ワーカプロセスが保持する状態（make_fn + params から各プロセスで一度だけ構築）
_W: dict = {}


def _init_worker(make_fn, params, fixed_point, axes, max_strides, converge_tol):
    _W["model"] = make_fn(params)
    _W["fp"] = np.asarray(fixed_point, dtype=float)
    _W["axes"] = axes
    _W["max_strides"] = max_strides
    _W["converge_tol"] = converge_tol


def _classify_point(point: tuple[float, float]) -> int:
    v0, v1 = point
    y0 = _W["fp"].copy()
    a0, a1 = _W["axes"]
    y0[a0] = v0
    y0[a1] = v1
    try:
        return classify_ic(
            _W["model"],
            y0,
            _W["fp"],
            max_strides=_W["max_strides"],
            converge_tol=_W["converge_tol"],
        )
    except Exception:
        # 野生の初期条件では solve_ivp が StrideError 以外で落ちうる → basin 外扱い
        return FELL


def basin_slice(
    make_fn: Callable,
    params,
    fixed_point: np.ndarray,
    *,
    axes: tuple[int, int] = (0, 1),
    half_widths: tuple[float, float],
    resolution: int,
    model_name: str,
    max_strides: int = 20,
    converge_tol: float = 1e-3,
    n_workers: int | None = None,
) -> BasinResult:
    """fixed_point 中心、axes 以外の断面座標を不動点値に固定して 2D 掃引。

    make_fn(params) で各ワーカが model を構築（multiprocessing の picklability 対策）。
    n_workers=1 なら同一プロセスで直列実行。
    """
    fp = np.asarray(fixed_point, dtype=float)
    a0, a1 = axes
    ax0_vals = np.linspace(fp[a0] - half_widths[0], fp[a0] + half_widths[0], resolution)
    ax1_vals = np.linspace(fp[a1] - half_widths[1], fp[a1] + half_widths[1], resolution)
    # row-major: 外側 ax1（縦）, 内側 ax0（横）
    points = [(v0, v1) for v1 in ax1_vals for v0 in ax0_vals]
    init_args = (make_fn, params, fp, axes, max_strides, converge_tol)

    if n_workers == 1:
        _init_worker(*init_args)
        codes = [_classify_point(pt) for pt in points]
    else:
        with Pool(processes=n_workers, initializer=_init_worker, initargs=init_args) as pool:
            codes = pool.map(_classify_point, points)

    grid = np.array(codes, dtype=int).reshape(resolution, resolution)
    basin_fraction = float(np.mean(grid == CONVERGED))
    return BasinResult(
        grid=grid,
        ax0_vals=ax0_vals,
        ax1_vals=ax1_vals,
        axes=(a0, a1),
        fixed_point=fp,
        basin_fraction=basin_fraction,
        model_name=model_name,
    )
