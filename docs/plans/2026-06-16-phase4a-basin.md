# Phase 4a: Basin of Attraction 可視化 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 受動歩行リミットサイクルの basin of attraction を Poincaré 断面の 2D スライスとして可視化し、5モデル統制比較で「円弧足が basin を広げるか」を示す。

**Architecture:** 新規 `src/crane/basin.py` が stride 写像の前進反復で各断面点を「収束/転倒/未決」に分類する。stride は副作用なしの純関数なので、グリッド各点を `multiprocessing` で並列分類する（model はワーカ内で picklable な make_fn + params dataclass から構築）。`viz.plot_basin` が描画、`scripts/basin_map.py`（単体）と `scripts/basin_compare.py`（5モデル比較）が実行口。

**Tech Stack:** numpy, scipy（既存）, matplotlib（既存）, multiprocessing（標準ライブラリ）。新依存なし。

**設計ドキュメント:** `docs/2026-06-16-basin-of-attraction-design.md`

---

## File Structure

- Create: `src/crane/basin.py` — 分類ロジック `classify_ic`、スライス計算 `basin_slice`、`BasinResult` dataclass、並列ワーカ。
- Modify: `src/crane/viz.py` — `plot_basin(result, out)` 追加。
- Modify: `src/crane/references.py` — Schwab & Wisse 2001 の simplest basin provenance 定数を追加。
- Create: `scripts/basin_map.py` — 1モデルの basin 地図 + meta.json。
- Create: `scripts/basin_compare.py` — 5モデル一括、比較 montage + basin_fraction 表。
- Create: `tests/test_basin_classify.py` — `classify_ic` の不変量。
- Create: `tests/test_basin_slice.py` — `basin_slice` の形状・並列一致・basin_fraction。
- Create: `tests/test_basin_stability_gate.py` — 内部整合ゲート（不動点近傍球が全 CONVERGED ＝ max|λ|<1 と整合）。
- Create: `tests/test_basin_simplest_literature.py` — Schwab & Wisse 照合ゲート。
- Create: `tests/test_viz_basin.py` — `plot_basin` の PNG 出力スモーク。

**重要な設計判断（設計 doc からの精緻化）:** 設計 doc では `basin_slice(model, ...)` としたが、`multiprocessing` は HybridModel（lambda を保持）を pickle できない。よって `basin_slice(make_fn, params, ...)` に変更し、ワーカ初期化子が各プロセスで `make_fn(params)` を一度だけ構築する。`make_fn`（モジュールレベル関数）と `params`（frozen dataclass）はどちらも picklable。

---

## Task 1: 分類ロジック `classify_ic`

**Files:**
- Create: `src/crane/basin.py`
- Test: `tests/test_basin_classify.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_basin_classify.py
import numpy as np

from crane.basin import CONVERGED, FELL, classify_ic
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane import references as ref


def _simplest_fp():
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(
        model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    )
    assert fp.converged
    return model, fp.y


def test_fixed_point_classifies_converged():
    model, y_star = _simplest_fp()
    assert classify_ic(model, y_star, y_star, max_strides=20) == CONVERGED


def test_far_point_classifies_fell():
    model, y_star = _simplest_fp()
    far = y_star + np.array([0.5, 0.5])  # basin 外（転倒）
    assert classify_ic(model, far, y_star, max_strides=20) == FELL
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_basin_classify.py -v`
Expected: FAIL（`ModuleNotFoundError: crane.basin`）

- [ ] **Step 3: Write minimal implementation**

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_basin_classify.py -v`
Expected: PASS（2 件）

- [ ] **Step 5: Commit**

```bash
git add src/crane/basin.py tests/test_basin_classify.py
git commit -m "feat: basin classify_ic (converged/fell/undecided via stride iteration)"
```

---

## Task 2: スライス計算 `basin_slice` + `BasinResult`（並列）

**Files:**
- Modify: `src/crane/basin.py`
- Test: `tests/test_basin_slice.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_basin_slice.py
import numpy as np

from crane.basin import CONVERGED, BasinResult, basin_slice
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane import references as ref


def _simplest_fp():
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(
        model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    )
    assert fp.converged
    return fp.y


def test_basin_slice_shape_and_center():
    y_star = _simplest_fp()
    res = basin_slice(
        make_simplest,
        SimplestParams(gamma=ref.GAMMA_REF),
        y_star,
        axes=(0, 1),
        half_widths=(0.002, 0.002),
        resolution=5,
        model_name="simplest",
        n_workers=1,
    )
    assert isinstance(res, BasinResult)
    assert res.grid.shape == (5, 5)
    # 奇数解像度なので中心セルは不動点そのもの → CONVERGED
    assert res.grid[2, 2] == CONVERGED
    assert 0.0 < res.basin_fraction <= 1.0


def test_basin_slice_serial_parallel_agree():
    y_star = _simplest_fp()
    kw = dict(
        axes=(0, 1),
        half_widths=(0.01, 0.01),
        resolution=5,
        model_name="simplest",
    )
    serial = basin_slice(
        make_simplest, SimplestParams(gamma=ref.GAMMA_REF), y_star, n_workers=1, **kw
    )
    parallel = basin_slice(
        make_simplest, SimplestParams(gamma=ref.GAMMA_REF), y_star, n_workers=2, **kw
    )
    assert np.array_equal(serial.grid, parallel.grid)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_basin_slice.py -v`
Expected: FAIL（`ImportError: cannot import name 'basin_slice'`）

- [ ] **Step 3: Write minimal implementation**

`src/crane/basin.py` に追記（`from __future__` 行の直後の import 群に `dataclass`, `Pool` を追加し、末尾に以下）:

```python
# --- ファイル先頭の import 群に追加 ---
from collections.abc import Callable
from dataclasses import dataclass
from multiprocessing import Pool

# --- ファイル末尾に追加 ---
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_basin_slice.py -v`
Expected: PASS（2 件）

- [ ] **Step 5: Commit**

```bash
git add src/crane/basin.py tests/test_basin_slice.py
git commit -m "feat: basin_slice 2D sweep with multiprocessing (picklable make_fn+params)"
```

---

## Task 3: 内部整合ゲート（固有値との接続）

**Files:**
- Test: `tests/test_basin_stability_gate.py`

不動点の十分小さい近傍球が全 CONVERGED であることを確認する。これは Phase 1 で
証明済みの `max|λ|<1`（simplest long-period の局所安定）と basin が矛盾しないことの
自動ゲート。

- [ ] **Step 1: Write the failing test**

```python
# tests/test_basin_stability_gate.py
import numpy as np

from crane.basin import CONVERGED
from crane.basin import basin_slice
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane import references as ref


def test_small_neighborhood_all_converges_simplest():
    """局所漸近安定（max|λ|=0.589<1）なら不動点の微小近傍は全て収束する。

    basin が固有値解析と整合することの内部ゲート。
    """
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(
        model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    )
    assert fp.converged
    assert np.max(np.abs(fp.eigenvalues)) < 1.0  # 前提（Phase 1 ゲート）

    # 不動点の |y*|·0.5% 程度の微小球
    hw = 0.005 * np.abs(fp.y)
    res = basin_slice(
        make_simplest,
        SimplestParams(gamma=ref.GAMMA_REF),
        fp.y,
        axes=(0, 1),
        half_widths=(float(hw[0]), float(hw[1])),
        resolution=5,
        model_name="simplest",
        n_workers=1,
    )
    assert np.all(res.grid == CONVERGED)
```

- [ ] **Step 2: Run test to verify it fails (then passes)**

Run: `uv run pytest tests/test_basin_stability_gate.py -v`
Expected: PASS（Task 1-2 実装済みなら即 PASS。これはゲートであり新規プロダクションコードは不要）

> 注: このタスクはテスト追加のみ。もし FAIL する場合は converge_tol/max_strides が
> 厳しすぎるか basin 判定にバグがある兆候なので、basin.py 側を調査する。

- [ ] **Step 3: Commit**

```bash
git add tests/test_basin_stability_gate.py
git commit -m "test: basin internal-consistency gate (small neighborhood converges, ties to max|lambda|<1)"
```

---

## Task 4: 描画 `plot_basin`

**Files:**
- Modify: `src/crane/viz.py`
- Test: `tests/test_viz_basin.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_viz_basin.py
import numpy as np

from crane.basin import BasinResult
from crane.viz import plot_basin


def test_plot_basin_writes_png(tmp_path):
    grid = np.array([[0, 1, 2], [0, 0, 1], [2, 1, 0]], dtype=int)
    res = BasinResult(
        grid=grid,
        ax0_vals=np.linspace(-0.1, 0.1, 3),
        ax1_vals=np.linspace(-0.1, 0.1, 3),
        axes=(0, 1),
        fixed_point=np.array([0.0, 0.0]),
        basin_fraction=float(np.mean(grid == 0)),
        model_name="dummy",
    )
    out = tmp_path / "basin.png"
    plot_basin(res, out)
    assert out.exists() and out.stat().st_size > 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_viz_basin.py -v`
Expected: FAIL（`ImportError: cannot import name 'plot_basin'`）

- [ ] **Step 3: Write minimal implementation**

`src/crane/viz.py` の末尾に追加（先頭の import に `from matplotlib.colors import ListedColormap` を追加）:

```python
# --- 先頭 import 群に追加 ---
from matplotlib.colors import ListedColormap  # noqa: E402

# --- ファイル末尾に追加（BasinResult は型注釈のみなので文字列注釈で循環 import 回避）---
def plot_basin(result, out: Path) -> None:
    """basin 分類グリッドを描画。CONVERGED=緑, FELL=赤, UNDECIDED=灰。"""
    cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])  # 0,1,2
    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [
        result.ax0_vals[0],
        result.ax0_vals[-1],
        result.ax1_vals[0],
        result.ax1_vals[-1],
    ]
    ax.imshow(
        result.grid,
        origin="lower",
        extent=extent,
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=2,
        interpolation="nearest",
    )
    a0, a1 = result.axes
    ax.plot(result.fixed_point[a0], result.fixed_point[a1], "k*", markersize=12)
    ax.set_xlabel(f"section[{a0}]")
    ax.set_ylabel(f"section[{a1}]")
    ax.set_title(
        f"{result.model_name}  basin_fraction={result.basin_fraction:.3f}"
    )
    fig.savefig(out, dpi=150)
    plt.close(fig)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_viz_basin.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/crane/viz.py tests/test_viz_basin.py
git commit -m "feat: plot_basin (converged/fell/undecided color map with fixed point marker)"
```

---

## Task 5: 単体スクリプト `scripts/basin_map.py`

**Files:**
- Create: `scripts/basin_map.py`

5モデルを名前で選び、不動点を求めて basin スライスを計算・描画・meta 保存する CLI。

- [ ] **Step 1: モデルレジストリとスクリプト本体を書く**

```python
# scripts/basin_map.py
"""Phase 4a: 単体モデルの basin of attraction 地図を生成。

Usage: uv run python scripts/basin_map.py <model> [--resolution 120] [--workers N]
  model: simplest | compass | kneed | rocker_compass | rocker_kneed
"""

import argparse
import json
import sys

import numpy as np

from crane.basin import CONVERGED, FELL, UNDECIDED, basin_slice
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.viz import plot_basin


def _registry():
    """name -> (make_fn, params, guess, axes, half_widths) を返す。

    half_widths は実装時に basin を窓に収めるよう較正（収束領域が窓端に
    接していたら広げる）。下記は初期値。
    """
    from crane import references as ref_s
    from crane import references_goswami as ref_c
    from crane import references_kneed as ref_k
    from crane import references_mcgeer as ref_rc
    from crane import references_mcgeer_knees as ref_rk
    from crane.models.compass import CompassParams, make_compass
    from crane.models.kneed import KneedParams, make_kneed
    from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
    from crane.models.rocker_kneed import RockerKneedParams, make_rocker_kneed
    from crane.models.simplest import SimplestParams, make_simplest

    return {
        "simplest": (
            make_simplest,
            SimplestParams(gamma=ref_s.GAMMA_REF),
            np.array([ref_s.LONG_PERIOD_THETA, ref_s.LONG_PERIOD_THETA_DOT]),
            (0, 1),
            (0.08, 0.08),
        ),
        "compass": (
            make_compass,
            CompassParams(
                m=ref_c.M_LEG, m_h=ref_c.M_HIP, a=ref_c.A, b=ref_c.B,
                gamma=ref_c.GAMMA_GAIT, g=ref_c.G,
            ),
            np.array(ref_c.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
        "kneed": (
            make_kneed,
            KneedParams(
                m_h=ref_k.M_HIP, m_t=ref_k.M_THIGH, m_s=ref_k.M_SHANK,
                l_t=ref_k.L_THIGH, l_s=ref_k.L_SHANK, b_t=ref_k.B_THIGH,
                b_s=ref_k.B_SHANK, gamma=ref_k.GAMMA_GAIT, g=ref_k.G,
            ),
            np.array(ref_k.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
        "rocker_compass": (
            make_rocker_compass,
            RockerCompassParams(
                m=ref_rc.M_LEG, m_h=ref_rc.M_HIP, c=ref_rc.C_HIP_TO_COM,
                rho=ref_rc.RHO_GYR, L=ref_rc.L_LEG, R=ref_rc.R_FOOT,
                gamma=ref_rc.GAMMA_GAIT, g=ref_rc.G,
            ),
            np.array(ref_rc.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
        "rocker_kneed": (
            make_rocker_kneed,
            RockerKneedParams(
                m_h=ref_rk.M_HIP, m_t=ref_rk.M_THIGH, m_s=ref_rk.M_SHANK,
                l_t=ref_rk.L_THIGH, l_s=ref_rk.L_SHANK, b_t=ref_rk.B_THIGH,
                b_s=ref_rk.B_SHANK, R=ref_rk.R_FOOT, gamma=ref_rk.GAMMA_GAIT, g=ref_rk.G,
            ),
            np.array(ref_rk.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
    }


def main() -> None:
    reg = _registry()
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=sorted(reg))
    parser.add_argument("--resolution", type=int, default=120)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    make_fn, params, guess, axes, half_widths = reg[args.model]
    model = make_fn(params)
    fp = find_limit_cycle(model, guess)
    if not fp.converged:
        sys.exit(f"ERROR: limit cycle not found for {args.model}")
    print(f"{args.model}: y*={fp.y}  |lambda|={np.abs(fp.eigenvalues)}")

    res = basin_slice(
        make_fn, params, fp.y,
        axes=axes, half_widths=half_widths, resolution=args.resolution,
        model_name=args.model, n_workers=args.workers,
    )
    run_dir = new_run_dir(f"basin_{args.model}")
    plot_basin(res, run_dir / "basin.png")
    meta = {
        "model": args.model,
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": np.abs(fp.eigenvalues).tolist(),
        "axes": list(axes),
        "half_widths": list(half_widths),
        "resolution": args.resolution,
        "basin_fraction": res.basin_fraction,
        "counts": {
            "converged": int(np.sum(res.grid == CONVERGED)),
            "fell": int(np.sum(res.grid == FELL)),
            "undecided": int(np.sum(res.grid == UNDECIDED)),
        },
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"basin_fraction={res.basin_fraction:.3f}  outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 低解像度スモーク実行で動作確認**

Run: `uv run python scripts/basin_map.py simplest --resolution 40 --workers 4`
Expected: `basin_fraction=...` が表示され、`data/runs/.../basin.png` と `meta.json` が生成される。

- [ ] **Step 3: ruff format + check**

Run: `uv run ruff format scripts/basin_map.py && uv run ruff check scripts/basin_map.py`
Expected: クリーン（CI は `ruff format --check` を回すので format も必須）

- [ ] **Step 4: Commit**

```bash
git add scripts/basin_map.py
git commit -m "feat: basin_map.py single-model basin map CLI (5-model registry)"
```

---

## Task 6: Schwab & Wisse 2001 provenance + simplest 文献ゲート

**Files:**
- Modify: `src/crane/references.py`
- Test: `tests/test_basin_simplest_literature.py`

simplest walker の basin は **Schwab & Wisse 2001 "Basin of Attraction of the
Simplest Walking Model"（Proc. ASME DETC2001）** が一次資料。実装時に取り寄せて
provenance 付きで定数を記録し、本実装の basin がそれと整合することをゲートにする。

- [ ] **Step 1: 一次資料を取り寄せて basin の特徴を記録**

WebSearch/WebFetch で Schwab & Wisse 2001 を探し、basin に関する定量・定性記述を読む。
`src/crane/references.py` に provenance 付きで追記する（**数値は必ず原典から取得。
記憶からは書かない**）:

```python
# --- src/crane/references.py 末尾に追加 ---
# --- Basin of attraction (Schwab & Wisse 2001) ---
# Schwab, A.L. & Wisse, M. (2001) "Basin of Attraction of the Simplest
# Walking Model", Proc. ASME Design Engineering Technical Conferences,
# DETC2001/VIB-21363. URL: <取得時に記入>  取得日: 2026-06-16
#
# 記録する量（原典から取得できたもののみ。取得できなければ None のまま定性ゲート）:
#   BASIN_PROVENANCE: 出典メモ
#   BASIN_QUALITATIVE: basin 形状の定性的特徴（薄い連結領域 + フラクタル境界 等）
#   BASIN_AREA_REF: 公表された basin 面積系の数値（あれば。なければ None）
#   BASIN_AREA_REF_TOL: 照合許容（図読み取りなら緩める）
BASIN_PROVENANCE: str = "Schwab & Wisse 2001 DETC2001/VIB-21363 (取得日 2026-06-16)"
BASIN_AREA_REF: float | None = None  # 原典から取得できれば数値、なければ None
BASIN_AREA_REF_TOL: float = 0.30
```

> 実装者へ: 原典に basin の面積（または cell-mapping のセル数比など本実装の
> `basin_fraction` と比較可能な量）があれば `BASIN_AREA_REF` に記入し、URL も埋める。
> なければ `None` のままとし、ゲートは定性一致（薄い連結領域 + 微細なフラクタル境界）+
> 内部整合（Task 3）に委ねる。**取得できない数値を捏造しない。**

- [ ] **Step 2: 文献ゲートのテストを書く**

```python
# tests/test_basin_simplest_literature.py
import numpy as np

from crane import references as ref
from crane.basin import CONVERGED, FELL, basin_slice
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle


def _simplest_basin(resolution=60):
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(
        model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    )
    assert fp.converged
    return basin_slice(
        make_simplest,
        SimplestParams(gamma=ref.GAMMA_REF),
        fp.y,
        axes=(0, 1),
        half_widths=(0.08, 0.08),
        resolution=resolution,
        model_name="simplest",
        n_workers=None,
    )


def test_simplest_basin_is_thin_connected_region():
    """定性ゲート: basin は薄い連結領域（窓全体ではない）かつ非空。

    Schwab & Wisse 2001 の basin 図と整合（薄く、窓を埋め尽くさない）。
    """
    res = _simplest_basin()
    frac = res.basin_fraction
    assert 0.0 < frac < 0.95  # 非空かつ窓を埋め尽くさない（＝薄い）
    # CONVERGED と FELL が両方存在（境界が窓内にある）
    assert np.any(res.grid == CONVERGED) and np.any(res.grid == FELL)


def test_simplest_basin_area_matches_reference_if_available():
    """定量ゲート: 原典に basin 面積数値があれば本実装と照合（なければ skip）。"""
    if ref.BASIN_AREA_REF is None:
        import pytest

        pytest.skip("Schwab & Wisse の basin 面積数値が未取得（定性ゲートに委譲）")
    res = _simplest_basin()
    rel = abs(res.basin_fraction - ref.BASIN_AREA_REF) / ref.BASIN_AREA_REF
    assert rel < ref.BASIN_AREA_REF_TOL
```

- [ ] **Step 3: Run tests**

Run: `uv run pytest tests/test_basin_simplest_literature.py -v`
Expected: `test_..._thin_connected_region` PASS、`test_..._area_matches_reference_if_available` PASS または SKIP（数値未取得時）

- [ ] **Step 4: ruff + commit**

```bash
uv run ruff format src/crane/references.py tests/test_basin_simplest_literature.py
uv run ruff check src/crane/references.py tests/test_basin_simplest_literature.py
git add src/crane/references.py tests/test_basin_simplest_literature.py
git commit -m "test: simplest basin literature gate (Schwab & Wisse 2001 provenance)"
```

---

## Task 7: 比較スクリプト `scripts/basin_compare.py`（方向性主張）

**Files:**
- Create: `scripts/basin_compare.py`

5モデルの basin を計算し、統制ペア（compass: point vs rocker、kneed: point vs rocker）を
並べた montage と basin_fraction 表を出力。**rocker の basin_fraction > point-foot** を
チェック・報告する（軸3 の方向性主張）。

- [ ] **Step 1: スクリプト本体を書く**

```python
# scripts/basin_compare.py
"""Phase 4a: 5モデルの basin を比較。montage + basin_fraction 表を出力。

統制ペア（point-foot vs rocker-foot）で「円弧足が basin を広げるか」を検証する。
Usage: uv run python scripts/basin_compare.py [--resolution 120] [--workers N]
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from crane.basin import basin_slice  # noqa: E402
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle  # noqa: E402

from basin_map import _registry  # noqa: E402  同 scripts/ 内


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=120)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    reg = _registry()
    order = ["simplest", "compass", "rocker_compass", "kneed", "rocker_kneed"]
    results = {}
    for name in order:
        make_fn, params, guess, axes, half_widths = reg[name]
        fp = find_limit_cycle(make_fn(params), guess)
        if not fp.converged:
            raise SystemExit(f"ERROR: limit cycle not found for {name}")
        res = basin_slice(
            make_fn, params, fp.y, axes=axes, half_widths=half_widths,
            resolution=args.resolution, model_name=name, n_workers=args.workers,
        )
        results[name] = res
        print(f"{name}: basin_fraction={res.basin_fraction:.3f}")

    # montage
    cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])
    fig, axs = plt.subplots(1, len(order), figsize=(4 * len(order), 4))
    for ax, name in zip(axs, order):
        r = results[name]
        ax.imshow(
            r.grid, origin="lower",
            extent=[r.ax0_vals[0], r.ax0_vals[-1], r.ax1_vals[0], r.ax1_vals[-1]],
            aspect="auto", cmap=cmap, vmin=0, vmax=2, interpolation="nearest",
        )
        a0, a1 = r.axes
        ax.plot(r.fixed_point[a0], r.fixed_point[a1], "k*", markersize=10)
        ax.set_title(f"{name}\nfrac={r.basin_fraction:.3f}")
    fig.tight_layout()
    run_dir = new_run_dir("basin_compare")
    fig.savefig(run_dir / "basin_compare.png", dpi=150)
    plt.close(fig)

    # 方向性主張: rocker > point
    table = {name: results[name].basin_fraction for name in order}
    claims = {
        "compass_rocker_gt_point": table["rocker_compass"] > table["compass"],
        "kneed_rocker_gt_point": table["rocker_kneed"] > table["kneed"],
    }
    (run_dir / "basin_fractions.json").write_text(
        json.dumps({"fractions": table, "claims": claims}, indent=2)
    )
    print("claims:", claims)
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 低解像度スモーク実行**

Run: `uv run python scripts/basin_compare.py --resolution 30 --workers 4`
Expected: 5モデルの `basin_fraction` が表示され、`basin_compare.png` と
`basin_fractions.json` が生成。`claims` の真偽が出力される（低解像度では暫定）。

> 注: `claims` が False の場合、それは「円弧足が basin を広げない」可能性 or
> 窓幅の較正不足を意味する。half_widths を見直し（統制ペアは同一窓幅で比較）、
> 解像度を上げて再評価する。最終判断はぽんぽこ殿の目視（Task 8）。

- [ ] **Step 3: ruff + commit**

```bash
uv run ruff format scripts/basin_compare.py && uv run ruff check scripts/basin_compare.py
git add scripts/basin_compare.py
git commit -m "feat: basin_compare.py 5-model montage + rocker-vs-point fraction claim"
```

---

## Task 8: 本番解像度の実行・記録・歩容判定ゲート

**Files:**
- Modify: `README.md`
- Modify: `GOALS.md`

- [ ] **Step 1: 本番解像度で比較を実行**

Run: `uv run python scripts/basin_compare.py --resolution 120 --workers <cores>`
（数分/モデル）。`basin_compare.png` と `basin_fractions.json` を得る。
窓が basin を切り落としていたら half_widths を `_registry` で較正して再実行。

- [ ] **Step 2: 結果を README に記録**

`README.md` に Phase 4a セクションを追加（実測 basin_fraction 表 + 方向性主張の結果 +
simplest フラクタル basin の所見）。**数値は Step 1 の実測値で埋める。**

- [ ] **Step 3: GOALS に Phase 4a セクションを追加**

`GOALS.md` の Phase 3.6 の後に Phase 4a を追加。ゲート（内部整合・Schwab/Wisse 照合・
方向性主張・目視）のチェックリストと実測値を記録。最後の項目は未チェックで残す:

```markdown
- [ ] **ぽんぽこ殿の目視判定**: basin_compare.png が「円弧足で basin が広がる」主張を支持して見える
```

また従来 Phase 4 の記述を 4a（本件・完了）/ 4b（物理エンジン種戻し・未着手）に分割。

- [ ] **Step 4: 全テスト + ruff**

Run: `uv run pytest -q && uv run ruff format --check . && uv run ruff check .`
Expected: 全 green、ruff クリーン

- [ ] **Step 5: Commit**

```bash
git add README.md GOALS.md
git commit -m "docs: record phase 4a basin results and visual judgment gate"
```

- [ ] **Step 6: ぽんぽこ殿の目視判定を依頼**

`basin_compare.png` を提示し、「円弧足で basin が広がって見えるか」の判定を仰ぐ。
合格後に GOALS のチェックを埋め、PR へ（finishing-a-development-branch）。

---

## Self-Review チェック結果

- **Spec coverage:** 設計 doc の A〜E を全タスクで網羅。A スコープ→全体、B モジュール→Task1,2,4、C 分類/比較→Task1,2,7、D ゲート（内部整合→Task3、文献→Task6、方向性→Task7、目視→Task8）、E 出力/命名→Task5,7,8。
- **Placeholder scan:** Schwab & Wisse の数値のみ「実装時に原典取得」とした。これは CLAUDE.md（記憶からの数値禁止・provenance 必須）に従う意図的な設計で、取得不能時は定性ゲートに委譲する明示フォールバックを用意（捏造しない）。
- **Type consistency:** `BasinResult`（grid, ax0_vals, ax1_vals, axes, fixed_point, basin_fraction, model_name）、`basin_slice(make_fn, params, fixed_point, *, axes, half_widths, resolution, model_name, ...)`、`classify_ic(model, y0, fixed_point, *, max_strides, converge_tol)`、定数 `CONVERGED/FELL/UNDECIDED=0/1/2` が全タスクで一貫。
- **設計 doc からの精緻化:** `basin_slice` を `model` 引数 → `make_fn + params` 引数に変更（multiprocessing picklability）。理由を File Structure 節に明記。
