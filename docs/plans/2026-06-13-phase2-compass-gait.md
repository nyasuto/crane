# Crane Phase 2: sympy 導出レイヤー + Compass Gait Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** sympy による運動方程式・衝突写像の自動導出レイヤーを構築し、Goswami 点足 compass gait のリミットサイクルを発見・文献検証し、slope continuation で period-doubling を観察する。

**Architecture:** 汎用 `HybridModel` プロトコルで stride/search を model 非依存に refactor。derive レイヤー（Euler-Lagrange + 角運動量衝突写像）で compass を記号導出し、**退化パラメータ (a=0, m→0) で Phase 1 の検証済み Garcia 実装と数値一致させることでレイヤーを検証**する。文献ゲートは Phase 1 と同じ provenance 規律。

**Tech Stack:** Python 3.12 / uv / numpy / scipy / sympy (新規依存、設計書で承認済み) / matplotlib / pytest / ruff

**Spec:** `docs/2026-06-12-crane-design.md` Phase 2 節。モデル系統 = 点足 Goswami 系（ぽんぽこ殿決定 2026-06-13）

---

## 物理モデル早見表（Goswami 点足 compass gait）

- 構成: 剛体2脚（点足）、脚長 l = a + b。各脚の質量 m は足から a（hip から b）の点に集中。hip に点質量 m_H。
- 座標: `x = [θ_st, θ_sw, θ̇_st, θ̇_sw]`。slope 法線基準の絶対角、進行方向 +x、stance 足原点の pinned パラメータ化。
- 重力: slope frame で方向 ĝ = (sin γ, −cos γ)（Phase 1 と同じ「坂は重力を傾けて表現」）。
- Strike 面: swing 足高さ 0 ⇔ `g(x) = θ_st + θ_sw = 0`（非自明側）。
- 受理条件（Phase 1 の教訓を継承）: `θ_st < 0` かつ swing 足が下降 `ẏ_sw = −l sin(θ_st)(θ̇_st + θ̇_sw) < 0` ⇔ θ_st<0 では `θ̇_st + θ̇_sw < 0`。
- 衝突（plastic, 角運動量保存2本）:
  1. **系全体**の角運動量を「新接地点（swing 足先）」回りで保存
  2. **trailing 脚**（旧 stance 脚 = 新 swing 脚）の角運動量を hip 回りで保存
  - pre は pre-pinned パラメータ化で、post は post-pinned パラメータ化（同形の式に swap 代入）で評価する。これが正しい定式化（下記 Task 4 にコードで明記）。
- Poincaré 断面 = 衝突直後。`θ_sw = −θ_st` が成立するため自由座標は **y = (θ_st, θ̇_st, θ̇_sw) の3次元**。
- Garcia 1998 への退化: a=0（脚質量を足に置く）, m/m_H → 0, l=1, g=1 で Phase 1 の simplest walker と厳密に一致するはず。**これがレイヤー検証ゲート**。
- Goswami 文献の標準パラメータと公表 gait 数値は Task 6 で原典から取得（このプランには書かない。Phase 1 と同じ規律）。

## File Structure

| パス | 責務 |
|---|---|
| `src/crane/model.py` | `HybridModel` プロトコル（新規） |
| `src/crane/stride.py` | model 非依存に refactor |
| `src/crane/search.py` | model 非依存に refactor + `n_strides`（period-2 用、Task 8） |
| `src/crane/models/simplest.py` | `make_simplest(p) -> HybridModel` 追加 |
| `src/crane/derive/lagrange.py` | Euler-Lagrange 自動導出（新規） |
| `src/crane/derive/impact.py` | 角運動量・衝突写像導出（新規） |
| `src/crane/models/compass.py` | Goswami compass（derive レイヤー使用、新規） |
| `src/crane/references_goswami.py` | Goswami 文献値（provenance 付き、Task 6 で記入） |
| `src/crane/viz.py` | `angles_of` パラメータで model 非依存化 |
| `scripts/walk_compass.py` | Phase 2 デモ |
| `scripts/bifurcation_compass.py` | slope continuation + 分岐図 |

---

### Task 1: HybridModel プロトコルへの refactor（挙動変更なし）

**Files:**
- Create: `src/crane/model.py`
- Modify: `src/crane/stride.py`, `src/crane/search.py`, `src/crane/models/simplest.py`, `tests/test_stride.py`, `tests/test_search.py`, `tests/test_scaling.py`, `scripts/walk_simplest.py`

- [ ] **Step 1: model.py を書く**

```python
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
```

- [ ] **Step 2: stride.py を model 非依存に書き換える**

`stride(p, x0)` → `stride(model: HybridModel, x0)`。simplest 固有の import・event 関数・受理条件を削除し、以下の形にする（既存の burn-in / scuffing / 重複除去 / t_max ガードのロジックは維持）:

```python
# src/crane/stride.py
"""stride 写像: 衝突直後状態から積分 → heel-strike → 衝突写像 → 次の衝突直後状態。

strike 面 g(x)=0 は1歩の間に複数回交差しうる（foot scuffing）。本物の
heel-strike の判定は model.strike_accept に委譲する（Phase 1 の教訓:
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
    x_strike: np.ndarray  # 衝突写像適用前（heel-strike 瞬間の状態）
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
    """1歩分の写像。副作用なし。"""

    def event(t, x):
        return model.strike_value(x)

    event.terminal = True
    event.direction = 0

    t0 = 0.0
    x = np.asarray(x0, dtype=float)
    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []

    def append(t_seg, x_seg):
        if ts:  # 2セグメント目以降は境界点の重複を落とす
            t_seg, x_seg = t_seg[1:], x_seg[:, 1:]
        ts.append(t_seg)
        xs.append(x_seg)

    while t0 < t_max:
        burn = solve_ivp(model.dynamics, (t0, t0 + T_BURN), x, rtol=rtol, atol=atol)
        append(burn.t, burn.y)
        t0, x = burn.t[-1], burn.y[:, -1]
        if t0 >= t_max:
            break

        sol = solve_ivp(model.dynamics, (t0, t_max), x, events=event, rtol=rtol, atol=atol)
        append(sol.t, sol.y)

        if sol.t_events[0].size == 0:
            break

        t_e = sol.t_events[0][0]
        x_e = sol.y_events[0][0]
        if model.strike_accept(x_e):
            return StrideResult(
                x_end=model.impact(x_e),
                x_strike=x_e,
                t_step=t_e,
                t=np.concatenate(ts),
                x=np.concatenate(xs, axis=1),
            )
        t0, x = t_e, x_e  # scuffing: 無視して続行

    raise StrideError(f"no heelstrike before t_max={t_max}")
```

- [ ] **Step 3: simplest.py に make_simplest を追加**

既存関数は変更しない。末尾に追加（ẏ_sw = sin(θ)·ġ の幾何コメントを stride.py から移す）:

```python
def make_simplest(p: SimplestParams) -> HybridModel:
    """パラメータ束縛済みの HybridModel を返す。

    受理条件: θ<0 かつ ġ>0。交差点では ẏ_sw = sin(θ)·ġ が成り立ち、
    θ<0 で ġ>0 ⇔ swing 足が下降して接地（scuff-dip 出口の上昇交差を除外）。
    """
    return HybridModel(
        dynamics=lambda t, x: dynamics(t, x, p),
        strike_value=lambda x: x[1] - 2.0 * x[0],
        strike_accept=lambda x: x[0] < 0.0 and (x[3] - 2.0 * x[2]) > 0.0,
        impact=heelstrike_map,
        lift=lambda y: lift(y[0], y[1]),
        project=lambda x: np.array([x[0], x[2]]),
    )
```

（import に `from crane.model import HybridModel` を追加）

- [ ] **Step 4: search.py を model 非依存に書き換える**

```python
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
```

`find_limit_cycle(model, y_guess, *, tol=1e-12, max_iter=30)` は現行ロジックそのまま、第1引数を `p` から `model` に置換（バックトラッキング含む）。

- [ ] **Step 5: テストとスクリプトの呼び出し側を更新**

- `tests/test_stride.py`: 冒頭を `MODEL = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))` とし、`stride(P, x0)` → `stride(MODEL, x0)` に全置換
- `tests/test_search.py` / `tests/test_scaling.py`: `find_limit_cycle(P, ...)` → `find_limit_cycle(make_simplest(P), ...)`（test_scaling は γ ごとに `make_simplest(SimplestParams(gamma=gamma))`）
- `scripts/walk_simplest.py`: 同様に model を作って渡す

- [ ] **Step 6: 全テストが green のまま、demo も動くことを確認**

Run: `cd /Users/yast/git/crane && uv run pytest -v && uv run ruff format . && uv run ruff check . && uv run python scripts/walk_simplest.py --strides 3`
Expected: 19 passed（挙動変更なし）、demo 完走

- [ ] **Step 7: Commit**

```bash
git add -A && git commit -m "refactor: generalize stride/search over HybridModel protocol"
```

---

### Task 2: derive/lagrange.py — Euler-Lagrange 自動導出

**Files:**
- Create: `src/crane/derive/__init__.py`, `src/crane/derive/lagrange.py`
- Test: `tests/test_derive_lagrange.py`
- Modify: `pyproject.toml`（sympy 追加）

- [ ] **Step 1: sympy を追加**

Run: `cd /Users/yast/git/crane && uv add sympy`
（設計書 Phase 2 で承認済みの依存）

- [ ] **Step 2: 失敗するテストを書く**

```python
# tests/test_derive_lagrange.py
import numpy as np
import sympy as sp

from crane.derive.lagrange import derive_qdd


def test_pendulum_eom():
    """単振子: θ̈ = −(g/l) sin θ を解析解として照合。"""
    th, w, m, l, g = sp.symbols("th w m l g", positive=True)
    T = m * l**2 * w**2 / 2
    V = -m * g * l * sp.cos(th)
    qdd = derive_qdd([th], [w], T, V)
    expected = -(g / l) * sp.sin(th)
    assert sp.simplify(qdd[0] - expected) == 0


def test_free_particle_2d():
    """自由粒子: ẍ=0, ÿ=−g。"""
    x, y, vx, vy, m, g = sp.symbols("x y vx vy m g", positive=True)
    T = m * (vx**2 + vy**2) / 2
    V = m * g * y
    qdd = derive_qdd([x, y], [vx, vy], T, V)
    assert sp.simplify(qdd[0]) == 0
    assert sp.simplify(qdd[1] + g) == 0


def test_lambdify_numeric():
    """lambdify した EOM が数値評価できる。"""
    th, w, g, l = sp.symbols("th w g l", positive=True)
    T = l**2 * w**2 / 2
    V = -g * l * sp.cos(th)
    qdd = derive_qdd([th], [w], T, V)
    f = sp.lambdify((th, w, g, l), qdd[0], "numpy")
    assert np.isclose(f(0.3, 0.0, 1.0, 1.0), -np.sin(0.3))
```

- [ ] **Step 3: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_derive_lagrange.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.derive'`

- [ ] **Step 4: lagrange.py を実装**

`src/crane/derive/__init__.py` は空ファイル。

```python
# src/crane/derive/lagrange.py
"""Euler-Lagrange 自動導出（scleronomic 系、q/q̇ は plain symbols）。"""

import sympy as sp


def derive_qdd(
    q: list[sp.Symbol],
    qd: list[sp.Symbol],
    T: sp.Expr,
    V: sp.Expr,
    *,
    simplify: bool = True,
) -> sp.Matrix:
    """M(q)·q̈ = ∂L/∂q − (∂p/∂q)·q̇ を解いて q̈ の式ベクトルを返す。

    p = ∂L/∂q̇（一般化運動量）。時間陽依存なし（scleronomic）前提:
    d/dt p = (∂p/∂q)·q̇ + (∂p/∂q̇)·q̈ と展開できる。
    """
    L = T - V
    p_vec = sp.Matrix([sp.diff(L, v) for v in qd])
    M = p_vec.jacobian(qd)  # = Hessian_q̇(T)
    rhs = sp.Matrix([sp.diff(L, qi) for qi in q]) - p_vec.jacobian(q) * sp.Matrix(qd)
    qdd = M.LUsolve(rhs)
    return sp.simplify(qdd) if simplify else qdd
```

- [ ] **Step 5: テスト確認 + Commit**

Run: `uv run pytest tests/test_derive_lagrange.py -v` → 3 PASSED、`uv run pytest -v` → 22 passed、ruff format/check

```bash
git add pyproject.toml uv.lock src/crane/derive/ tests/test_derive_lagrange.py
git commit -m "feat: Euler-Lagrange symbolic derivation layer (sympy)"
```

---

### Task 3: derive/impact.py — 角運動量と衝突写像の導出

**Files:**
- Create: `src/crane/derive/impact.py`
- Test: `tests/test_derive_impact.py`

- [ ] **Step 1: 失敗するテストを書く**

解析解: 原点 pivot の質点（質量 m=2, r=(0,1), 角速度 θ̇=1, v=(1,0)）が
新 pivot P=(1,0) に乗り移る plastic impact。新 pivot 回り角運動量保存:
L_pre = m·cross(r−P, v) = 2·((−1)·0 − 1·1) = −2。post は半径² = |r−P|² = 2 の
回転: L_post = m·ρ²·ψ̇ = 4ψ̇ ⇒ ψ̇ = −0.5。

```python
# tests/test_derive_impact.py
import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum, cross2d


def test_cross2d():
    a = sp.Matrix([1, 0])
    b = sp.Matrix([0, 2])
    assert cross2d(a, b) == 2


def test_angular_momentum_point_mass():
    """質点の pivot 回り角運動量 L = m l² θ̇。"""
    th, w, m, l = sp.symbols("th w m l", positive=True)
    r = l * sp.Matrix([sp.sin(th), sp.cos(th)])
    L = angular_momentum([(m, r)], [th], [w], sp.Matrix([0, 0]))
    assert sp.simplify(L - (-m * l**2 * w)) == 0 or sp.simplify(L - m * l**2 * w) == 0
    # 符号は回転規約依存。絶対値が m l² w であることを数値で確認:
    f = sp.lambdify((th, w, m, l), L, "numpy")
    assert np.isclose(abs(f(0.3, 1.0, 2.0, 1.5)), 2.0 * 1.5**2)


def test_pivot_transfer_analytic():
    """pivot 乗り移り衝突の解析解 ψ̇ = −0.5 を再現。"""
    th, w = sp.symbols("th w")
    m_val = 2.0
    r = sp.Matrix([sp.sin(th), sp.cos(th)])  # l=1
    P = sp.Matrix([1, 0])
    L_pre = angular_momentum([(m_val, r)], [th], [w], P)
    # post: P 回りの回転、角度 psi、半径 rho
    psi, wp = sp.symbols("psi wp")
    rho = sp.sqrt(2)
    r_post = P + rho * sp.Matrix([-sp.sin(psi), sp.cos(psi)])
    L_post = angular_momentum([(m_val, r_post)], [psi], [wp], P)
    # 保存則を wp について解く
    sol = sp.solve(sp.Eq(L_post, L_pre), wp)[0]
    f = sp.lambdify((th, w, psi), sol, "numpy")
    # th=0 (r=(0,1)), w=1; post 配置は同一点: P + sqrt(2)(−sin ψ, cos ψ) = (0,1) → ψ=π/4
    assert np.isclose(f(0.0, 1.0, np.pi / 4), -0.5)
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_derive_impact.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: impact.py を実装**

```python
# src/crane/derive/impact.py
"""点質量系の角運動量と、保存則に基づく plastic impact 写像の導出部品。"""

import sympy as sp


def cross2d(a: sp.Matrix, b: sp.Matrix) -> sp.Expr:
    """2D 外積（スカラー）。"""
    return a[0] * b[1] - a[1] * b[0]


def angular_momentum(
    bodies: list[tuple[sp.Expr, sp.Matrix]],
    q: list[sp.Symbol],
    qd: list[sp.Symbol],
    about: sp.Matrix,
) -> sp.Expr:
    """点質量系の点 about 回りの角運動量。bodies = [(mass, r(q))]。

    v_i = (∂r_i/∂q)·q̇ で速度を導出する。
    """
    total = sp.S.Zero
    for mass, r in bodies:
        v = r.jacobian(q) * sp.Matrix(qd)
        total += mass * cross2d(r - about, v)
    return sp.expand(total)
```

- [ ] **Step 4: テスト確認 + Commit**

Run: `uv run pytest tests/test_derive_impact.py -v` → 3 PASSED、全体 25 passed、ruff

```bash
git add src/crane/derive/impact.py tests/test_derive_impact.py
git commit -m "feat: angular momentum derivation for impact maps"
```

---

### Task 4: models/compass.py — Goswami compass の記号導出

**Files:**
- Create: `src/crane/models/compass.py`
- Test: `tests/test_compass.py`

- [ ] **Step 1: 失敗するテストを書く（不変量ベース）**

```python
# tests/test_compass.py
import numpy as np
from scipy.integrate import solve_ivp

from crane.models.compass import CompassParams, energy, kinetic_energy, make_compass

P = CompassParams(m=5.0, m_h=10.0, a=0.5, b=0.5, gamma=0.05)
MODEL = make_compass(P)


def test_equilibrium_aligned_with_gravity():
    """両脚が重力方向 (絶対角 γ) に揃った静止状態では加速度ゼロ。"""
    x = np.array([P.gamma, P.gamma, 0.0, 0.0])
    xdot = MODEL.dynamics(0.0, x)
    assert np.allclose(xdot, [0.0, 0.0, 0.0, 0.0], atol=1e-12)


def test_swing_phase_conserves_energy():
    """連続相は保存系: E = T + V が一定（受動歩行の本質）。"""
    x0 = np.array([0.2, -0.3, -1.0, 0.5])
    sol = solve_ivp(MODEL.dynamics, (0.0, 0.5), x0, rtol=1e-11, atol=1e-13, dense_output=True)
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)


def test_impact_dissipates_kinetic_energy():
    """plastic 衝突は運動エネルギーを増やさない（一般に厳密減少）。"""
    x_pre = np.array([-0.15, 0.15, -1.2, -0.8])  # strike 面上 (θ_st+θ_sw=0)、下降接地
    assert abs(MODEL.strike_value(x_pre)) < 1e-15
    assert MODEL.strike_accept(x_pre)
    x_post = MODEL.impact(x_pre)
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_impact_swaps_leg_labels():
    """衝突で配置は脚ラベル交換: θ_st⁺=θ_sw⁻, θ_sw⁺=θ_st⁻。"""
    x_pre = np.array([-0.15, 0.15, -1.2, -0.8])
    x_post = MODEL.impact(x_pre)
    assert x_post[0] == x_pre[1]
    assert x_post[1] == x_pre[0]


def test_impact_zero_velocity_is_fixed():
    """静止状態の衝突は静止のまま（線形写像の自明解）。"""
    x_pre = np.array([-0.15, 0.15, 0.0, 0.0])
    x_post = MODEL.impact(x_pre)
    assert np.allclose(x_post[2:], [0.0, 0.0], atol=1e-12)


def test_section_lift_project_roundtrip():
    y = np.array([0.2, -1.1, 0.3])
    x = MODEL.lift(y)
    assert np.isclose(x[1], -x[0])
    assert np.allclose(MODEL.project(x), y)
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_compass.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.models.compass'`

- [ ] **Step 3: compass.py を実装**

```python
# src/crane/models/compass.py
"""Goswami 点足 compass gait。derive レイヤーで記号導出し lambdify。

座標: x = [θ_st, θ_sw, θ̇_st, θ̇_sw]（slope 法線基準の絶対角、進行 +x）。
stance 足原点の pinned パラメータ化。strike 面: g(x) = θ_st + θ_sw = 0。
受理: θ_st < 0 かつ swing 足下降（ẏ_sw = −l sinθ_st (θ̇_st+θ̇_sw) < 0、
θ_st<0 では θ̇_st+θ̇_sw < 0 と等価）。

衝突（plastic）: pre は pre-pinned 系で、post は post-pinned 系（脚ラベル
交換済みの同形の式）で評価した角運動量を等置する:
  1. 系全体: post の pivot 回り = pre の swing 足先回り（同一物理点 = 新接地点）
  2. trailing 脚: post の swing 脚 hip 回り = pre の stance 脚 hip 回り（同一物理脚）
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum
from crane.derive.lagrange import derive_qdd
from crane.model import HybridModel


@dataclass(frozen=True)
class CompassParams:
    m: float  # 脚質量 [kg]
    m_h: float  # hip 質量 [kg]
    a: float  # 足 → 脚質量点 距離 [m]
    b: float  # 脚質量点 → hip 距離 [m]
    gamma: float  # 斜面角 [rad]
    g: float = 9.81

    @property
    def l(self) -> float:
        return self.a + self.b


@lru_cache(maxsize=1)
def _build():
    """記号導出して lambdify した関数群を返す（モジュール内キャッシュ）。"""
    th_st, th_sw, w_st, w_sw = sp.symbols("th_st th_sw w_st w_sw")
    m, m_h, a, b, gamma, g = sp.symbols("m m_h a b gamma g", positive=True)
    length = a + b
    q, qd = [th_st, th_sw], [w_st, w_sw]

    def down(theta):
        """絶対角 theta の脚に沿って hip から下向きの単位ベクトル。"""
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    hip = sp.Matrix([-length * sp.sin(th_st), length * sp.cos(th_st)])
    p_st = hip + b * down(th_st)  # stance 脚質量点（足から a）
    p_sw = hip + b * down(th_sw)  # swing 脚質量点（hip から b）
    foot_sw = hip + length * down(th_sw)

    def vel(r):
        return r.jacobian(q) * sp.Matrix(qd)

    T = (
        m_h * vel(hip).dot(vel(hip))
        + m * vel(p_st).dot(vel(p_st))
        + m * vel(p_sw).dot(vel(p_sw))
    ) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])  # slope frame の重力方向
    V = -g * (m_h * hip.dot(g_dir) + m * p_st.dot(g_dir) + m * p_sw.dot(g_dir))

    qdd = derive_qdd(q, qd, T, V, simplify=False)
    args = (th_st, th_sw, w_st, w_sw, m, m_h, a, b, gamma, g)
    f_qdd = sp.lambdify(args, [qdd[0], qdd[1]], "numpy")
    f_energy = sp.lambdify(args, T + V, "numpy")
    f_kinetic = sp.lambdify(args, T, "numpy")

    # --- 衝突写像（docstring の定式化） ---
    bodies = [(m_h, hip), (m, p_st), (m, p_sw)]
    origin = sp.Matrix([0, 0])
    L_sys_pivot = angular_momentum(bodies, q, qd, origin)
    L_sys_swfoot = angular_momentum(bodies, q, qd, foot_sw)
    L_st_hip = angular_momentum([(m, p_st)], q, qd, hip)
    L_sw_hip = angular_momentum([(m, p_sw)], q, qd, hip)

    wp_st, wp_sw = sp.symbols("wp_st wp_sw")  # post 速度（post ラベル）
    swap = {th_st: th_sw, th_sw: th_st, w_st: wp_st, w_sw: wp_sw}
    eqs = [
        L_sys_pivot.subs(swap, simultaneous=True) - L_sys_swfoot,
        L_sw_hip.subs(swap, simultaneous=True) - L_st_hip,
    ]
    A, rhs = sp.linear_eq_to_matrix(eqs, [wp_st, wp_sw])
    qd_post = A.LUsolve(rhs)
    f_impact = sp.lambdify(args, [qd_post[0], qd_post[1]], "numpy")

    return f_qdd, f_energy, f_kinetic, f_impact


def _args(x, p: CompassParams):
    return (*x, p.m, p.m_h, p.a, p.b, p.gamma, p.g)


def dynamics(t: float, x, p: CompassParams):
    f_qdd, _, _, _ = _build()
    qdd = f_qdd(*_args(x, p))
    return [x[2], x[3], qdd[0], qdd[1]]


def energy(x, p: CompassParams) -> float:
    _, f_energy, _, _ = _build()
    return float(f_energy(*_args(x, p)))


def kinetic_energy(x, p: CompassParams) -> float:
    _, _, f_kinetic, _ = _build()
    return float(f_kinetic(*_args(x, p)))


def heelstrike_map(x: np.ndarray, p: CompassParams) -> np.ndarray:
    """衝突写像 + 脚ラベル交換。post 速度は post ラベルで返る。"""
    _, _, _, f_impact = _build()
    wp = f_impact(*_args(x, p))
    return np.array([x[1], x[0], float(wp[0]), float(wp[1])])


def make_compass(p: CompassParams) -> HybridModel:
    return HybridModel(
        dynamics=lambda t, x: dynamics(t, x, p),
        strike_value=lambda x: x[0] + x[1],
        strike_accept=lambda x: x[0] < 0.0 and (x[2] + x[3]) < 0.0,
        impact=lambda x: heelstrike_map(x, p),
        lift=lambda y: np.array([y[0], -y[0], y[1], y[2]]),
        project=lambda x: np.array([x[0], x[2], x[3]]),
    )
```

- [ ] **Step 4: テスト確認 + Commit**

Run: `uv run pytest tests/test_compass.py -v` → 6 PASSED、全体 31 passed、ruff format/check

（FAIL 時の切り分け: equilibrium が落ちる → V の重力方向の符号、energy drift → derive_qdd、
KE 増加 → 衝突の swap 定式化を疑う。各関数は独立に検証可能。）

```bash
git add src/crane/models/compass.py tests/test_compass.py
git commit -m "feat: Goswami compass gait via symbolic derivation layer"
```

---

### Task 5: Garcia 退化ゲート — レイヤーを Phase 1 実装で検証

**Files:**
- Test: `tests/test_compass_garcia_limit.py`

compass を a=0（脚質量を足に配置）, m=1e-9, m_h=1, l=1, g=1 に退化させると
Garcia 1998 simplest walker に一致するはず。Phase 1 の検証済み実装（文献照合 +
不動点6桁一致済み）を 'truth' として、導出レイヤー全体を end-to-end 検証する。

座標変換: Garcia (θ, φ) ↔ compass 絶対角 (θ_st, θ_sw) は θ=θ_st, φ=θ_st−θ_sw。
速度も同形。strike 面も一致（φ=2θ ⇔ θ_sw=−θ_st）。

- [ ] **Step 1: テストを書く**

```python
# tests/test_compass_garcia_limit.py
import numpy as np

from crane import references as ref
from crane.models import simplest
from crane.models.compass import CompassParams, make_compass
from crane.search import find_limit_cycle

P_DEG = CompassParams(m=1e-9, m_h=1.0, a=0.0, b=1.0, gamma=ref.GAMMA_REF, g=1.0)
MODEL = make_compass(P_DEG)
P_GARCIA = simplest.SimplestParams(gamma=ref.GAMMA_REF)


def to_garcia(x_abs):
    """compass 絶対角 → Garcia (θ, φ) 座標。"""
    return np.array([x_abs[0], x_abs[0] - x_abs[1], x_abs[2], x_abs[2] - x_abs[3]])


def to_abs(x_g):
    return np.array([x_g[0], x_g[0] - x_g[1], x_g[2], x_g[2] - x_g[3]])


def test_dynamics_reduces_to_garcia():
    """m→0, a=0 で連続相が Garcia EOM に一致（ランダム20状態）。"""
    rng = np.random.default_rng(0)
    for _ in range(20):
        x_g = rng.uniform([-0.3, -0.6, -2.0, -2.0], [0.3, 0.6, 2.0, 2.0])
        xdot_g = np.asarray(simplest.dynamics(0.0, x_g, P_GARCIA))
        xdot_abs = np.asarray(MODEL.dynamics(0.0, to_abs(x_g)))
        # 加速度を Garcia 座標に変換: θ̈=ẅ_st, φ̈=ẅ_st−ẅ_sw
        got = np.array([xdot_abs[2], xdot_abs[2] - xdot_abs[3]])
        assert np.allclose(got, xdot_g[2:], rtol=0, atol=1e-5)


def test_impact_reduces_to_garcia():
    """m→0, a=0 で衝突写像が Garcia の cos(2θ) 写像に一致。"""
    rng = np.random.default_rng(1)
    for _ in range(10):
        theta = rng.uniform(-0.3, -0.1)
        theta_dot = rng.uniform(-0.5, -0.1)
        phi_dot = rng.uniform(-0.5, 0.5)
        x_g_pre = np.array([theta, 2 * theta, theta_dot, phi_dot])
        expected = simplest.heelstrike_map(x_g_pre)
        got = to_garcia(MODEL.impact(to_abs(x_g_pre)))
        assert np.allclose(got, expected, rtol=0, atol=1e-5)


def test_fixed_point_reduces_to_garcia():
    """退化 compass の不動点・固有値が Phase 1 の結果に一致。"""
    x_g = simplest.lift(ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT)
    x_abs = to_abs(x_g)
    y_guess = np.array([x_abs[0], x_abs[2], x_abs[3]])
    fp = find_limit_cycle(MODEL, y_guess)
    assert fp.converged
    # Phase 1 実測の真の不動点 (0.2003109, -0.1998325) と照合
    assert np.isclose(fp.y[0], 0.2003109, atol=1e-5)
    assert np.isclose(fp.y[1], -0.1998325, atol=1e-5)
    mags = np.sort(np.abs(fp.eigenvalues))
    # 3D 断面の multiplier: 上位2つは Garcia の複素対 |λ|=0.5891、
    # 第3方向は m→0 で swing 速度が slave されるため ≈ 0
    assert np.allclose(mags[1:], [0.5891, 0.5891], atol=2e-3)
    assert mags[0] < 1e-3
```

- [ ] **Step 2: テスト実行**

Run: `uv run pytest tests/test_compass_garcia_limit.py -v`
Expected: 3 PASSED（Task 4 まで揃っていれば実装は不要。FAIL なら derive レイヤー
または compass の定式化にバグ — どのテストが落ちるかで連続相/衝突/全体を切り分け、
原因を調査して報告。**勝手に許容値を緩めない**）

- [ ] **Step 3: Commit**

```bash
git add tests/test_compass_garcia_limit.py
git commit -m "test: compass degenerates to verified Garcia implementation (layer gate)"
```

---

### Task 6: Goswami 文献値の取得と記録

**Files:**
- Create: `src/crane/references_goswami.py`
- Test: `tests/test_references_goswami.py`
- Modify: `CLAUDE.md`（references の記述を `src/crane/references*.py` に一般化）

**この Task は Web 調査を含む。** WebSearch / WebFetch で以下を入手・照合する:

1. **一次候補**: Goswami, Thuilot, Espiau — INRIA Research Report RR-2996 (1996)
   "Compass-Like Biped Robot Part I: Stability and Bifurcation of Passive Gaits"
   （inria.fr の HAL で公開されているはず）および/or IJRR 1998
   "A Study of the Passive Gait of a Compass-Like Biped Robot: Symmetry and Chaos"
2. 取得する情報:
   - 標準パラメータセット（m, m_H, a, b, 脚長、g）
   - passive gait が存在する slope と、リミットサイクルの公表数値
     （strike 時の角度・角速度、step period、固有値など、印刷されているもの全て）
   - period-doubling 開始 slope とカオス到達 slope（分岐図の記述）
   - **衝突写像・EOM の定式化を本プランの早見表・Task 4 実装と照合**
     （特に角運動量保存2本の取り方）。相違があれば報告（勝手に直さない）
3. 数表がなく図のみの量は「図からの読み取り値（digitized, 精度±X%）」と
   正直にラベルして記録する（Phase 1 の O(γ) 導出値と同じ流儀）

```python
# tests/test_references_goswami.py（先に書く: TDD）
from crane import references_goswami as ref


def test_goswami_references_are_filled():
    """文献値が provenance 付きで記入済み（isinstance 検査は Ellipsis 対策）。"""
    assert ref.PROVENANCE.startswith("Goswami")
    assert "http" in ref.PROVENANCE
    assert isinstance(ref.M_LEG, float)
    assert isinstance(ref.M_HIP, float)
    assert isinstance(ref.A, float)
    assert isinstance(ref.B, float)
    assert isinstance(ref.GAMMA_GAIT, float)  # 公表 gait のある slope [rad]
    assert 0.0 < ref.GAMMA_GAIT < 0.2
```

`references_goswami.py` の必須フィールド: `PROVENANCE`, `M_LEG`, `M_HIP`, `A`, `B`,
`G`, `GAMMA_GAIT`, gait 照合量（論文の記載形式に合わせ、例えば
`STEP_PERIOD: float | None`, `THETA_STRIKE: float | None`,
`SECTION_GUESS: tuple[float, float, float]`（断面座標での初期推定、出所コメント必須）,
`GAIT_TOLERANCE: float`（読み取り精度に応じた照合許容、根拠コメント必須）,
`FIRST_PERIOD_DOUBLING_GAMMA: float | None`, `CHAOS_GAMMA: float | None`）。
論文に無い量は None + 理由コメント。

- [ ] **Step 1: ゲートテストを書き、FAIL 確認**
- [ ] **Step 2: Web 調査を実施し references_goswami.py を実数で記入**
- [ ] **Step 3: CLAUDE.md の「src/crane/references.py」を「src/crane/references*.py」に修正**
- [ ] **Step 4: テスト + ruff + Commit**

```bash
git add src/crane/references_goswami.py tests/test_references_goswami.py CLAUDE.md
git commit -m "feat: record Goswami compass published values with provenance"
```

（EOM/衝突定式化の照合結果を commit body に `Formulation cross-check: <内容>` で記録）

---

### Task 7: Compass リミットサイクル発見 + 文献ゲート

**Files:**
- Test: `tests/test_compass_cycle.py`

- [ ] **Step 1: ゲートテストを書く**

```python
# tests/test_compass_cycle.py
import numpy as np

from crane import references_goswami as ref
from crane.models.compass import CompassParams, make_compass
from crane.search import find_limit_cycle

P = CompassParams(m=ref.M_LEG, m_h=ref.M_HIP, a=ref.A, b=ref.B, gamma=ref.GAMMA_GAIT, g=ref.G)
MODEL = make_compass(P)


def test_compass_limit_cycle_exists_at_published_slope():
    """Phase 2 ゲート: 公表パラメータ・slope でリミットサイクルが存在し安定。"""
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    assert np.max(np.abs(fp.eigenvalues)) < 1.0  # 公表 gait は安定（文献記述）


def test_compass_gait_matches_published_quantities():
    """公表されている gait 量と照合（許容は references_goswami の根拠付き値）。"""
    from crane.stride import stride

    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    result = stride(MODEL, MODEL.lift(fp.y))
    if ref.STEP_PERIOD is not None:
        assert np.isclose(result.t_step, ref.STEP_PERIOD, rtol=ref.GAIT_TOLERANCE)
    if ref.THETA_STRIKE is not None:
        assert np.isclose(abs(result.x_strike[0]), abs(ref.THETA_STRIKE), rtol=ref.GAIT_TOLERANCE)
```

- [ ] **Step 2: テスト実行（探索デバッグ込み）**

Run: `uv run pytest tests/test_compass_cycle.py -v`

`SECTION_GUESS` から収束しない場合のフォールバック（決定論的に）:
断面座標の粗グリッド（例: θ_st ∈ [0.05, 0.45] step 0.02、θ̇_st ∈ [-2.5, -0.5] step 0.1、
θ̇_sw ∈ [-2.0, 1.0] step 0.2）で `‖S(y)−y‖` 最小の点を求めて Newton の seed にする
小スクリプトを書いて実行し、見つかった seed を references_goswami.py の
SECTION_GUESS に「グリッド探索による補正値 (原典の図読み取りを起点)」とコメント
付きで更新する。それでも見つからなければ BLOCKED 報告（衝突定式化 or パラメータ解釈
の疑いを添えて）。

- [ ] **Step 3: 実測値の報告**

収束した不動点 y*、固有値、step period、strike 角度を必ず報告に含める
（Phase 3 / 分岐図の seed になる）。

- [ ] **Step 4: Commit**

```bash
git add tests/test_compass_cycle.py src/crane/references_goswami.py
git commit -m "feat: compass limit cycle discovery gated on Goswami published gait"
```

---

### Task 8: slope continuation + 分岐図（period-doubling 観察）

**Files:**
- Modify: `src/crane/search.py`（`n_strides` パラメータ追加）
- Create: `scripts/bifurcation_compass.py`
- Test: `tests/test_search_n_strides.py`, `tests/test_compass_continuation.py`

- [ ] **Step 1: search.py に n_strides を追加（period-2 探索用）**

`poincare_map(model, y, *, n_strides: int = 1)`: stride を n 回合成。
`find_limit_cycle(model, y_guess, *, tol=1e-12, max_iter=30, n_strides: int = 1)`:
内部の poincare_map / _try_poincare / _jacobian 呼び出しに n_strides を貫通させる。

テスト（simplest を使った安価で強い検証）:

```python
# tests/test_search_n_strides.py
import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle

MODEL = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))


def test_period1_is_also_fixed_point_of_s2():
    """period-1 不動点は S² の不動点でもあり、multiplier は2乗になる。"""
    fp1 = find_limit_cycle(MODEL, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    fp2 = find_limit_cycle(MODEL, fp1.y, n_strides=2)
    assert fp2.converged
    assert np.allclose(fp2.y, fp1.y, atol=1e-9)
    assert np.isclose(
        np.max(np.abs(fp2.eigenvalues)), np.max(np.abs(fp1.eigenvalues)) ** 2, rtol=1e-3
    )
```

- [ ] **Step 2: continuation テスト（公表安定域のみゲート）**

```python
# tests/test_compass_continuation.py
import numpy as np

from crane import references_goswami as ref
from crane.models.compass import CompassParams, make_compass
from crane.search import find_limit_cycle


def test_stable_family_continues_below_published_gait_slope():
    """公表 slope から下向きの単調 continuation で安定 family が続く。"""
    y = np.array(ref.SECTION_GUESS)
    gammas = np.linspace(ref.GAMMA_GAIT, ref.GAMMA_GAIT * 0.6, 5)
    for gamma in gammas:
        p = CompassParams(m=ref.M_LEG, m_h=ref.M_HIP, a=ref.A, b=ref.B, gamma=float(gamma), g=ref.G)
        fp = find_limit_cycle(make_compass(p), y)
        assert fp.converged, f"no convergence at gamma={gamma}"
        assert np.max(np.abs(fp.eigenvalues)) < 1.0, f"unstable at gamma={gamma}"
        y = fp.y
```

- [ ] **Step 3: bifurcation_compass.py を書く**

```python
# scripts/bifurcation_compass.py
"""slope continuation で compass gait family を追跡し分岐図を描く。

Usage: uv run python scripts/bifurcation_compass.py [--gamma-max-deg 6.0] [--step-deg 0.05]
"""

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from crane import references_goswami as ref  # noqa: E402
from crane.models.compass import CompassParams, make_compass  # noqa: E402
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle  # noqa: E402


def make_params(gamma: float) -> CompassParams:
    return CompassParams(m=ref.M_LEG, m_h=ref.M_HIP, a=ref.A, b=ref.B, gamma=gamma, g=ref.G)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma-max-deg", type=float, default=6.0)
    parser.add_argument("--step-deg", type=float, default=0.05)
    args = parser.parse_args()

    run_dir = new_run_dir("bifurcation_compass")
    rows = []

    # period-1 branch: 公表 slope から上向きに continuation
    y = np.array(ref.SECTION_GUESS)
    gamma = ref.GAMMA_GAIT
    step = np.deg2rad(args.step_deg)
    y2 = None  # period-2 branch の seed
    while gamma <= np.deg2rad(args.gamma_max_deg):
        fp = find_limit_cycle(make_compass(make_params(gamma)), y)
        if not fp.converged:
            print(f"gamma={np.rad2deg(gamma):.3f}deg: period-1 lost")
            break
        lam = float(np.max(np.abs(fp.eigenvalues)))
        rows.append([gamma, 1, fp.y[0], lam])
        print(f"gamma={np.rad2deg(gamma):.3f}deg p1 theta*={fp.y[0]:.5f} max|l|={lam:.4f}")
        if lam > 1.0 and y2 is None:
            y2 = fp.y.copy()  # 不安定化を検知 → ここから period-2 を試す
        y = fp.y
        gamma += step

    # period-2 branch（period-1 不安定化点から）
    if y2 is not None:
        gamma2 = gamma - step
        y = y2 * (1.0 + 1e-3)  # 対称性を破る微小摂動
        while gamma2 <= np.deg2rad(args.gamma_max_deg):
            fp = find_limit_cycle(make_compass(make_params(gamma2)), y, n_strides=2)
            if not fp.converged:
                print(f"gamma={np.rad2deg(gamma2):.3f}deg: period-2 lost")
                break
            lam = float(np.max(np.abs(fp.eigenvalues)))
            rows.append([gamma2, 2, fp.y[0], lam])
            print(f"gamma={np.rad2deg(gamma2):.3f}deg p2 theta*={fp.y[0]:.5f} max|l|={lam:.4f}")
            y = fp.y
            gamma2 += step

    with (run_dir / "bifurcation.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gamma_rad", "period", "theta_star", "max_abs_lambda"])
        writer.writerows(rows)

    data = np.array([[r[0], r[1], r[2], r[3]] for r in rows])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    for period, marker in [(1, "o"), (2, "s")]:
        sel = data[data[:, 1] == period]
        if sel.size:
            stable = sel[sel[:, 3] < 1.0]
            unstable = sel[sel[:, 3] >= 1.0]
            ax1.plot(np.rad2deg(stable[:, 0]), stable[:, 2], marker, ms=3, label=f"p{period} stable")
            if unstable.size:
                ax1.plot(
                    np.rad2deg(unstable[:, 0]), unstable[:, 2], marker, ms=3, mfc="none",
                    label=f"p{period} unstable",
                )
    ax1.set_ylabel("theta* [rad]")
    ax1.legend()
    sel1 = data[data[:, 1] == 1]
    ax2.plot(np.rad2deg(sel1[:, 0]), sel1[:, 3], "-")
    ax2.axhline(1.0, color="k", lw=0.5)
    ax2.set_xlabel("slope [deg]")
    ax2.set_ylabel("max|lambda| (period-1)")
    fig.savefig(run_dir / "bifurcation.png", dpi=150)
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: 実行・観察報告**

Run: `uv run python scripts/bifurcation_compass.py`
報告に含める: period-1 が不安定化する γ（文献の period-doubling 開始値と比較）、
period-2 branch が見つかったか、その γ 範囲。**文献値との一致/乖離は報告のみ**
（ゲートにしない — 分岐点の精密照合は研究観察の領域）。

- [ ] **Step 5: 全テスト + ruff + Commit**

```bash
git add src/crane/search.py tests/test_search_n_strides.py tests/test_compass_continuation.py scripts/bifurcation_compass.py
git commit -m "feat: slope continuation and period-doubling bifurcation diagram"
```

---

### Task 9: viz の model 非依存化 + walk_compass.py デモ

**Files:**
- Modify: `src/crane/viz.py`, `tests/test_viz.py`
- Create: `scripts/walk_compass.py`

- [ ] **Step 1: viz.py を一般化**

1. `link_points_abs(theta_st: float, psi_sw: float, foot_x: float)` を新設
   （絶対角版。実装は現 link_points の中身そのまま、`psi = theta - phi` の行を
   `psi = psi_sw` に）。既存 `link_points(theta, phi, foot_x)` は
   `return link_points_abs(theta, theta - phi, foot_x)` に委譲（後方互換）。
2. `animate_walk(strides, gamma, out, fps=30, angles_of=None)`:
   `angles_of: Callable[[float, float], tuple[float, float]] | None` —
   状態の最初の2成分 (q0, q1) から (stance 絶対角, swing 絶対角) を返す。
   `None` なら simplest 規約 `lambda q0, q1: (q0, q0 - q1)`。
   フレーム生成と foot_x 前進（x_strike 利用）の `link_points` 呼び出しを
   `link_points_abs(*angles_of(q0, q1), foot_x)` に置換。

テスト追加（tests/test_viz.py 末尾）:

```python
def test_link_points_abs_matches_relative_convention():
    """絶対角版と (θ, φ) 版の整合。"""
    from crane.viz import link_points, link_points_abs

    hip1, swing1 = link_points(0.2, 0.5, foot_x=1.0)
    hip2, swing2 = link_points_abs(0.2, 0.2 - 0.5, foot_x=1.0)
    assert np.allclose(hip1, hip2)
    assert np.allclose(swing1, swing2)
```

- [ ] **Step 2: walk_compass.py を書く**

```python
# scripts/walk_compass.py
"""Phase 2 デモ: compass gait のリミットサイクル発見 → 多歩シミュ → アニメ生成。

Usage: uv run python scripts/walk_compass.py [--gamma-deg <published>] [--strides 30] [--perturb 0.005]
"""

import argparse
import json
import sys

import numpy as np

from crane import references_goswami as ref
from crane.models.compass import CompassParams, make_compass
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_walk, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma-deg", type=float, default=np.rad2deg(ref.GAMMA_GAIT))
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.005)
    args = parser.parse_args()

    gamma = np.deg2rad(args.gamma_deg)
    p = CompassParams(m=ref.M_LEG, m_h=ref.M_HIP, a=ref.A, b=ref.B, gamma=gamma, g=ref.G)
    model = make_compass(p)
    run_dir = new_run_dir(f"compass_g{args.gamma_deg:g}deg")

    fp = find_limit_cycle(model, np.array(ref.SECTION_GUESS))
    if not fp.converged:
        for i, (y, r) in enumerate(fp.history):
            print(f"  newton[{i}] y={y} residual={r:.3e}")
        sys.exit("ERROR: limit cycle not found")
    lam = np.abs(fp.eigenvalues)
    print(f"converged: y*={fp.y}  |lambda|={lam}")
    if lam.max() > 1.0:
        print(f"WARNING: cycle is UNSTABLE (max|lambda|={lam.max():.3f}) — walk will fall")

    y0 = fp.y * (1.0 + args.perturb)
    x = model.lift(y0)
    strides = []
    for i in range(args.strides):
        try:
            result = stride(model, x)
        except StrideError as e:
            print(f"stride {i}: FELL ({e})")
            break
        deviation = float(np.linalg.norm(model.project(result.x_end) - fp.y))
        print(f"stride {i}: t={result.t_step:.4f} deviation={deviation:.3e}")
        strides.append(result)
        x = result.x_end

    plot_phase_portrait(strides, run_dir / "phase_portrait.png")
    animate_walk(strides, gamma, run_dir / "walk.mp4", angles_of=lambda q0, q1: (q0, q1))

    meta = {
        "params": {"m": p.m, "m_h": p.m_h, "a": p.a, "b": p.b, "gamma": p.gamma, "g": p.g},
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": lam.tolist(),
        "stable": bool(lam.max() < 1.0),
        "n_strides_completed": len(strides),
        "perturb": args.perturb,
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

（注: meta の `stable` フィールドと不安定警告は issue #2 の解決をここで先取りする形。
walk_simplest.py への同様の追加は issue #2 側で別途。）

- [ ] **Step 3: 実行確認 + 全テスト + ruff + Commit**

Run: `uv run python scripts/walk_compass.py` → mp4/png/meta 生成、deviation 減衰を確認。
`uv run pytest -v` 全 green、ruff format/check。

```bash
git add src/crane/viz.py tests/test_viz.py scripts/walk_compass.py
git commit -m "feat: model-agnostic animation and compass walk demo"
```

---

### Task 10: ドキュメント更新 + 歩容判定ゲート

**Files:**
- Modify: `GOALS.md`, `README.md`

- [ ] **Step 1: GOALS.md に Phase 2 セクションを記入**

「Phase 2: sympy 導出レイヤー + Compass Gait (McGeer 1990)」の節を実績で置き換える
（モデル系統は Goswami 点足に決定済みなので見出しも修正）。チェックリスト:
derive レイヤー（lagrange/impact）、Garcia 退化ゲート一致、Goswami 文献値記録、
compass リミットサイクル発見+文献照合、continuation+分岐図、
最後に `- [ ] **ぽんぽこ殿の歩容判定**: walk.mp4 が「歩行」に見える`（未チェックで残す）。
Phase 2 の知見・実測数値（不動点、固有値、period-doubling 観察結果）を
「Phase 1 で得た知見」と同じ形式で記録。

- [ ] **Step 2: README.md に Phase 2 結果を追記**

Phase 1 結果の下に Phase 2 節: compass の不動点・安定性・分岐図の要約と
`walk_compass.py` / `bifurcation_compass.py` の使い方。

- [ ] **Step 3: Commit**

```bash
git add GOALS.md README.md
git commit -m "docs: record phase 2 results and gait judgment gate"
```

- [ ] **Step 4: 歩容判定ゲート（人間ゲート — STOP）**

`data/runs/<ts>_compass_*/walk.mp4` と分岐図 `bifurcation.png` をぽんぽこ殿に提示。
**判定が出るまで Phase 3 計画に進まない。**

---

## Self-Review Notes

- 設計書 Phase 2 の要件カバー: sympy レイヤー (Task 2-3)、simplest 再導出検証
  (Task 5 — 設計書の「sympy 層で simplest を再導出して Phase 1 と一致確認」を
  compass の退化パラメータで実現、コードパス1本で DRY)、compass リミットサイクル
  (Task 7)、continuation + period-doubling 観察 (Task 8)、歩容判定 (Task 10) — 全カバー
- モデル系統は McGeer rocker foot ではなく Goswami 点足（ぽんぽこ殿決定、
  GOALS.md の Phase 2 見出しを Task 10 で修正）
- 文献数値はプランに書かない（Task 6 の research が取得、provenance 必須）—
  Phase 1 と同じ規律。プラン内の数値は Phase 1 で実測済みの検証値
  (0.2003109, −0.1998325, 0.5891) のみ
- 型一貫性: HybridModel のフィールド名 (dynamics/strike_value/strike_accept/
  impact/lift/project) を Task 1, 4, 9 で統一確認済み。find_limit_cycle(model, y,
  n_strides) の シグネチャは Task 1 → Task 8 で後方互換に拡張
- リスク: (a) compass の衝突 swap 定式化 — Task 4 の不変量テスト + Task 5 の
  Garcia 退化ゲートの二段防御; (b) Goswami の数表欠如 — digitized 値 + 根拠付き
  許容で Phase 1 と同じ対処; (c) sympy lambdify の性能 — lru_cache でビルド1回、
  実行時は純 numpy
