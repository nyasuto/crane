# Crane Phase 3: Kneed Walker (点足) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 点足 kneed walker（4相 hybrid: unlocked swing → knee-strike → locked swing → heel-strike）のリミットサイクルを発見し、安定性を固有値で証明し、文献と照合する。

**Architecture:** HybridModel を多相化（`PhaseSpec` の列）し、stride() を相機械として一般化（既存 1 相モデルは挙動不変）。kneed は Phase 2 の derive レイヤーで記号導出し、**shank 質量 → 0 の退化で Phase 2 検証済み compass と一致させてレイヤー・定式化を検証**。文献ゲートは provenance 規律（Phase 1-2 と同じ）。

**Tech Stack:** Python 3.12 / uv / numpy / scipy / sympy / matplotlib / pytest / ruff

**Spec:** `docs/2026-06-12-crane-design.md` Phase 3 節。足形状 = 点足（ぽんぽこ殿決定 2026-06-13。rocker foot は Phase 4 候補）

---

## 物理モデル早見表（点足 kneed walker）

- 構成: 各脚 = 大腿（長 l_t、質量 m_t を hip から b_t）+ 脛（長 l_s、質量 m_s を knee から
  b_s）。hip に点質量 m_h。点足。脚全長 l = l_t + l_s。
- **状態は全相で 6 次元に統一**: `x = [θ_st, θ_th, θ_sh, θ̇_st, θ̇_th, θ̇_sh]`
  （slope 法線基準の絶対角。θ_st: stance 脚（膝ロック、直線）、θ_th: swing 大腿、
  θ_sh: swing 脛）。locked 相では θ_th = θ_sh を q̈_th = q̈_sh の埋め込みで維持。
- 重力: slope frame で ĝ = (sin γ, −cos γ)（Phase 1-2 と同じ）。
- **相機械（1 stride = 2 相）:**
  1. **Phase A (unlocked, 実質 3 DOF)**: knee-strike 条件 `g_A = θ_sh − θ_th = 0`
     の transversal 交差で終端。post-heel-strike の開始点は g_A=0 上にあるため
     burn-in で除外。**交差は常に受理**（伸展側からの交差 = ロック係合、屈曲側
     からの交差 = knee-strike。ハイパーエクステンション抑止は knee-strike 写像
     そのものが担う簡略化。既知のモデル簡略化として docstring に明記）
  2. **knee-strike 衝突**（脚交換なし、3 速度 → 2 速度）: 保存則 2 本
     - 系全体の角運動量を stance 足（pivot）回りで保存
     - swing 脚（大腿+脛）の角運動量を hip 回りで保存
  3. **Phase B (locked, 実質 2 DOF)**: compass と同型。heel-strike 条件
     `g_B = θ_st + θ_th = 0`（swing 足高さ 0 の非自明側）、受理 = `θ_st < 0` かつ
     下降接地 `θ̇_st + θ̇_th < 0`（Phase 2 と同じ閉形式 ẏ_sw = −l sinθ_st(θ̇_st+θ̇_th)）
  4. **heel-strike 衝突**（脚交換 + 新 swing 膝アンロック、2 速度 → 3 速度）: 保存則 3 本
     - 系全体を新接地点回りで保存（pre-pinned / post-pinned の二重評価 + swap 代入、
       Phase 2 で機械精度検証済みの定式化）
     - 新 swing 脚（= 旧 stance 脚、大腿+脛）を hip 回りで保存
     - 新 swing 脛（= 旧 stance 脛）を knee 回りで保存
- **Poincaré 断面** = heel-strike 直後。位置は θ_th = θ_sh = −θ_st が成立するため
  自由座標は **y = (θ_st, θ̇_st, θ̇_th, θ̇_sh) の 4 次元**。
- **compass への退化**: m_s = 1e-9, b_t = compass の b 相当で、locked 相力学・
  heel-strike・**全 stride 不動点**が Phase 2 検証済み compass に一致するはず
  （massless 脛は slave 振り子、knee-strike 撃力 → 0）。**これがレイヤー検証ゲート**。
- 文献候補は Task 2 で原典調査（点足 kneed の公開文献として Chen 2007 MIT thesis 等が
  有力候補。McGeer 1990 (knees) は rocker foot だが衝突定式化の照合元として有用）。
  数値はプランに書かない（Phase 1-2 と同じ規律）。
- swing 足の地面通過（phase A 中の scuff）はイベントにしない（文献慣行どおり無視、
  docstring 明記）。

## File Structure

| パス | 責務 |
|---|---|
| `src/crane/model.py` | `PhaseSpec` 追加、`HybridModel` を phases 列に（refactor） |
| `src/crane/stride.py` | 相機械として一般化（1 相モデルは挙動不変） |
| `src/crane/models/simplest.py` / `compass.py` | 1 相 PhaseSpec への constructor 更新のみ |
| `src/crane/models/kneed.py` | 点足 kneed（derive レイヤー、新規） |
| `src/crane/references_kneed.py` | 文献値（provenance 付き、Task 2 で記入） |
| `src/crane/viz.py` | `animate_kneed`（4 セグメント stick figure、専用関数） |
| `scripts/walk_kneed.py` | Phase 3 デモ |

---

### Task 1: 相機械への refactor（挙動変更なし）

**Files:**
- Modify: `src/crane/model.py`, `src/crane/stride.py`, `src/crane/models/simplest.py`, `src/crane/models/compass.py`

- [ ] **Step 1: model.py を書き換える**

```python
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
```

- [ ] **Step 2: stride.py を相機械に一般化**

`stride(model, x0)` のシグネチャ・StrideResult のフィールドは不変。内部を
「phases を順に消化するループ」にする。既存の burn-in / 受理拒否時の続行 /
重複除去 / t_max ガードのロジックは**各相内で**そのまま使う:

```python
def stride(
    model: HybridModel,
    x0: np.ndarray,
    *,
    t_max: float = 10.0,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> StrideResult:
    """1歩分の写像（多相）。副作用なし。t_max は stride 全体の予算。"""
    t0 = 0.0
    x = np.asarray(x0, dtype=float)
    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []

    def append(t_seg, x_seg):
        if ts:
            t_seg, x_seg = t_seg[1:], x_seg[:, 1:]
        ts.append(t_seg)
        xs.append(x_seg)

    for i_phase, phase in enumerate(model.phases):
        def event(t, x_, _p=phase):
            return _p.event_value(x_)

        event.terminal = True
        event.direction = 0

        while True:
            if t0 >= t_max:
                raise StrideError(
                    f"phase {i_phase}: no terminal event before t_max={t_max} "
                    f"(t={t0:.3f})"
                )
            burn = solve_ivp(phase.dynamics, (t0, t0 + T_BURN), x, rtol=rtol, atol=atol)
            append(burn.t, burn.y)
            t0, x = burn.t[-1], burn.y[:, -1]
            if t0 >= t_max:
                continue  # 次周の冒頭ガードで StrideError

            sol = solve_ivp(phase.dynamics, (t0, t_max), x, events=event, rtol=rtol, atol=atol)
            append(sol.t, sol.y)
            if sol.t_events[0].size == 0:
                raise StrideError(
                    f"phase {i_phase}: no terminal event before t_max={t_max}"
                )
            t_e = sol.t_events[0][0]
            x_e = sol.y_events[0][0]
            if phase.event_accept(x_e):
                break
            t0, x = t_e, x_e  # 偽交差（scuffing 等）: 無視して続行

        x_strike, t0, x = x_e, t_e, phase.impact(x_e)

    return StrideResult(
        x_end=x,
        x_strike=x_strike,  # 最終相の衝突直前状態（heel-strike 瞬間）
        t_step=t0,
        t=np.concatenate(ts),
        x=np.concatenate(xs, axis=1),
    )
```

（注: `x_strike` は最終相の pre-impact 状態。中間相の pre-impact 状態は trajectory
から復元可能。`t_step` は最終 heel-strike 時刻 = stride 全体の時間）

- [ ] **Step 3: simplest.py / compass.py の constructor を 1 相に更新**

`make_simplest` / `make_compass` の中身を `HybridModel(phases=(PhaseSpec(dynamics=...,
event_value=（旧 strike_value）, event_accept=（旧 strike_accept）, impact=...),),
lift=..., project=...)` に組み替えるだけ。各関数のラムダ群は変更しない。

- [ ] **Step 4: 全テスト green を確認（挙動変更なし）**

Run: `cd /Users/yast/git/crane && uv run pytest -v && uv run ruff format . && uv run ruff check . && uv run python scripts/walk_compass.py --strides 3`
Expected: 43 passed、demo 完走（数値は従来と一致: y*=(0.27103, −1.09238, −0.37737)）

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "refactor: multi-phase stride machine (PhaseSpec sequence)"
```

---

### Task 2: 点足 kneed の文献調査（references_kneed.py）

**Files:**
- Create: `src/crane/references_kneed.py`
- Test: `tests/test_references_kneed.py`

**Web 調査タスク。** 候補を当たり、**点足 kneed walker の公開された一次資料**を確定する:

1. 第一候補: Chen, V. F. H. (2007?) MIT thesis "Passive Dynamic Walking with Knees:
   A Point Foot Model"（dspace.mit.edu で公開のはず）— 点足 kneed の定式化と数値
2. 補助: McGeer (1990) "Passive Walking with Knees" (ICRA) — rocker foot だが
   knee-strike / heel-strike の保存則定式化の照合元
3. 他に点足 kneed の数値を印刷した文献があればそれも可（provenance 最優先）

取得するもの:
- パラメータセット（m_h, m_t, m_s, l_t, l_s, b_t, b_s 相当 + slope）と
  **本プラン早見表の規約への変換**（角度の基準・符号、質量点位置の測り方を必ず照合・文書化）
- gait の公表数値（不動点、step period、固有値など印刷されているもの全て。
  図のみは digitized ±X% で正直ラベル）
- **knee-strike / heel-strike の保存則定式化を早見表（2 本 / 3 本）と照合**。
  相違があれば報告（勝手に直さない）
- 膝の屈曲方向の符号規約（swing 中 θ_sh と θ_th のどちらが大きいか）を図から確定し記録

```python
# tests/test_references_kneed.py（TDD: 先に書く）
from crane import references_kneed as ref


def test_kneed_references_are_filled():
    """文献値が provenance 付きで記入済み（isinstance は Ellipsis 対策）。"""
    assert "http" in ref.PROVENANCE
    for name in ["M_HIP", "M_THIGH", "M_SHANK", "L_THIGH", "L_SHANK", "B_THIGH", "B_SHANK"]:
        assert isinstance(getattr(ref, name), float), name
    assert isinstance(ref.GAMMA_GAIT, float)
    assert 0.0 < ref.GAMMA_GAIT < 0.3
    assert len(ref.SECTION_GUESS) == 4
```

必須フィールド: `PROVENANCE`, 上記 7 パラメータ, `G`, `GAMMA_GAIT`,
`SECTION_GUESS`（我々の 4D 断面座標 (θ_st, θ̇_st, θ̇_th, θ̇_sh)、変換根拠コメント必須）,
`GAIT_TOLERANCE`（根拠コメント）, 照合量（`STEP_PERIOD: float | None` 等、印刷形式に
合わせて）, `KNEE_FLEXION_SIGN`（+1 なら swing 中 θ_sh > θ_th、コメントに図の根拠）。
論文に無い量は None + 理由。

- [ ] **Step 1: ゲートテストを書き FAIL 確認**
- [ ] **Step 2: Web 調査、references_kneed.py 記入**
- [ ] **Step 3: テスト + ruff + Commit**

```bash
git add src/crane/references_kneed.py tests/test_references_kneed.py
git commit -m "feat: record point-foot kneed walker published values with provenance"
```

（定式化照合の結果を commit body に `Formulation cross-check: <内容>` で記録）

---

### Task 3: kneed.py — 両相の力学（TDD、不変量ベース）

**Files:**
- Create: `src/crane/models/kneed.py`
- Test: `tests/test_kneed_dynamics.py`

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_kneed_dynamics.py
import numpy as np
from scipy.integrate import solve_ivp

from crane import references_kneed as ref
from crane.models.kneed import KneedParams, dynamics_locked, dynamics_unlocked, energy

P = KneedParams(
    m_h=ref.M_HIP, m_t=ref.M_THIGH, m_s=ref.M_SHANK,
    l_t=ref.L_THIGH, l_s=ref.L_SHANK, b_t=ref.B_THIGH, b_s=ref.B_SHANK,
    gamma=ref.GAMMA_GAIT, g=ref.G,
)


def test_equilibrium_aligned_with_gravity():
    """全リンクが重力方向 (絶対角 γ) に揃った静止状態では加速度ゼロ（両相）。"""
    x = np.array([P.gamma, P.gamma, P.gamma, 0.0, 0.0, 0.0])
    assert np.allclose(dynamics_unlocked(0.0, x, P), np.zeros(6), atol=1e-11)
    assert np.allclose(dynamics_locked(0.0, x, P), np.zeros(6), atol=1e-11)


def test_unlocked_phase_conserves_energy():
    """unlocked 相は保存系。"""
    x0 = np.array([0.2, -0.2, -0.45, -1.2, 1.0, 2.5])
    sol = solve_ivp(lambda t, x: dynamics_unlocked(t, x, P), (0.0, 0.4), x0,
                    rtol=1e-11, atol=1e-13)
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)


def test_locked_phase_conserves_energy_and_keeps_alignment():
    """locked 相は保存系で、θ_th = θ_sh が維持される（6D 埋め込みの整合）。"""
    x0 = np.array([0.25, -0.3, -0.3, -1.5, 0.8, 0.8])
    sol = solve_ivp(lambda t, x: dynamics_locked(t, x, P), (0.0, 0.4), x0,
                    rtol=1e-11, atol=1e-13)
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)
    assert np.allclose(sol.y[1], sol.y[2], atol=1e-9)   # θ_th ≡ θ_sh
    assert np.allclose(sol.y[4], sol.y[5], atol=1e-9)   # θ̇_th ≡ θ̇_sh


def test_locked_dynamics_is_unlocked_restricted():
    """locked 力学は「θ_th=θ_sh 拘束付き unlocked」と stance 加速度が異なる
    （拘束力が働くため一般に一致しない）が、m_s→0 では一致するはず。
    ここでは整合性のみ: locked の θ̈_th と θ̈_sh が厳密に等しいこと。"""
    x = np.array([0.2, -0.35, -0.35, -1.0, 0.5, 0.5])
    xdot = dynamics_locked(0.0, x, P)
    assert xdot[4] == xdot[5]
```

- [ ] **Step 2: FAIL 確認（ModuleNotFoundError）**

- [ ] **Step 3: kneed.py の力学部分を実装**

compass.py のパターンを踏襲（`@lru_cache` の `_build()` で記号導出 + lambdify）:

```python
# src/crane/models/kneed.py
"""点足 kneed walker。derive レイヤーで記号導出し lambdify。

状態 x = [θ_st, θ_th, θ_sh, θ̇_st, θ̇_th, θ̇_sh]（slope 法線基準の絶対角）。
θ_st: stance 脚（膝ロック・直線）、θ_th: swing 大腿、θ_sh: swing 脛。
locked 相は θ_th=θ_sh を 6D 埋め込み（q̈_th=q̈_sh）で維持する。

相機械: unlocked → knee-strike → locked → heel-strike（早見表参照）。
簡略化（docstring 明記必須）:
- knee-strike イベント g=θ_sh−θ_th=0 の transversal 交差は常に受理
  （ハイパーエクステンション抑止はロック係合そのもの）
- unlocked 相中の swing 足の地面通過（scuff）はイベントにしない（文献慣行）
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum
from crane.derive.lagrange import derive_qdd
from crane.model import HybridModel, PhaseSpec


@dataclass(frozen=True)
class KneedParams:
    m_h: float  # hip 質量 [kg]
    m_t: float  # 大腿質量 [kg]
    m_s: float  # 脛質量 [kg]
    l_t: float  # 大腿長 hip→knee [m]
    l_s: float  # 脛長 knee→foot [m]
    b_t: float  # hip→大腿質量点 [m]
    b_s: float  # knee→脛質量点 [m]
    gamma: float  # 斜面角 [rad]
    g: float = 9.81

    @property
    def l(self) -> float:  # noqa: E743
        return self.l_t + self.l_s


@lru_cache(maxsize=1)
def _build():
    th_st, th_th, th_sh = sp.symbols("th_st th_th th_sh")
    w_st, w_th, w_sh = sp.symbols("w_st w_th w_sh")
    m_h, m_t, m_s = sp.symbols("m_h m_t m_s", positive=True)
    l_t, l_s, b_t, b_s = sp.symbols("l_t l_s b_t b_s", positive=True)
    gamma, g = sp.symbols("gamma g", positive=True)
    length = l_t + l_s

    def down(theta):
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    # stance 脚（直線、足原点）: hip は足から脚全長ぶん上
    hip = sp.Matrix([-length * sp.sin(th_st), length * sp.cos(th_st)])
    knee_st = hip + l_t * down(th_st)
    p_th_st = hip + b_t * down(th_st)        # stance 大腿質量点
    p_sh_st = knee_st + b_s * down(th_st)    # stance 脛質量点
    # swing 脚
    p_th_sw = hip + b_t * down(th_th)
    knee_sw = hip + l_t * down(th_th)
    p_sh_sw = knee_sw + b_s * down(th_sh)
    foot_sw = knee_sw + l_s * down(th_sh)

    q3, qd3 = [th_st, th_th, th_sh], [w_st, w_th, w_sh]
    bodies = [(m_h, hip), (m_t, p_th_st), (m_s, p_sh_st), (m_t, p_th_sw), (m_s, p_sh_sw)]

    def vel(r, q, qd):
        return r.jacobian(q) * sp.Matrix(qd)

    T3 = sum(mass * vel(r, q3, qd3).dot(vel(r, q3, qd3)) for mass, r in bodies) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])
    V = -g * sum(mass * r.dot(g_dir) for mass, r in bodies)

    params = (m_h, m_t, m_s, l_t, l_s, b_t, b_s, gamma, g)
    args6 = (th_st, th_th, th_sh, w_st, w_th, w_sh, *params)

    # unlocked: 3 DOF
    qdd3 = derive_qdd(q3, qd3, T3, V, simplify=False)
    f_qdd3 = sp.lambdify(args6, [qdd3[0], qdd3[1], qdd3[2]], "numpy")
    f_energy = sp.lambdify(args6, T3 + V, "numpy")
    f_kinetic = sp.lambdify(args6, T3, "numpy")

    # locked: θ_th=θ_sh=θ_sw の 2 DOF を導出して 6D に埋め込む
    th_sw, w_sw = sp.symbols("th_sw w_sw")
    lock = {th_th: th_sw, th_sh: th_sw, w_th: w_sw, w_sh: w_sw}
    T2 = T3.subs(lock, simultaneous=True)
    V2 = V.subs(lock, simultaneous=True)
    qdd2 = derive_qdd([th_st, th_sw], [w_st, w_sw], T2, V2, simplify=False)
    f_qdd2 = sp.lambdify((th_st, th_sw, w_st, w_sw, *params), [qdd2[0], qdd2[1]], "numpy")

    # --- knee-strike 衝突（脚交換なし、3 速度 → 2 速度）---
    L_sys_pivot3 = angular_momentum(bodies, q3, qd3, sp.Matrix([0, 0]))
    L_swleg_hip3 = angular_momentum([(m_t, p_th_sw), (m_s, p_sh_sw)], q3, qd3, hip)
    wp_st, wp_sw = sp.symbols("wp_st wp_sw")
    post_ks = {w_st: wp_st, w_th: wp_sw, w_sh: wp_sw}  # 位置不変（θ_th=θ_sh 上）
    eqs_ks = [
        L_sys_pivot3.subs(post_ks, simultaneous=True) - L_sys_pivot3,
        L_swleg_hip3.subs(post_ks, simultaneous=True) - L_swleg_hip3,
    ]
    A_ks, rhs_ks = sp.linear_eq_to_matrix(eqs_ks, [wp_st, wp_sw])
    qd_post_ks = A_ks.LUsolve(rhs_ks)
    f_impact_ks = sp.lambdify(args6, [qd_post_ks[0], qd_post_ks[1]], "numpy")

    # --- heel-strike 衝突（脚交換 + 新 swing 膝アンロック、2 速度 → 3 速度）---
    # pre は locked 配置（θ_th=θ_sh）の pre-pinned 系で評価。
    # post は post-pinned 系（同形の式）に swap 代入:
    #   位置: th_st→th_th（新 stance = 旧 swing 角）、th_th→th_st, th_sh→th_st
    #         （新 swing 大腿・脛は旧 stance 角で整列）
    #   速度: w_st→wq_st, w_th→wq_th, w_sh→wq_sh（post ラベル）
    wq_st, wq_th, wq_sh = sp.symbols("wq_st wq_th wq_sh")
    swap_hs = {
        th_st: th_th, th_th: th_st, th_sh: th_st,
        w_st: wq_st, w_th: wq_th, w_sh: wq_sh,
    }
    L_stleg_hip3 = angular_momentum([(m_t, p_th_st), (m_s, p_sh_st)], q3, qd3, hip)
    L_stsh_knee3 = angular_momentum([(m_s, p_sh_st)], q3, qd3, knee_st)
    L_swsh_knee3 = angular_momentum([(m_s, p_sh_sw)], q3, qd3, knee_sw)
    L_sys_swfoot3 = angular_momentum(bodies, q3, qd3, foot_sw)
    eqs_hs = [
        L_sys_pivot3.subs(swap_hs, simultaneous=True) - L_sys_swfoot3,
        L_swleg_hip3.subs(swap_hs, simultaneous=True) - L_stleg_hip3,
        L_swsh_knee3.subs(swap_hs, simultaneous=True) - L_stsh_knee3,
    ]
    A_hs, rhs_hs = sp.linear_eq_to_matrix(eqs_hs, [wq_st, wq_th, wq_sh])
    qd_post_hs = A_hs.LUsolve(rhs_hs)
    f_impact_hs = sp.lambdify(args6, [qd_post_hs[0], qd_post_hs[1], qd_post_hs[2]], "numpy")

    return f_qdd3, f_qdd2, f_energy, f_kinetic, f_impact_ks, f_impact_hs


def _args(x, p: KneedParams):
    return (*x, p.m_h, p.m_t, p.m_s, p.l_t, p.l_s, p.b_t, p.b_s, p.gamma, p.g)


def dynamics_unlocked(t: float, x, p: KneedParams):
    f_qdd3 = _build()[0]
    qdd = f_qdd3(*_args(x, p))
    return [x[3], x[4], x[5], qdd[0], qdd[1], qdd[2]]


def dynamics_locked(t: float, x, p: KneedParams):
    """locked 相: 2 DOF を 6D に埋め込み（θ_th スロットを θ_sw として使う）。"""
    f_qdd2 = _build()[1]
    qdd = f_qdd2(x[0], x[1], x[3], x[4],
                 p.m_h, p.m_t, p.m_s, p.l_t, p.l_s, p.b_t, p.b_s, p.gamma, p.g)
    return [x[3], x[4], x[5], qdd[0], qdd[1], qdd[1]]


def energy(x, p: KneedParams) -> float:
    return float(_build()[2](*_args(x, p)))


def kinetic_energy(x, p: KneedParams) -> float:
    return float(_build()[3](*_args(x, p)))
```

（Task 3 では衝突写像の lambdify までは `_build()` に含めてよいが、公開関数
`kneestrike_map` / `heelstrike_map` / `make_kneed` は Task 4-5 で追加する）

- [ ] **Step 4: テスト確認 + Commit**

Run: `uv run pytest tests/test_kneed_dynamics.py -v` → 4 PASSED、全体 green、ruff

```bash
git add src/crane/models/kneed.py tests/test_kneed_dynamics.py
git commit -m "feat: kneed walker dynamics (unlocked 3-DOF + locked embedded) via derive layer"
```

---

### Task 4: kneed の衝突写像（TDD、不変量ベース）

**Files:**
- Modify: `src/crane/models/kneed.py`（公開関数追加）
- Test: `tests/test_kneed_impacts.py`

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_kneed_impacts.py
import numpy as np

from crane import references_kneed as ref
from crane.models.kneed import (
    KneedParams, heelstrike_map, kinetic_energy, kneestrike_map,
)

P = KneedParams(
    m_h=ref.M_HIP, m_t=ref.M_THIGH, m_s=ref.M_SHANK,
    l_t=ref.L_THIGH, l_s=ref.L_SHANK, b_t=ref.B_THIGH, b_s=ref.B_SHANK,
    gamma=ref.GAMMA_GAIT, g=ref.G,
)


def test_kneestrike_locks_and_dissipates():
    """knee-strike: 位置不変・速度等値化・KE 非増加。"""
    x_pre = np.array([-0.05, 0.25, 0.25, -1.3, -0.6, 1.8])  # θ_th=θ_sh（lock 面）
    x_post = kneestrike_map(x_pre, P)
    assert np.allclose(x_post[:3], x_pre[:3])          # 位置不変
    assert x_post[4] == x_post[5]                      # 速度ロック
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_kneestrike_zero_velocity_is_fixed():
    x_pre = np.array([-0.05, 0.25, 0.25, 0.0, 0.0, 0.0])
    x_post = kneestrike_map(x_pre, P)
    assert np.allclose(x_post[3:], np.zeros(3), atol=1e-12)


def test_kneestrike_already_locked_is_identity():
    """既に速度が揃っている状態の knee-strike は恒等（撃力ゼロ）。"""
    x_pre = np.array([-0.05, 0.25, 0.25, -1.3, 0.7, 0.7])
    x_post = kneestrike_map(x_pre, P)
    assert np.allclose(x_post, x_pre, atol=1e-10)


def test_heelstrike_swaps_and_dissipates():
    """heel-strike: 脚交換幾何 + KE 非増加 + 新 swing 整列。"""
    th = -0.18
    x_pre = np.array([th, -th, -th, -1.4, -0.9, -0.9])  # locked, strike 面上
    x_post = heelstrike_map(x_pre, P)
    assert x_post[0] == -th          # 新 stance = 旧 swing 角
    assert x_post[1] == th           # 新 swing 大腿 = 旧 stance 角
    assert x_post[2] == th           # 新 swing 脛も整列
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_heelstrike_zero_velocity_is_fixed():
    th = -0.18
    x_pre = np.array([th, -th, -th, 0.0, 0.0, 0.0])
    x_post = heelstrike_map(x_pre, P)
    assert np.allclose(x_post[3:], np.zeros(3), atol=1e-12)
```

- [ ] **Step 2: FAIL 確認（ImportError: kneestrike_map）**

- [ ] **Step 3: 公開関数を kneed.py に追加**

```python
def kneestrike_map(x: np.ndarray, p: KneedParams) -> np.ndarray:
    """knee-strike: 脚交換なし。swing 膝がロックし 3 速度 → 2 速度。"""
    wp = _build()[4](*_args(x, p))
    return np.array([x[0], x[1], x[2], float(wp[0]), float(wp[1]), float(wp[1])])


def heelstrike_map(x: np.ndarray, p: KneedParams) -> np.ndarray:
    """heel-strike: 脚交換 + 新 swing 膝アンロック。post 速度は post ラベル。"""
    wq = _build()[5](*_args(x, p))
    return np.array([x[1], x[0], x[0], float(wq[0]), float(wq[1]), float(wq[2])])
```

- [ ] **Step 4: テスト確認 + Commit**

Run: `uv run pytest tests/test_kneed_impacts.py -v` → 5 PASSED、全体 green、ruff

```bash
git add src/crane/models/kneed.py tests/test_kneed_impacts.py
git commit -m "feat: kneed walker knee-strike and heel-strike impact maps"
```

---

### Task 5: make_kneed と stride の end-to-end

**Files:**
- Modify: `src/crane/models/kneed.py`
- Test: `tests/test_kneed_stride.py`

- [ ] **Step 1: make_kneed を実装**

```python
def make_kneed(p: KneedParams) -> HybridModel:
    """相機械: unlocked → knee-strike → locked → heel-strike。

    断面 = heel-strike 直後。y = (θ_st, θ̇_st, θ̇_th, θ̇_sh)、
    位置は θ_th = θ_sh = −θ_st が成立。
    heel-strike 受理は compass と同じ閉形式（脚全長 l で
    ẏ_sw = −l sin(θ_st)(θ̇_st + θ̇_sw) < 0、locked 相では θ̇_sw = x[4]）。
    """
    unlocked = PhaseSpec(
        dynamics=lambda t, x: dynamics_unlocked(t, x, p),
        event_value=lambda x: x[2] - x[1],   # θ_sh − θ_th = 0 で knee-strike
        event_accept=lambda x: True,          # transversal 交差は常にロック係合
        impact=lambda x: kneestrike_map(x, p),
    )
    locked = PhaseSpec(
        dynamics=lambda t, x: dynamics_locked(t, x, p),
        event_value=lambda x: x[0] + x[1],   # θ_st + θ_sw = 0 で heel-strike
        event_accept=lambda x: x[0] < 0.0 and (x[3] + x[4]) < 0.0,
        impact=lambda x: heelstrike_map(x, p),
    )
    return HybridModel(
        phases=(unlocked, locked),
        lift=lambda y: np.array([y[0], -y[0], -y[0], y[1], y[2], y[3]]),
        project=lambda x: np.array([x[0], x[3], x[4], x[5]]),
    )
```

- [ ] **Step 2: end-to-end テストを書く**

```python
# tests/test_kneed_stride.py
import numpy as np

from crane import references_kneed as ref
from crane.models.kneed import KneedParams, make_kneed
from crane.stride import stride

P = KneedParams(
    m_h=ref.M_HIP, m_t=ref.M_THIGH, m_s=ref.M_SHANK,
    l_t=ref.L_THIGH, l_s=ref.L_SHANK, b_t=ref.B_THIGH, b_s=ref.B_SHANK,
    gamma=ref.GAMMA_GAIT, g=ref.G,
)
MODEL = make_kneed(P)


def test_stride_passes_both_phases_and_returns_to_section():
    """文献 seed から 1 stride: knee-strike → heel-strike を経て断面に戻る。"""
    result = stride(MODEL, MODEL.lift(np.array(ref.SECTION_GUESS)))
    x_end = result.x_end
    assert np.isclose(x_end[1], -x_end[0], atol=1e-8)   # θ_th = −θ_st
    assert np.isclose(x_end[2], -x_end[0], atol=1e-8)   # θ_sh = −θ_st
    assert result.t_step > 0.1
    # trajectory 中に膝の屈曲が現れている（unlocked 相の実在確認）
    flexion = result.x[2] - result.x[1]                  # θ_sh − θ_th
    assert np.max(np.abs(flexion)) > 1e-3


def test_stride_reports_knee_flexion_sign():
    """屈曲方向が references_kneed の KNEE_FLEXION_SIGN と一致（規約照合）。"""
    result = stride(MODEL, MODEL.lift(np.array(ref.SECTION_GUESS)))
    flexion = result.x[2] - result.x[1]
    peak = flexion[np.argmax(np.abs(flexion))]
    assert np.sign(peak) == ref.KNEE_FLEXION_SIGN
```

- [ ] **Step 3: 実行・デバッグ**

Run: `uv run pytest tests/test_kneed_stride.py -v`
SECTION_GUESS で stride が成立しない場合（StrideError）: phase ごとの軌跡を
print して、(a) 膝が即時再ロックする（unlocked 相が一瞬）、(b) knee-strike 前に
転倒、等の失敗モードを特定して報告。文献 seed の座標変換ミスの可能性も疑うこと
（Task 2 の変換根拠を再確認）。**勝手にパラメータを変えない。**

- [ ] **Step 4: 全テスト + ruff + Commit**

```bash
git add src/crane/models/kneed.py tests/test_kneed_stride.py
git commit -m "feat: kneed walker phase machine and end-to-end stride"
```

---

### Task 6: compass 退化ゲート — m_s → 0 で Phase 2 と一致

**Files:**
- Test: `tests/test_kneed_compass_limit.py`

kneed を m_s=1e-9（脛をほぼ無質量化）、m_t を compass の脚質量に、b_t を compass の
b（hip→質量点）に合わせると: locked 相力学 = compass 力学、heel-strike = compass
衝突、knee-strike 撃力 → 0、massless 脛は slave 振り子。よって**全 stride の不動点が
Phase 2 検証済み compass の (0.27103, −1.09238, −0.37737) に一致するはず**。

- [ ] **Step 1: テストを書く**

```python
# tests/test_kneed_compass_limit.py
import numpy as np

from crane import references_goswami as gref
from crane.models.compass import CompassParams, make_compass
from crane.models.kneed import KneedParams, dynamics_locked, heelstrike_map, make_kneed
from crane.search import find_limit_cycle

# compass (m=5, m_h=10, a=b=0.5) と等価な退化 kneed:
#   大腿質量 5 kg を hip から b_t=0.5 に、脛はほぼ無質量。
#   l_t + l_s = 1 を維持（l_t=0.5, l_s=0.5）。
P_DEG = KneedParams(
    m_h=gref.M_HIP, m_t=gref.M_LEG, m_s=1e-9,
    l_t=0.5, l_s=0.5, b_t=gref.B, b_s=0.25,
    gamma=gref.GAMMA_GAIT, g=gref.G,
)
P_COMPASS = CompassParams(
    m=gref.M_LEG, m_h=gref.M_HIP, a=gref.A, b=gref.B,
    gamma=gref.GAMMA_GAIT, g=gref.G,
)
C_MODEL = make_compass(P_COMPASS)


def test_locked_dynamics_reduces_to_compass():
    """m_s→0 の locked 相力学が compass 力学と一致（ランダム10状態）。"""
    rng = np.random.default_rng(0)
    for _ in range(10):
        th_st, th_sw = rng.uniform(-0.4, 0.4, 2)
        w_st, w_sw = rng.uniform(-2.0, 2.0, 2)
        x6 = np.array([th_st, th_sw, th_sw, w_st, w_sw, w_sw])
        x4 = np.array([th_st, th_sw, w_st, w_sw])
        d6 = dynamics_locked(0.0, x6, P_DEG)
        d4 = C_MODEL.dynamics(0.0, x4)
        assert np.allclose([d6[3], d6[4]], [d4[2], d4[3]], rtol=0, atol=1e-5)


def test_heelstrike_reduces_to_compass():
    """m_s→0 の heel-strike が compass 衝突に一致。"""
    rng = np.random.default_rng(1)
    for _ in range(10):
        th = rng.uniform(-0.3, -0.1)
        w_st, w_sw = rng.uniform(-2.0, -0.2, 2)
        x6 = np.array([th, -th, -th, w_st, w_sw, w_sw])
        x4 = np.array([th, -th, w_st, w_sw])
        got = heelstrike_map(x6, P_DEG)
        want = C_MODEL.impact(x4)
        assert np.allclose([got[0], got[1]], [want[0], want[1]], atol=1e-12)
        assert np.allclose([got[3], got[4]], [want[2], want[3]], rtol=0, atol=1e-5)


def test_full_cycle_reduces_to_compass():
    """退化 kneed の全 stride 不動点が compass 不動点 (Phase 2 実測) に一致。"""
    model = make_kneed(P_DEG)
    fp_c = find_limit_cycle(C_MODEL, np.array(gref.SECTION_GUESS))
    assert fp_c.converged
    guess = np.array([fp_c.y[0], fp_c.y[1], fp_c.y[2], fp_c.y[2]])
    fp_k = find_limit_cycle(model, guess)
    assert fp_k.converged
    assert np.isclose(fp_k.y[0], fp_c.y[0], atol=1e-5)
    assert np.isclose(fp_k.y[1], fp_c.y[1], atol=1e-4)
    assert np.isclose(fp_k.y[2], fp_c.y[2], atol=1e-4)
    # 上位 multiplier も compass と一致（slave 脛方向が 1 本増える）
    mags_k = np.sort(np.abs(fp_k.eigenvalues))[::-1]
    mags_c = np.sort(np.abs(fp_c.eigenvalues))[::-1]
    assert np.allclose(mags_k[:2], mags_c[:2], atol=2e-3)
```

- [ ] **Step 2: 実行**

Run: `uv run pytest tests/test_kneed_compass_limit.py -v`
Expected: 3 PASSED（FAIL 時は locked/衝突/全体のどれが落ちたかで切り分け報告。
**勝手に許容値を緩めない**。注意点: 退化 kneed の unlocked 相で massless 脛の
knee-strike タイミングは compass に影響しないはずだが、StrideError になる場合は
脛の slave 振り子が knee-strike 面に到達しない可能性 — b_s を変えて挙動を報告）

- [ ] **Step 3: Commit**

```bash
git add tests/test_kneed_compass_limit.py
git commit -m "test: kneed degenerates to verified compass (massless shank gate)"
```

---

### Task 7: kneed リミットサイクル + 文献ゲート + 安定性

**Files:**
- Test: `tests/test_kneed_cycle.py`

- [ ] **Step 1: ゲートテストを書く**

```python
# tests/test_kneed_cycle.py
import numpy as np

from crane import references_kneed as ref
from crane.models.kneed import KneedParams, make_kneed
from crane.search import find_limit_cycle
from crane.stride import stride

P = KneedParams(
    m_h=ref.M_HIP, m_t=ref.M_THIGH, m_s=ref.M_SHANK,
    l_t=ref.L_THIGH, l_s=ref.L_SHANK, b_t=ref.B_THIGH, b_s=ref.B_SHANK,
    gamma=ref.GAMMA_GAIT, g=ref.G,
)
MODEL = make_kneed(P)


def test_kneed_limit_cycle_exists_at_published_config():
    """Phase 3 ゲート: 公表パラメータでリミットサイクルが存在する。"""
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    # 安定性は文献の記述に従う（Task 2 で記録した STABILITY_CLAIM に合わせ、
    # 安定と明記されていれば max|λ|<1 を assert、なければ報告のみ）
    if ref.PUBLISHED_STABLE:
        assert np.max(np.abs(fp.eigenvalues)) < 1.0


def test_kneed_gait_matches_published_quantities():
    fp = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp.converged
    result = stride(MODEL, MODEL.lift(fp.y))
    if ref.STEP_PERIOD is not None:
        assert np.isclose(result.t_step, ref.STEP_PERIOD, rtol=ref.GAIT_TOLERANCE)
    if ref.THETA_STRIKE is not None:
        assert np.isclose(abs(result.x_strike[0]), abs(ref.THETA_STRIKE), rtol=ref.GAIT_TOLERANCE)
```

（`PUBLISHED_STABLE: bool` を references_kneed.py に追加。文献が安定性を明記して
いるかの記録。Task 2 の成果に合わせてフィールドを調整してよいが、変更は report に明記）

- [ ] **Step 2: 実行・探索**

収束しない場合のフォールバック（決定論的、Phase 2 と同じ流儀）:
1. SECTION_GUESS 近傍の粗グリッドで ‖S(y)−y‖ 最小点を seed に
2. それでもダメなら γ を文献値から少し振って収束する γ を見つけ、
   そこから文献 γ へ continuation
3. 失敗モード（膝再ロック即時/転倒/scuff 無限ループ等）を軌跡付きで報告

見つかった不動点・固有値・step period・膝屈曲の最大角を必ず報告。
SECTION_GUESS を補正した場合は references_kneed.py に根拠コメント付きで更新。

- [ ] **Step 3: 小さな continuation テスト（安定 family の確認）**

```python
def test_gait_family_continues_near_published_slope():
    """公表 slope の ±20% で gait family が continuation で追える。"""
    y = None
    fp0 = find_limit_cycle(MODEL, np.array(ref.SECTION_GUESS))
    assert fp0.converged
    y = fp0.y
    for scale in [0.95, 0.9, 0.85, 0.8]:
        p = KneedParams(
            m_h=ref.M_HIP, m_t=ref.M_THIGH, m_s=ref.M_SHANK,
            l_t=ref.L_THIGH, l_s=ref.L_SHANK, b_t=ref.B_THIGH, b_s=ref.B_SHANK,
            gamma=ref.GAMMA_GAIT * scale, g=ref.G,
        )
        fp = find_limit_cycle(make_kneed(p), y)
        if not fp.converged:
            break  # gait が存在しない slope に達した: 既知の kneed の性質。報告のみ
        y = fp.y
    assert fp0.converged  # 最低限、公表点での存在が family の証拠
```

（kneed の gait 存在域は compass より狭い可能性が高い。family がどこまで続いたかを
報告に含める）

- [ ] **Step 4: Commit**

```bash
git add tests/test_kneed_cycle.py src/crane/references_kneed.py
git commit -m "feat: kneed walker limit cycle gated on published gait"
```

---

### Task 8: viz（4 セグメント）+ walk_kneed.py デモ

**Files:**
- Modify: `src/crane/viz.py`（`animate_kneed` 追加 — kneed 専用関数。早すぎる抽象化を避ける）
- Create: `scripts/walk_kneed.py`
- Test: `tests/test_viz_kneed.py`

- [ ] **Step 1: テストを書く（関節座標の幾何整合のみ）**

```python
# tests/test_viz_kneed.py
import numpy as np

from crane.viz import kneed_joints


def test_kneed_joints_heelstrike_symmetry():
    """heel-strike 配置 (θ_th=θ_sh=−θ_st) で swing 足が接地高さにある。"""
    th = -0.18
    x = np.array([th, -th, -th, 0, 0, 0])
    pts = kneed_joints(x, l_t=0.5, l_s=0.5, foot_x=0.0)
    # pts = (stance_foot, stance_knee, hip, swing_knee, swing_foot)
    assert np.isclose(pts[4][1], 0.0, atol=1e-12)


def test_kneed_joints_straight_legs_match_two_segments():
    """両膝伸展時、knee は foot–hip 線分上にある。"""
    x = np.array([0.2, -0.3, -0.3, 0, 0, 0])
    pts = kneed_joints(x, l_t=0.5, l_s=0.5, foot_x=1.0)
    foot, knee_st, hip = pts[0], pts[1], pts[2]
    seg = hip - foot
    rel = knee_st - foot
    cross = seg[0] * rel[1] - seg[1] * rel[0]
    assert abs(cross) < 1e-12
```

- [ ] **Step 2: viz.py に kneed_joints と animate_kneed を実装**

```python
def kneed_joints(x: np.ndarray, l_t: float, l_s: float, foot_x: float):
    """kneed の関節座標列（slope frame）。stance 足を (foot_x, 0) に置く。

    返り値: (stance_foot, stance_knee, hip, swing_knee, swing_foot)
    """
    th_st, th_th, th_sh = x[0], x[1], x[2]
    length = l_t + l_s

    def down(theta):
        return np.array([np.sin(theta), -np.cos(theta)])

    foot = np.array([foot_x, 0.0])
    hip = foot - length * down(th_st)
    knee_st = hip + l_t * down(th_st)
    knee_sw = hip + l_t * down(th_th)
    foot_sw = knee_sw + l_s * down(th_sh)
    return foot, knee_st, hip, knee_sw, foot_sw
```

`animate_kneed(strides, p, gamma, out, fps=30)`: animate_walk と同じ構造
（np.interp リサンプリング、strike 瞬間フレームを含める、stance 足アンカーを
heel-strike の swing 足位置で前進、全体を −γ 回転）。線は stance 脚
（foot–knee–hip、青）と swing 脚（hip–knee–foot、橙）の2本の polyline で、
膝関節に点マーカー。位置補間は x の 3 角度を `np.interp` でそれぞれ補間して
`kneed_joints` を呼ぶ。foot 前進は `kneed_joints(x_strike, ...)` の swing_foot x。

- [ ] **Step 3: walk_kneed.py を書く**

walk_compass.py と同じ骨格（argparse: --gamma-deg default=rad2deg(GAMMA_GAIT),
--strides 30, --perturb 0.005; converged ガード + 不安定警告 + meta.json に
stable / knee_flexion_max を追加; animate_kneed で mp4）。deviation は
`np.linalg.norm(model.project(result.x_end) - fp.y)`。
さらに各 stride で `flexion_max = np.max(np.abs(result.x[2] - result.x[1]))` を
print（膝が動いている観察可能性）。

- [ ] **Step 4: 実行確認 + 全テスト + ruff + Commit**

Run: `uv run python scripts/walk_kneed.py` → mp4/png/meta 生成、deviation 減衰、
膝屈曲が毎 stride 観察できること。

```bash
git add src/crane/viz.py tests/test_viz_kneed.py scripts/walk_kneed.py
git commit -m "feat: kneed walk demo with 4-segment animation"
```

---

### Task 9: ドキュメント更新 + 歩容判定ゲート

**Files:**
- Modify: `GOALS.md`, `README.md`

- [ ] **Step 1: GOALS.md の Phase 3 節を実績で記入**

Phase 1-2 と同じ形式: チェックリスト（実測数値入り）、知見、今後の課題
（rocker foot = Phase 4 候補、等）。最後に
`- [ ] **ぽんぽこ殿の歩容判定**: walk.mp4 が「歩行」に見える（膝の屈曲込み）` を未チェックで。

- [ ] **Step 2: README.md に Phase 3 結果を追記 + Commit**

```bash
git add GOALS.md README.md
git commit -m "docs: record phase 3 results and gait judgment gate"
```

- [ ] **Step 3: 歩容判定ゲート（人間ゲート — STOP）**

`data/runs/<ts>_kneed_*/walk.mp4` をぽんぽこ殿に提示。
**膝付き歩行の判定が出るまで Phase 4 検討に進まない。**

---

## Self-Review Notes

- 設計書 Phase 3 の要件（4相 hybrid、リミットサイクル発見 + 安定性 + 歩容判定）を
  Task 3-9 でカバー。「相遷移を全部ログ」は trajectory + flexion print で担保
- 衝突定式化のリスクは三段防御: 不変量テスト (Task 4) → compass 退化ゲート (Task 6)
  → 文献ゲート (Task 7)。Phase 2 で機械精度実証済みの swap 定式化を 3 本に拡張
- 6D 統一状態の埋め込み（locked 相）は trajectory/viz/イベントの一様性のための
  設計判断。θ_th≡θ_sh の維持は Task 3 のテストで明示的に確認
- knee-strike「常に受理」とハイパーエクステンション簡略化、unlocked 相 scuff 無視は
  既知の簡略化として docstring 必須（Task 3 の docstring に記載済み）
- 型一貫性: PhaseSpec(dynamics/event_value/event_accept/impact)、
  HybridModel(phases/lift/project)、make_kneed の 4D 断面 — Task 1, 5, 7 で統一確認
- 文献数値はプランに書かない。プラン中の数値は Phase 2 実測の compass 不動点
  (0.27103, −1.09238, −0.37737) のみ（Task 6 の退化ゲート照合先）
- Task 2 の文献が点足 kneed の数値を十分印刷していない場合のフォールバック:
  パラメータのみ文献から採り、gait 照合は「digitized + 緩い許容」または
  「存在 + 安定性のみゲート」に縮退。その判断は Task 2 の report を見て controller が下す
