# Crane Phase 0–1: Simplest Walker Limit Cycle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Garcia 1998 "The Simplest Walking Model" のリミットサイクルを Poincaré shooting で発見し、文献値と照合し、歩容アニメーションを生成する。

**Architecture:** Hybrid dynamics（連続相 = `solve_ivp` + event 検出、離散相 = 衝突写像）で stride 写像 S を純関数として構築。Poincaré 断面（heel-strike 直後、reduced 2D 座標）上で S の不動点をニュートン法で解き、数値 Jacobian の固有値で安定性を判定する。

**Tech Stack:** Python 3.12 / uv / numpy / scipy / matplotlib / pytest / ruff

**Spec:** `docs/2026-06-12-crane-design.md`

---

## 物理モデル早見表（Garcia et al. 1998, m/M → 0 極限）

実装者向けの要約。**Task 3 で原論文と照合し、相違があればこの表と実装を修正すること。**

- 状態ベクトル: `x = [θ, φ, θ̇, φ̇]`
  - θ: stance 脚の角度（slope 法線から、無次元時間スケール済み）
  - φ: stance 脚から測った swing 脚の相対角度
- 運動方程式（脚長 l=1, 重力 g=1 にスケール済み）:
  - `θ̈ = sin(θ − γ)`
  - `φ̈ = sin(θ − γ) + θ̇² sin φ − cos(θ − γ) sin φ`
- Heel-strike 条件: `φ − 2θ = 0`（かつ θ < 0、θ > 0 側の交差は foot scuffing として無視）
- Heel-strike 写像（脚交換込み）:
  - `θ⁺ = −θ⁻`
  - `φ⁺ = −2θ⁻`
  - `θ̇⁺ = cos(2θ⁻) θ̇⁻`
  - `φ̇⁺ = cos(2θ⁻)(1 − cos(2θ⁻)) θ̇⁻`
- Poincaré 断面 = heel-strike 直後。このとき `φ = 2θ` かつ `φ̇ = (1 − cos 2θ) θ̇` が
  自動的に成立するため、断面上の自由座標は **y = (θ, θ̇) の 2 次元**。
- 既知の性質: γ = 0.009 に long-period / short-period の 2 解。long-period は小傾斜で
  安定（θ* ≈ 0.2 オーダー）、θ* ∝ γ^(1/3) スケーリング。**正確な数値は Task 3 で取得。**

## File Structure

| パス | 責務 |
|---|---|
| `pyproject.toml` | uv プロジェクト定義、ruff / pytest 設定 |
| `CLAUDE.md` | Claude Code への指示書（Heron 流儀の引き継ぎ） |
| `src/crane/models/simplest.py` | Garcia model: params / dynamics / 衝突写像 / lift |
| `src/crane/references.py` | 文献の公表数値（provenance 付き、Task 3 で記入） |
| `src/crane/stride.py` | stride 写像（積分 + イベント + 衝突の合成） |
| `src/crane/search.py` | Poincaré map / Newton shooting / 固有値解析 |
| `src/crane/viz.py` | phase portrait / stick-figure アニメーション |
| `src/crane/runs.py` | `data/runs/<timestamp>_<tag>/` 出力ディレクトリ管理 |
| `scripts/walk_simplest.py` | Phase 1 デモ: 発見→多歩シミュ→アニメ mp4 |
| `tests/test_*.py` | 各モジュールのテスト（文献値ゲート含む） |

---

### Task 1: プロジェクト雛形（uv + ruff + pytest + src layout）

**Files:**
- Create: `pyproject.toml`, `.gitignore`, `CLAUDE.md`, `README.md`, `src/crane/__init__.py`

- [ ] **Step 1: pyproject.toml を書く**

```toml
[project]
name = "crane"
version = "0.1.0"
description = "Passive dynamic walking limit cycle discovery (math-first revival of Heron)"
requires-python = ">=3.12"
dependencies = [
    "numpy>=2.0",
    "scipy>=1.14",
    "matplotlib>=3.9",
]

[dependency-groups]
dev = [
    "pytest>=8.0",
    "ruff>=0.8",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/crane"]

[tool.ruff]
line-length = 100

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 2: .gitignore を書く**

```gitignore
.venv/
__pycache__/
*.egg-info/
.pytest_cache/
.ruff_cache/
data/
```

- [ ] **Step 3: CLAUDE.md を書く**

```markdown
# Claude Codeへの指示書 (Crane)

受動歩行機械のリミットサイクルを数学的に発見する研究プロジェクト。
Heron (~/git/heron) の後継。設計は docs/2026-06-12-crane-design.md を参照。

## 言語と呼びかけ
- 回答は日本語で。ユーザーは「ぽんぽこ殿」と呼ぶ
- 技術的な深さは妥協しない

## 開発スタイル
- 1タスク1イテレーション1コミット。コミットメッセージは英語 Conventional Commits
- パッケージ管理は uv、フォーマッタ/リンタは ruff
- 新しい依存を入れる時は事前に許可を取る
- dataclass + 型ヒント積極活用、ログ・収束履歴を惜しまない

## 検証の流儀
- テストカバレッジではなく「文献の公表数値との一致」が各 Phase のゲート
- 文献値は必ず provenance (URL, 取得日) 付きで src/crane/references.py に記録
- Claude の記憶からの数値をゲートに使わない（Heron 教訓）
- 歩容の最終判定はぽんぽこ殿の動画観察（観察フェーズ）

## やってほしくないこと
- 大きめフレームワークの導入、過剰な抽象化、ファイルの増やしすぎ
- 受動歩行の分野知識を勝手に類推して仮定しない（ぽんぽこ殿が研究者）
```

- [ ] **Step 4: README.md スタブと src/crane/__init__.py を書く**

`README.md`:

```markdown
# Crane

受動歩行機械のリミットサイクルを hybrid dynamics + Poincaré shooting で
数学的に発見する研究プロジェクト。Heron の後継。

- 設計: `docs/2026-06-12-crane-design.md`
- セットアップ: `uv sync`
- テスト: `uv run pytest`

詳細は Phase 1 完了時に拡充予定。
```

`src/crane/__init__.py`: 空ファイル。

- [ ] **Step 5: uv sync で環境構築し、ruff を通す**

Run: `cd /Users/yast/git/crane && uv sync && uv run ruff check . && uv run pytest`
Expected: `uv sync` 成功、ruff エラー 0 件、pytest は exit code 5（"no tests ran" — この段階では正常）

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml .gitignore CLAUDE.md README.md src/ uv.lock
git commit -m "chore: scaffold crane project (uv + ruff + pytest, src layout)"
```

---

### Task 2: Phase 0 ゲート — solve_ivp イベント検出スモークと runs ヘルパー

**Files:**
- Create: `src/crane/runs.py`
- Test: `tests/test_phase0.py`

- [ ] **Step 1: イベント検出の失敗しないことを確認するテストを書く**

調和振動子 `ẍ = −x`, `x(0)=1, ẋ(0)=0` の最初のゼロ交差は解析的に `t = π/2`。
event 検出がこの精度で動くことが Phase 0 のゲート。

```python
# tests/test_phase0.py
import numpy as np
from scipy.integrate import solve_ivp

from crane.runs import new_run_dir


def test_event_detection_quarter_period():
    """solve_ivp の event 検出が解析解 t=π/2 と一致する。"""

    def f(t, x):
        return [x[1], -x[0]]

    def zero_cross(t, x):
        return x[0]

    zero_cross.terminal = True
    zero_cross.direction = -1

    sol = solve_ivp(f, (0.0, 10.0), [1.0, 0.0], events=zero_cross, rtol=1e-10, atol=1e-12)
    t_event = sol.t_events[0][0]
    assert abs(t_event - np.pi / 2) < 1e-8


def test_new_run_dir(tmp_path):
    """タイムスタンプ付き run ディレクトリが作られ、重複しない。"""
    d1 = new_run_dir("smoke", base=tmp_path)
    assert d1.is_dir()
    assert "smoke" in d1.name
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_phase0.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.runs'`

- [ ] **Step 3: runs.py を実装**

```python
# src/crane/runs.py
"""data/runs/<timestamp>_<tag>/ 出力ディレクトリ管理。"""

from datetime import datetime
from pathlib import Path


def new_run_dir(tag: str, base: Path | str = "data/runs") -> Path:
    """タイムスタンプ付き run ディレクトリを作って返す。上書きしない。"""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = Path(base) / f"{stamp}_{tag}"
    path.mkdir(parents=True, exist_ok=False)
    return path
```

- [ ] **Step 4: テストを実行して全部通ることを確認**

Run: `uv run pytest tests/test_phase0.py -v`
Expected: 2 PASSED

- [ ] **Step 5: Commit**

```bash
git add src/crane/runs.py tests/test_phase0.py
git commit -m "feat: phase 0 gate - solve_ivp event detection smoke + run dir helper"
```

---

### Task 3: 文献値の取得と記録（references.py）

**Files:**
- Create: `src/crane/references.py`
- Test: `tests/test_references.py`

**この Task は Web 調査を含む。** WebSearch / WebFetch で Garcia, Chatterjee, Ruina,
Coleman (1998) "The Simplest Walking Model: Stability, Complexity, and Scaling"
(ASME Journal of Biomechanical Engineering 120(2)) の本文または信頼できる転載
（著者サイト ruina.tam.cornell.edu 等）を取得し、以下を行う:

1. **運動方程式と heel-strike 写像を本プランの「物理モデル早見表」と照合**。
   相違があれば早見表と Task 4 のコードを修正し、相違内容をコミットメッセージに記録
2. γ = 0.009 における long-period / short-period 両解の不動点数値を取得
3. 安定性に関する公表値（固有値、または安定な γ 範囲、period-doubling 開始点）を取得
4. θ* ∝ γ^(1/3) スケーリング則の記述を確認

- [ ] **Step 1: references.py の枠組みとゲートテストを書く**

```python
# tests/test_references.py
from crane import references as ref


def test_references_are_filled():
    """文献値が provenance 付きで記入済みであること（Phase 1 ゲートの前提）。

    isinstance 検査なのは、未記入の Ellipsis (...) が `is not None` を
    すり抜けるため。
    """
    assert ref.PROVENANCE.startswith("Garcia")
    assert "http" in ref.PROVENANCE
    assert ref.GAMMA_REF == 0.009
    assert isinstance(ref.LONG_PERIOD_THETA, float)
    assert isinstance(ref.LONG_PERIOD_THETA_DOT, float)
    assert 0.0 < ref.LONG_PERIOD_THETA < 0.5
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_references.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.references'`

- [ ] **Step 3: Web 調査を実施し、references.py を実数で記入**

ファイル形式（数値は必ず論文から転記する。下記 `...` は調査で得た実数に置き換え、
記載桁数も論文どおりにする）:

```python
# src/crane/references.py
"""Garcia et al. 1998 の公表数値。Phase 1 検証ゲートの照合対象。

取得した数値の桁数は論文の記載どおり。論文に無い項目は None とし、
その旨をコメントに残す。
"""

PROVENANCE = (
    "Garcia, Chatterjee, Ruina, Coleman (1998), "
    "'The Simplest Walking Model: Stability, Complexity, and Scaling', "
    "ASME J. Biomech. Eng. 120(2). Retrieved 2026-06-XX from <URL>"
)

GAMMA_REF = 0.009

# γ=0.009 の long-period gait 不動点 (heel-strike 直後の断面座標)
LONG_PERIOD_THETA: float = ...      # stance angle θ* [rad]
LONG_PERIOD_THETA_DOT: float = ...  # θ̇*

# γ=0.009 の short-period gait 不動点
SHORT_PERIOD_THETA: float | None = ...
SHORT_PERIOD_THETA_DOT: float | None = ...

# 安定性: 論文が報告する形式に合わせて記録（固有値 or 安定 γ 範囲）
LONG_PERIOD_EIGENVALUE_ABS: tuple[float, ...] | None = ...
STABLE_GAMMA_MAX: float | None = ...  # long-period gait が安定な γ の上限

# スケーリング則 θ* ∝ γ^(1/3)
SCALING_EXPONENT = 1.0 / 3.0
```

- [ ] **Step 4: テストを実行して通ることを確認**

Run: `uv run pytest tests/test_references.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/crane/references.py tests/test_references.py
git commit -m "feat: record Garcia 1998 published values with provenance"
```

（早見表との相違が見つかった場合は body に `EOM cross-check: <内容>` を追記）

---

### Task 4: models/simplest.py — 運動方程式と衝突写像（TDD）

**Files:**
- Create: `src/crane/models/__init__.py`, `src/crane/models/simplest.py`
- Test: `tests/test_simplest.py`

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_simplest.py
import numpy as np

from crane.models.simplest import SimplestParams, dynamics, heelstrike_map, lift


P = SimplestParams(gamma=0.009)


def test_equilibrium_state_has_zero_acceleration():
    """θ=γ, φ=0, 静止状態では加速度ゼロ（直立平衡）。"""
    x = [P.gamma, 0.0, 0.0, 0.0]
    xdot = dynamics(0.0, x, P)
    assert np.allclose(xdot, [0.0, 0.0, 0.0, 0.0], atol=1e-15)


def test_heelstrike_geometry_is_exact():
    """θ⁺=−θ⁻, φ⁺=−2θ⁻ は厳密（脚交換の鏡映）。"""
    x_minus = np.array([-0.2, -0.4, -0.25, -0.01])
    x_plus = heelstrike_map(x_minus)
    assert x_plus[0] == 0.2
    assert x_plus[1] == 0.4


def test_heelstrike_velocity_map():
    """θ̇⁺ = cos(2θ⁻)θ̇⁻, φ̇⁺ = cos(2θ⁻)(1−cos(2θ⁻))θ̇⁻。"""
    theta, theta_dot = -0.2, -0.25
    x_plus = heelstrike_map(np.array([theta, 2 * theta, theta_dot, -0.01]))
    c = np.cos(2 * theta)
    assert np.isclose(x_plus[2], c * theta_dot)
    assert np.isclose(x_plus[3], c * (1 - c) * theta_dot)


def test_post_impact_state_is_on_section():
    """衝突直後は φ=2θ かつ φ̇=(1−cos2θ)θ̇ が成立 = lift と整合。"""
    x_plus = heelstrike_map(np.array([-0.2, -0.4, -0.25, -0.01]))
    expected = lift(x_plus[0], x_plus[2])
    assert np.allclose(x_plus, expected)
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_simplest.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.models'`

- [ ] **Step 3: simplest.py を実装**

`src/crane/models/__init__.py` は空ファイル。

```python
# src/crane/models/simplest.py
"""Garcia et al. 1998 'The Simplest Walking Model' (m/M → 0 極限、l=g=1 スケール)。

状態 x = [θ, φ, θ̇, φ̇]。θ: stance 脚角度（slope 法線基準）、
φ: stance 脚から測った swing 脚相対角度。
heel-strike 条件は φ − 2θ = 0 (θ < 0 側)。
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SimplestParams:
    gamma: float  # slope angle [rad]


def dynamics(t: float, x, p: SimplestParams):
    """連続相（swing phase）の運動方程式。"""
    theta, phi, theta_dot, _phi_dot = x
    theta_dd = np.sin(theta - p.gamma)
    phi_dd = theta_dd + theta_dot**2 * np.sin(phi) - np.cos(theta - p.gamma) * np.sin(phi)
    return [theta_dot, _phi_dot, theta_dd, phi_dd]


def heelstrike_map(x: np.ndarray) -> np.ndarray:
    """heel-strike の衝突写像（角運動量保存）+ 脚交換。"""
    theta, _phi, theta_dot, _phi_dot = x
    c = np.cos(2.0 * theta)
    return np.array([-theta, -2.0 * theta, c * theta_dot, c * (1.0 - c) * theta_dot])


def lift(theta: float, theta_dot: float) -> np.ndarray:
    """Poincaré 断面座標 y=(θ, θ̇) を全状態へ持ち上げる（衝突直後の拘束）。"""
    return np.array(
        [theta, 2.0 * theta, theta_dot, (1.0 - np.cos(2.0 * theta)) * theta_dot]
    )
```

- [ ] **Step 4: テストを実行して通ることを確認**

Run: `uv run pytest tests/test_simplest.py -v`
Expected: 4 PASSED

- [ ] **Step 5: Commit**

```bash
git add src/crane/models/ tests/test_simplest.py
git commit -m "feat: simplest walker EOM, heelstrike map, section lift (Garcia 1998)"
```

---

### Task 5: stride.py — 1歩写像（積分 + scuffing ガード + 衝突）

**Files:**
- Create: `src/crane/stride.py`
- Test: `tests/test_stride.py`

**設計メモ:** 衝突直後の状態は event 関数 `g = φ − 2θ` の零点上にあり、また mid-swing
に θ > 0 での偽交差（foot scuffing）がある。対策:
1. 各積分区間の冒頭 `t_burn = 1e-3`（無次元時間。1歩 ≈ 3.8 に対し十分小）は
   event なしで積分してから event 付き積分に入る
2. event 発火時に `θ ≥ 0` なら scuffing とみなし、その地点から再積分（同じ burn-in）

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_stride.py
import numpy as np
import pytest

from crane import references as ref
from crane.models.simplest import SimplestParams, lift
from crane.stride import StrideError, stride


P = SimplestParams(gamma=ref.GAMMA_REF)


def test_stride_from_reference_guess_returns_to_section():
    """文献不動点近傍から1歩進むと、終状態が断面 (φ=2θ) 上にある。"""
    x0 = lift(ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT)
    result = stride(P, x0)
    assert result.t_step > 1.0
    assert np.isclose(result.x_end[1], 2.0 * result.x_end[0], atol=1e-8)


def test_stride_records_trajectory():
    """軌跡が時系列で記録される（観察可能性）。"""
    x0 = lift(ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT)
    result = stride(P, x0)
    assert result.t.shape[0] == result.x.shape[1]
    assert result.t[0] == 0.0
    assert np.all(np.diff(result.t) >= 0)


def test_stride_raises_when_no_heelstrike():
    """heel-strike に至らない初期条件では StrideError。"""
    x0 = lift(0.001, 0.0)  # ほぼ直立静止 → 歩かない
    with pytest.raises(StrideError):
        stride(P, x0)
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_stride.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.stride'`

- [ ] **Step 3: stride.py を実装**

```python
# src/crane/stride.py
"""stride 写像: 衝突直後状態から積分 → heel-strike → 衝突写像 → 次の衝突直後状態。"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from crane.models.simplest import SimplestParams, dynamics, heelstrike_map

T_BURN = 1e-3  # event 無効の冒頭区間（衝突直後 g=0 と scuffing 直後の再発火対策）


class StrideError(RuntimeError):
    """heel-strike に到達しなかった（転倒・停止など）。"""


@dataclass
class StrideResult:
    x_end: np.ndarray  # 衝突写像適用後（次の断面上の状態）
    x_strike: np.ndarray  # 衝突写像適用前（heel-strike 瞬間の状態）
    t_step: float
    t: np.ndarray  # 軌跡時刻列
    x: np.ndarray  # 軌跡状態列 shape (4, N)


def _heelstrike_event(t, x, p):
    return x[1] - 2.0 * x[0]


_heelstrike_event.terminal = True
_heelstrike_event.direction = 0


def stride(
    p: SimplestParams,
    x0: np.ndarray,
    *,
    t_max: float = 10.0,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> StrideResult:
    """1歩分の写像。副作用なし。"""
    t0 = 0.0
    x = np.asarray(x0, dtype=float)
    ts: list[np.ndarray] = []
    xs: list[np.ndarray] = []

    while t0 < t_max:
        # burn-in: event なしで微小区間を進める
        burn = solve_ivp(dynamics, (t0, t0 + T_BURN), x, args=(p,), rtol=rtol, atol=atol)
        ts.append(burn.t)
        xs.append(burn.y)
        t0, x = burn.t[-1], burn.y[:, -1]

        sol = solve_ivp(
            dynamics,
            (t0, t_max),
            x,
            args=(p,),
            events=_heelstrike_event,
            rtol=rtol,
            atol=atol,
        )
        ts.append(sol.t)
        xs.append(sol.y)

        if sol.t_events[0].size == 0:
            break  # t_max まで heel-strike なし

        t_e = sol.t_events[0][0]
        x_e = sol.y_events[0][0]
        if x_e[0] < 0.0:  # 本物の heel-strike（stance 脚が鉛直を越えている）
            return StrideResult(
                x_end=heelstrike_map(x_e),
                x_strike=x_e,
                t_step=t_e,
                t=np.concatenate(ts),
                x=np.concatenate(xs, axis=1),
            )
        t0, x = t_e, x_e  # scuffing: 無視して続行

    raise StrideError(f"no heelstrike before t_max={t_max}")
```

- [ ] **Step 4: テストを実行して通ることを確認**

Run: `uv run pytest tests/test_stride.py -v`
Expected: 3 PASSED

（注: `test_stride_from_reference_guess_returns_to_section` が FAIL する場合、
文献不動点からの1歩が成立していない = EOM か衝突写像の転記ミスの可能性が高い。
Task 3 の照合結果に立ち返って早見表との相違を探すこと。）

- [ ] **Step 5: Commit**

```bash
git add src/crane/stride.py tests/test_stride.py
git commit -m "feat: stride map with heelstrike event and scuffing guard"
```

---

### Task 6: search.py — Poincaré shooting と固有値解析

**Files:**
- Create: `src/crane/search.py`
- Test: `tests/test_search.py`

- [ ] **Step 1: 失敗するテストを書く（文献値ゲート本体）**

```python
# tests/test_search.py
import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams
from crane.search import find_limit_cycle


P = SimplestParams(gamma=ref.GAMMA_REF)


def test_long_period_fixed_point_matches_garcia1998():
    """Phase 1 ゲート: 不動点が文献値と一致（相対誤差 1e-3 以内）。"""
    guess = np.array([ref.LONG_PERIOD_THETA * 1.05, ref.LONG_PERIOD_THETA_DOT * 0.95])
    fp = find_limit_cycle(P, guess)
    assert fp.converged
    assert np.isclose(fp.y[0], ref.LONG_PERIOD_THETA, rtol=1e-3)
    assert np.isclose(fp.y[1], ref.LONG_PERIOD_THETA_DOT, rtol=1e-3)


def test_long_period_gait_is_stable():
    """γ=0.009 の long-period gait は漸近安定 (max|λ| < 1)。"""
    fp = find_limit_cycle(P, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    assert fp.converged
    assert np.max(np.abs(fp.eigenvalues)) < 1.0
    if ref.LONG_PERIOD_EIGENVALUE_ABS is not None:
        ours = np.sort(np.abs(fp.eigenvalues))
        theirs = np.sort(np.array(ref.LONG_PERIOD_EIGENVALUE_ABS))
        assert np.allclose(ours, theirs, rtol=0.05)


def test_convergence_history_is_logged():
    """収束履歴が残る（観察可能性）。残差は単調減少に近い挙動を示す。"""
    guess = np.array([ref.LONG_PERIOD_THETA * 1.05, ref.LONG_PERIOD_THETA_DOT * 0.95])
    fp = find_limit_cycle(P, guess)
    residuals = [r for _, r in fp.history]
    assert len(residuals) >= 2
    assert residuals[-1] < residuals[0]
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_search.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.search'`

- [ ] **Step 3: search.py を実装**

```python
# src/crane/search.py
"""Poincaré 写像の不動点探索（Newton shooting）と安定性解析。"""

from dataclasses import dataclass, field

import numpy as np

from crane.models.simplest import SimplestParams, lift
from crane.stride import stride


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


def find_limit_cycle(
    p: SimplestParams,
    y_guess: np.ndarray,
    *,
    tol: float = 1e-12,
    max_iter: int = 30,
) -> FixedPoint:
    """Newton 法で S(y) = y を解く。収束履歴を必ず残す。"""
    y = np.asarray(y_guess, dtype=float).copy()
    history: list[tuple[np.ndarray, float]] = []

    for _ in range(max_iter):
        residual = poincare_map(p, y) - y
        norm = float(np.linalg.norm(residual))
        history.append((y.copy(), norm))
        if norm < tol:
            J = _jacobian(p, y)
            return FixedPoint(y=y, eigenvalues=np.linalg.eigvals(J), converged=True, history=history)
        J = _jacobian(p, y)
        y = y - np.linalg.solve(J - np.eye(y.size), residual)

    return FixedPoint(y=y, eigenvalues=None, converged=False, history=history)
```

- [ ] **Step 4: テストを実行して通ることを確認**

Run: `uv run pytest tests/test_search.py -v`
Expected: 3 PASSED

（FAIL する場合の切り分け: (a) 収束しない → 初期推定と Newton step を history で確認、
(b) 収束するが文献値とズレる → EOM / 衝突写像 / 断面の定義の相違。Task 3 へ立ち返る）

- [ ] **Step 5: Commit**

```bash
git add src/crane/search.py tests/test_search.py
git commit -m "feat: Poincare shooting with eigenvalue stability analysis

Gate: fixed point matches Garcia 1998 published values"
```

---

### Task 7: γ^(1/3) スケーリング則の検証（continuation の芽）

**Files:**
- Test: `tests/test_scaling.py`

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_scaling.py
import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams
from crane.search import find_limit_cycle


def test_stance_angle_scales_as_gamma_one_third():
    """θ* ∝ γ^(1/3) (Garcia 1998)。continuation で γ を変えながら追跡。"""
    gammas = [0.004, 0.006, 0.009, 0.012]
    thetas = []
    y = np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    # γ=0.009 から出発し、近い γ へ順に continuation
    for gamma in sorted(gammas, key=lambda g: abs(g - ref.GAMMA_REF)):
        fp = find_limit_cycle(SimplestParams(gamma=gamma), y)
        assert fp.converged, f"no convergence at gamma={gamma}"
        y = fp.y  # 次の γ の初期推定に使う
        thetas.append((gamma, fp.y[0]))

    thetas.sort()
    log_g = np.log([g for g, _ in thetas])
    log_t = np.log([t for _, t in thetas])
    slope = np.polyfit(log_g, log_t, 1)[0]
    assert abs(slope - ref.SCALING_EXPONENT) < 0.05
```

- [ ] **Step 2: テストを実行する**

Run: `uv run pytest tests/test_scaling.py -v`
Expected: PASS（実装は Task 6 までで揃っている。FAIL なら物理が壊れている兆候なので
原因を調査し、必要なら Task 3-6 を見直す）

- [ ] **Step 3: Commit**

```bash
git add tests/test_scaling.py
git commit -m "test: verify theta* ~ gamma^(1/3) scaling law via continuation"
```

---

### Task 8: viz.py — phase portrait と stick-figure アニメーション

**Files:**
- Create: `src/crane/viz.py`
- Test: `tests/test_viz.py`

- [ ] **Step 1: 失敗するテストを書く**

アニメは目視判定が本番なので、テストはファイル生成と座標変換の整合のみ。

```python
# tests/test_viz.py
import numpy as np

from crane.viz import link_points


def test_link_points_heelstrike_symmetry():
    """heel-strike 瞬間 (φ=2θ) には swing 足が接地高さ (y=0) にある。"""
    theta = -0.2
    hip, swing = link_points(theta, 2 * theta, foot_x=0.0)
    assert np.isclose(hip[1], np.cos(theta))
    assert np.isclose(swing[1], 0.0, atol=1e-12)


def test_link_points_upright():
    """θ=0, φ=0 で両脚とも鉛直、swing 足は stance 足と重なる。"""
    hip, swing = link_points(0.0, 0.0, foot_x=1.5)
    assert np.allclose(hip, [1.5, 1.0])
    assert np.allclose(swing, [1.5, 0.0])
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_viz.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crane.viz'`

- [ ] **Step 3: viz.py を実装**

```python
# src/crane/viz.py
"""phase portrait と stick-figure アニメーション（slope frame 描画）。"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from crane.stride import StrideResult


def link_points(theta: float, phi: float, foot_x: float) -> tuple[np.ndarray, np.ndarray]:
    """slope frame での hip / swing 足の座標 (脚長 l=1)。

    stance 足を (foot_x, 0) に置く。脚角 ψ の足先は hip + (sin ψ, −cos ψ)。
    stance: ψ=θ、swing: ψ=θ−φ。
    """
    hip = np.array([foot_x - np.sin(theta), np.cos(theta)])
    psi = theta - phi
    swing = hip + np.array([np.sin(psi), -np.cos(psi)])
    return hip, swing


def plot_phase_portrait(strides: list[StrideResult], out: Path) -> None:
    """θ–θ̇ 平面の軌道。リミットサイクルなら閉軌道が重なる。"""
    fig, ax = plt.subplots(figsize=(6, 5))
    for s in strides:
        ax.plot(s.x[0], s.x[2], lw=0.8)
    ax.set_xlabel("theta [rad]")
    ax.set_ylabel("theta_dot [rad/s]")
    ax.set_title(f"Phase portrait ({len(strides)} strides)")
    fig.savefig(out, dpi=150)
    plt.close(fig)


def animate_walk(strides: list[StrideResult], gamma: float, out: Path, fps: int = 30) -> None:
    """stick-figure アニメ。slope frame で描き、全体を −γ 回転して坂を見せる。"""
    rot = np.array([[np.cos(-gamma), -np.sin(-gamma)], [np.sin(-gamma), np.cos(-gamma)]])

    # 各 stride の stance 足アンカーを heel-strike 位置で更新しながらフレーム列を作る
    frames: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    foot_x = 0.0
    dt_frame = 1.0 / fps
    for s in strides:
        t_resampled = np.arange(0.0, s.t[-1], dt_frame)
        for ti in t_resampled:
            idx = int(np.searchsorted(s.t, ti))
            idx = min(idx, s.t.size - 1)
            theta, phi = s.x[0, idx], s.x[1, idx]
            hip, swing = link_points(theta, phi, foot_x)
            foot = np.array([foot_x, 0.0])
            frames.append((rot @ foot, rot @ hip, rot @ swing))
        # heel-strike: swing 足の着地点が次の stance アンカー
        _, swing_end = link_points(s.x_strike[0], s.x_strike[1], foot_x)
        foot_x = float(swing_end[0])

    fig, ax = plt.subplots(figsize=(8, 4))
    (stance_line,) = ax.plot([], [], "o-", lw=2, color="tab:blue")
    (swing_line,) = ax.plot([], [], "o-", lw=2, color="tab:orange")
    span = max(foot_x + 2.0, 4.0)
    ground = rot @ np.array([[-1.0, span], [0.0, 0.0]])
    ax.plot(ground[0], ground[1], "k-", lw=1)
    ax.set_xlim(-1.0, span)
    ax.set_ylim(-span * np.sin(gamma) - 0.5, 1.5)
    ax.set_aspect("equal")

    def update(i: int):
        foot, hip, swing = frames[i]
        stance_line.set_data([foot[0], hip[0]], [foot[1], hip[1]])
        swing_line.set_data([hip[0], swing[0]], [hip[1], swing[1]])
        return stance_line, swing_line

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=True)
    try:
        anim.save(out, writer=animation.FFMpegWriter(fps=fps))
    except (FileNotFoundError, RuntimeError):
        # ffmpeg が無い環境では GIF にフォールバック
        out = out.with_suffix(".gif")
        anim.save(out, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
```

- [ ] **Step 4: テストを実行して通ることを確認**

Run: `uv run pytest tests/test_viz.py -v`
Expected: 2 PASSED

- [ ] **Step 5: Commit**

```bash
git add src/crane/viz.py tests/test_viz.py
git commit -m "feat: phase portrait and stick-figure walk animation"
```

---

### Task 9: scripts/walk_simplest.py — Phase 1 デモスクリプト

**Files:**
- Create: `scripts/walk_simplest.py`

- [ ] **Step 1: スクリプトを書く**

```python
# scripts/walk_simplest.py
"""Phase 1 デモ: リミットサイクル発見 → 摂動から多歩シミュ → アニメ生成。

Usage: uv run python scripts/walk_simplest.py [--gamma 0.009] [--strides 30] [--perturb 0.01]
"""

import argparse
import json

import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams, lift
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_walk, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma", type=float, default=ref.GAMMA_REF)
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.01, help="不動点への相対摂動")
    args = parser.parse_args()

    p = SimplestParams(gamma=args.gamma)
    run_dir = new_run_dir(f"simplest_g{args.gamma:g}")

    fp = find_limit_cycle(p, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    print(f"converged={fp.converged}  y*={fp.y}  |lambda|={np.abs(fp.eigenvalues)}")
    for i, (y, r) in enumerate(fp.history):
        print(f"  newton[{i}] y={y} residual={r:.3e}")

    # 不動点からわずかに摂動した初期条件で多歩シミュ（basin 内なら収束するはず）
    y0 = fp.y * (1.0 + args.perturb)
    x = lift(y0[0], y0[1])
    strides = []
    for i in range(args.strides):
        try:
            result = stride(p, x)
        except StrideError as e:
            print(f"stride {i}: FELL ({e})")
            break
        deviation = float(np.linalg.norm([result.x_end[0], result.x_end[2]] - fp.y))
        print(f"stride {i}: t={result.t_step:.4f} deviation={deviation:.3e}")
        strides.append(result)
        x = result.x_end

    plot_phase_portrait(strides, run_dir / "phase_portrait.png")
    animate_walk(strides, p.gamma, run_dir / "walk.mp4")

    meta = {
        "params": {"gamma": p.gamma},
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": np.abs(fp.eigenvalues).tolist(),
        "n_strides_completed": len(strides),
        "perturb": args.perturb,
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 実行して出力を確認**

Run: `uv run python scripts/walk_simplest.py`
Expected:
- `converged=True`、固有値の絶対値が全て < 1
- 30 strides 完走、deviation が単調減少傾向（リミットサイクルへの収束）
- `data/runs/<ts>_simplest_g0.009/` に `walk.mp4`, `phase_portrait.png`, `meta.json`

- [ ] **Step 3: ruff と全テストの最終確認**

Run: `uv run ruff format . && uv run ruff check . && uv run pytest -v`
Expected: format 差分なし or 自動修正、check エラー 0、全テスト PASS

- [ ] **Step 4: Commit**

```bash
git add scripts/walk_simplest.py
git commit -m "feat: phase 1 demo script - limit cycle discovery to walk animation"
```

---

### Task 10: ドキュメント更新と歩容判定ゲート

**Files:**
- Create: `GOALS.md`
- Modify: `README.md`

- [ ] **Step 1: GOALS.md を書く**

```markdown
# ゴール設定 (Crane)

設計は docs/2026-06-12-crane-design.md。各 Phase のゲートを満たすまで次へ進まない。

## Phase 0: 雛形と数値基盤
- [x] uv / ruff / pytest / src layout
- [x] solve_ivp イベント検出が解析解と 1e-8 で一致

## Phase 1: Simplest Walker (Garcia 1998)
- [x] EOM・衝突写像の実装と論文との照合
- [x] 文献値の provenance 付き記録 (references.py)
- [x] Poincaré shooting で不動点発見、文献値と一致 (rtol 1e-3)
- [x] 固有値で漸近安定性を確認 (max|λ| < 1)
- [x] θ* ∝ γ^(1/3) スケーリング則を再現
- [ ] **ぽんぽこ殿の歩容判定**: data/runs の walk.mp4 が「歩行」に見える

## Phase 2: sympy 導出レイヤー + Compass Gait (McGeer 1990)
Phase 1 ゲート通過後に計画策定。

## Phase 3: Kneed Walker (McGeer 1990)
Phase 2 完了後。
```

（チェック状態は実際のタスク完了状況に合わせて記入すること）

- [ ] **Step 2: README.md を拡充**

セットアップ、`scripts/walk_simplest.py` の使い方、Heron との関係、
docs/ への参照を含める。Phase 1 の結果数値（発見した不動点・固有値）を記載。

- [ ] **Step 3: Commit**

```bash
git add GOALS.md README.md
git commit -m "docs: add GOALS with phase gates, expand README with phase 1 results"
```

- [ ] **Step 4: 歩容判定ゲート（人間ゲート — STOP）**

`data/runs/<ts>_simplest_g0.009/walk.mp4` をぽんぽこ殿に提示し、
「これは歩行に見えるか」の判定を仰ぐ。**判定が出るまで Phase 2 計画に進まない。**

---

## Self-Review Notes

- スペックの Phase 0 / Phase 1 ゲートはそれぞれ Task 2 / Task 3-10 がカバー
- 文献値を Claude の記憶からゲートに使わない方針は Task 3 の Web 調査 +
  provenance 必須で担保（早見表はあくまで実装の出発点で、Task 3 で照合される）
- 型・命名の一貫性: `SimplestParams` / `stride()` / `StrideResult` / `lift()` /
  `find_limit_cycle()` / `FixedPoint` を全タスクで統一済み
- `references.py` の `...` は研究タスクの出力欄であり、Task 3 Step 3 で実数記入が
  必須（ゲートテストが None / Ellipsis を弾く）
