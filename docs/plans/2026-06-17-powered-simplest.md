# Phase 5a: 動力付き simplest walker（Kuo 2002）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** simplest walker に pre-emptive 撃力 push-off を足し、平地 γ=0 を push-off 駆動で歩くリミットサイクルを発見・検証する。

**Architecture:** simplest.py をミラーした新 `powered_simplest.py`（simplest.py は不変）。連続相は simplest の dynamics を γ 経由で再利用、heel-strike だけを合成写像（push-off→衝突→脚交換）に差し替え。push-off 撃力 P が制御パラメータで、P=0 で受動 simplest に厳密退化。既存の find_limit_cycle / stride / animate_walk を再利用。

**Tech Stack:** numpy, scipy, matplotlib（すべて既存）。新依存なし。

**設計ドキュメント:** `docs/2026-06-17-powered-simplest-design.md`

---

## 既存コードの要点

- `src/crane/models/simplest.py`: `SimplestParams(gamma)`、`dynamics(t, x, p)`（連続相、`theta_dd=sin(theta-gamma)`）、
  `heelstrike_map(x)`（受動衝突: `c=cos(2θ)` として `[-θ, -2θ, c·θ̇, c(1-c)·θ̇]`）、`lift(theta, theta_dot)`、
  `make_simplest(p)`。状態 `x=[θ, φ, θ̇, φ̇]`、断面 `y=(θ, θ̇)`、project=`[x0, x2]`。
- `crane.model.HybridModel`, `PhaseSpec`。
- `crane.search.find_limit_cycle(model, y_guess) -> FixedPoint(.y, .eigenvalues, .converged)`。
- `crane.stride.stride(model, x) -> StrideResult(.x_end, .x_strike, .t, .x, ...)`、`StrideError`。
- `crane.viz.animate_walk`, `plot_phase_portrait`（model 非依存、strides を取る）。`crane.runs.new_run_dir`。
- `crane import references as ref`: `ref.GAMMA_REF=0.009`, `ref.LONG_PERIOD_THETA=0.199529`, `ref.LONG_PERIOD_THETA_DOT=-0.198983`。

**push-off 写像（設計で導出済み、reduced 座標 M=l=g=1）**: `c=cos(2θ), s=sin(2θ)`、
`θ̇⁺ = c·θ̇ + s·P`、返り値 `[-θ, -2θ, θ̇⁺, (1-c)·θ̇⁺]`。P=0 で受動に厳密一致。
push-off 仕事 = P²/2（脚軸方向の初速 0 に撃力 P を与えるため、KE 増分 = P²/2）。

## File Structure
- Create: `src/crane/references_kuo.py` — Kuo 2002 provenance ＋公表数値。
- Create: `src/crane/models/powered_simplest.py` — push-off 合成衝突写像＋factory（simplest をミラー）。
- Create: `scripts/walk_powered.py` — 平地動力サイクル発見→多歩→アニメ（既存 viz 再利用）。
- Modify: `README.md`, `GOALS.md` — Phase 5a セクション。
- Tests: `test_powered_simplest_impact.py`, `test_powered_simplest_degenerate.py`,
  `test_powered_simplest_cycle.py`, `test_references_kuo.py`。

---

## Task 1: powered_simplest.py（push-off 合成衝突写像 ＋ factory）

**Files:**
- Create: `src/crane/models/powered_simplest.py`
- Test: `tests/test_powered_simplest_impact.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_powered_simplest_impact.py
import numpy as np

from crane.models.simplest import heelstrike_map as passive_heelstrike
from crane.models.powered_simplest import (
    PoweredSimplestParams,
    make_powered_simplest,
    powered_heelstrike_map,
)


def test_pushoff_zero_equals_passive_map():
    # P=0 で受動 heelstrike_map に厳密一致（退化）
    for x in [
        np.array([-0.2, -0.4, -0.2, 0.1]),
        np.array([-0.15, -0.3, -0.25, 0.05]),
    ]:
        np.testing.assert_allclose(
            powered_heelstrike_map(x, 0.0), passive_heelstrike(x), atol=1e-14
        )


def test_pushoff_adds_sin2theta_P_to_stance_velocity():
    # θ̇⁺(P) − θ̇⁺(0) = sin(2θ)·P（push-off 項）
    x = np.array([-0.2, -0.4, -0.2, 0.1])
    P = 0.1
    out0 = powered_heelstrike_map(x, 0.0)
    outP = powered_heelstrike_map(x, P)
    s = np.sin(2.0 * x[0])
    assert np.isclose(outP[2] - out0[2], s * P, atol=1e-14)
    # 4 番目（swing 速度）は (1-c)·θ̇⁺ の関係を保つ
    c = np.cos(2.0 * x[0])
    assert np.isclose(outP[3], (1.0 - c) * outP[2], atol=1e-14)


def test_factory_builds_model_with_passive_dynamics():
    p = PoweredSimplestParams(gamma=0.0, push_off=0.1)
    model = make_powered_simplest(p)
    assert len(model.phases) == 1
    # lift/project が受動と同形
    y = np.array([0.2, -0.2])
    x = model.lift(y)
    np.testing.assert_allclose(model.project(x), y, atol=1e-14)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_powered_simplest_impact.py -v`
Expected: FAIL（`ModuleNotFoundError: crane.models.powered_simplest`）

- [ ] **Step 3: Write the implementation**

```python
# src/crane/models/powered_simplest.py
"""Kuo 2002 の動力付き simplest walker。

simplest walker（Garcia 1998）に pre-emptive 撃力 push-off を足し、平地 γ=0 を歩く。
連続相は simplest.dynamics を γ 経由で再利用し、heel-strike だけを
「後脚軸方向の撃力 P → 角運動量保存衝突 → 脚交換」の合成写像に差し替える。

push-off 写像（reduced 座標, M=l=g=1）: c=cos(2θ), s=sin(2θ),
    θ̇⁺ = c·θ̇ + s·P
P=0 で受動 heelstrike_map に厳密退化。push-off 仕事 = P²/2。
"""

from dataclasses import dataclass

import numpy as np

from crane.model import HybridModel, PhaseSpec
from crane.models.simplest import SimplestParams, dynamics, lift


@dataclass(frozen=True)
class PoweredSimplestParams:
    gamma: float  # slope angle [rad]（平地は 0）
    push_off: float  # pre-emptive 撃力 push-off の大きさ P（P=0 で受動）


def powered_heelstrike_map(x: np.ndarray, push_off: float) -> np.ndarray:
    """合成衝突写像: 後脚軸撃力 push_off → 角運動量保存衝突 → 脚交換。

    P=0 で simplest.heelstrike_map に厳密一致。
    """
    theta, _phi, theta_dot, _phi_dot = x
    c = np.cos(2.0 * theta)
    s = np.sin(2.0 * theta)
    td_plus = c * theta_dot + s * push_off
    return np.array([-theta, -2.0 * theta, td_plus, (1.0 - c) * td_plus])


def make_powered_simplest(p: PoweredSimplestParams) -> HybridModel:
    """パラメータ束縛済みの HybridModel。連続相は simplest と同一（γ 経由）。"""
    sp = SimplestParams(gamma=p.gamma)
    return HybridModel(
        phases=(
            PhaseSpec(
                dynamics=lambda t, x: dynamics(t, x, sp),
                event_value=lambda x: x[1] - 2.0 * x[0],
                event_accept=lambda x: x[0] < 0.0 and (x[3] - 2.0 * x[2]) > 0.0,
                impact=lambda x: powered_heelstrike_map(x, p.push_off),
            ),
        ),
        lift=lambda y: lift(y[0], y[1]),
        project=lambda x: np.array([x[0], x[2]]),
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_powered_simplest_impact.py -v`
Expected: PASS（3 tests）

- [ ] **Step 5: ruff then commit**

```bash
uv run ruff format src/crane/models/powered_simplest.py tests/test_powered_simplest_impact.py
uv run ruff check src/crane/models/powered_simplest.py tests/test_powered_simplest_impact.py
git add src/crane/models/powered_simplest.py tests/test_powered_simplest_impact.py
git commit -m "$(cat <<'EOF'
feat: powered simplest walker push-off impact map (P=0 reduces to passive)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: push-off→0 退化ゲート（最強検証）

**Files:**
- Test: `tests/test_powered_simplest_degenerate.py`

powered model が P=0 で Phase 1 検証済み受動 simplest walker に全 stride・不動点・固有値で一致することを担保。

- [ ] **Step 1: Write the test**

```python
# tests/test_powered_simplest_degenerate.py
import numpy as np

from crane import references as ref
from crane.models.powered_simplest import PoweredSimplestParams, make_powered_simplest
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane.stride import stride


def test_single_stride_matches_passive_at_pushoff_zero():
    gamma = ref.GAMMA_REF
    passive = make_simplest(SimplestParams(gamma=gamma))
    powered = make_powered_simplest(PoweredSimplestParams(gamma=gamma, push_off=0.0))
    x0 = passive.lift(np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    r_passive = stride(passive, x0)
    r_powered = stride(powered, x0)
    np.testing.assert_allclose(r_powered.x_end, r_passive.x_end, atol=1e-12)


def test_limit_cycle_matches_phase1_at_pushoff_zero():
    gamma = ref.GAMMA_REF
    powered = make_powered_simplest(PoweredSimplestParams(gamma=gamma, push_off=0.0))
    fp = find_limit_cycle(
        powered, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    )
    assert fp.converged
    # Phase 1 検証済み long-period 不動点・固有値に一致
    np.testing.assert_allclose(fp.y, [0.2003109, -0.1998325], atol=1e-5)
    assert np.max(np.abs(fp.eigenvalues)) < 1.0  # max|λ|≈0.589 安定
```

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/test_powered_simplest_degenerate.py -v`
Expected: PASS（2 tests）。

FAIL する場合は push-off 写像か factory にバグ。**緩和せず**原因調査（assert を弱めない）。
Phase 1 の long-period 不動点は README 記載 (0.2003109, −0.1998325)。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format tests/test_powered_simplest_degenerate.py
uv run ruff check tests/test_powered_simplest_degenerate.py
git add tests/test_powered_simplest_degenerate.py
git commit -m "$(cat <<'EOF'
test: powered simplest push-off->0 degenerate gate (matches Phase 1 passive)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Kuo 2002 provenance（references_kuo.py）

**Files:**
- Create: `src/crane/references_kuo.py`
- Test: `tests/test_references_kuo.py`

- [ ] **Step 1: 一次資料を取り寄せて記録**

WebSearch/WebFetch で **Kuo, A.D. (2002) "Energetics of actively powered locomotion using the simplest
walking model", J. Biomech. Eng. 124(1):113–120** を当たり、push-off 駆動 simplest walker について
公表されている関係（push-off 仕事 ∝ push-off²、力学コスト ∝ 速度²、最適 push-off の効率など）と、
照合可能な具体値があれば取得する。`src/crane/references_kuo.py` に provenance 付きで記録:

```python
# src/crane/references_kuo.py
"""Kuo (2002) "Energetics of actively powered locomotion using the simplest
walking model", J. Biomech. Eng. 124(1):113-120 の公表値。Phase 5a 検証ゲート。

URL: <実際に参照した URL>  取得日: 2026-06-17
"""

# push-off 駆動 simplest walker の公表知見（原典から取得できたもののみ。
# 取得できなければ None とし、定性関係＋エネルギー収支ゲートに委譲。捏造禁止）。
PROVENANCE: str = "<参照 URL と確認した記述の要約、または取得失敗の明記>"

# push-off 仕事 ∝ push-off²（本実装は厳密に P²/2）。Kuo の力学コスト∝速度²関係:
COT_SCALES_WITH_SPEED_SQUARED: bool = True  # 定性関係（原典の記述を確認して設定）

# 照合可能な数値（あれば。窓/条件付き。なければ None）:
PUBLISHED_VALUE: float | None = None
PUBLISHED_VALUE_DESC: str = "<何の値か。なければ空>"
```

**記憶からの数値はゲートに使わない（CLAUDE.md）。** Kuo 本文が paywall 等で取得不能なら、その旨を
PROVENANCE に明記し、数値は None のままとする（ゲートは退化＋エネルギー収支に委譲）。

- [ ] **Step 2: Write the provenance test**

```python
# tests/test_references_kuo.py
from crane import references_kuo as ref


def test_kuo_provenance_present():
    assert isinstance(ref.PROVENANCE, str) and len(ref.PROVENANCE) > 0
    # 数値を主張するなら provenance に URL か取得失敗の明記があること
    assert ("http" in ref.PROVENANCE) or ("未取得" in ref.PROVENANCE) or ("取得失敗" in ref.PROVENANCE)
```

- [ ] **Step 3: Run + ruff + commit**

```bash
uv run pytest tests/test_references_kuo.py -v
uv run ruff format src/crane/references_kuo.py tests/test_references_kuo.py
uv run ruff check src/crane/references_kuo.py tests/test_references_kuo.py
git add src/crane/references_kuo.py tests/test_references_kuo.py
git commit -m "$(cat <<'EOF'
docs: record Kuo 2002 provenance for powered simplest walker

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: 平地リミットサイクル ＋ 安定性 ＋ エネルギー収支ゲート

**Files:**
- Test: `tests/test_powered_simplest_cycle.py`

- [ ] **Step 1: Write the test**

```python
# tests/test_powered_simplest_cycle.py
import numpy as np

from crane.models.powered_simplest import PoweredSimplestParams, make_powered_simplest
from crane.search import find_limit_cycle
from crane.stride import stride


def _ke(theta_dot):
    # simplest walker の KE = (1/2)θ̇²（hip 点質量 M=1, 脚 massless, l=1）
    return 0.5 * theta_dot**2


def _find_level_cycle():
    """平地 γ=0 の動力サイクルを探す。受動不動点近傍を seed に P を選ぶ。"""
    # 受動 long-period 不動点 (0.2003, -0.1998) 近傍から探索。
    # P を小さくスキャンして最初に収束する平地サイクルを取る。
    for push_off in [0.06, 0.07, 0.08, 0.09, 0.10, 0.12]:
        model = make_powered_simplest(PoweredSimplestParams(gamma=0.0, push_off=push_off))
        fp = find_limit_cycle(model, np.array([0.18, -0.18]))
        if fp.converged:
            return push_off, model, fp
    raise AssertionError("平地動力サイクルが見つからない")


def test_level_ground_limit_cycle_exists_and_stable():
    push_off, model, fp = _find_level_cycle()
    assert fp.converged
    assert np.max(np.abs(fp.eigenvalues)) < 1.0  # 安定
    # 平地（γ=0）でも前進する: θ̇* < 0（前進方向）
    assert fp.y[1] < 0.0


def test_energy_balance_pushoff_work_equals_collision_loss():
    """平地サイクル上で push-off 仕事 P²/2 = 一歩衝突損失。"""
    push_off, model, fp = _find_level_cycle()
    result = stride(model, model.lift(fp.y))
    theta_strike = result.x_strike[0]
    td_strike = result.x_strike[2]
    td_end = result.x_end[2]
    # 衝突損失 = (push-off 後 KE) − (衝突後 KE)
    ke_after_pushoff = 0.5 * (td_strike**2 + push_off**2)
    collision_loss = ke_after_pushoff - _ke(td_end)
    pushoff_work = 0.5 * push_off**2
    assert np.isclose(pushoff_work, collision_loss, atol=1e-3)
```

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/test_powered_simplest_cycle.py -v`
Expected: PASS（2 tests）。実 solve_ivp 積分のため数秒。

もし `_find_level_cycle` が全 P で収束しない場合: seed を受動不動点 `[0.2003, -0.1998]` に変える、
または P の範囲を広げる（`[0.04, ..., 0.16]`）。それでも収束しなければ continuation
（γ=0.009 の受動サイクルから γ を 0 へ小刻みに下げ、各 γ で前解を seed に push_off を増やす）を
`_find_level_cycle` 内に実装する。**エネルギー収支 assert は緩めない**（バグの検出器）。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format tests/test_powered_simplest_cycle.py
uv run ruff check tests/test_powered_simplest_cycle.py
git add tests/test_powered_simplest_cycle.py
git commit -m "$(cat <<'EOF'
test: powered simplest level-ground limit cycle + energy balance gate

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: walk_powered.py ＋ 記録 ＋ 歩容判定

**Files:**
- Create: `scripts/walk_powered.py`
- Modify: `README.md`, `GOALS.md`

- [ ] **Step 1: walk_powered.py を書く**

```python
# scripts/walk_powered.py
"""Phase 5a デモ: 動力付き simplest walker の平地リミットサイクル → 多歩シミュ → アニメ。

Usage: uv run python scripts/walk_powered.py [--push-off 0.08] [--strides 30] [--perturb 0.01]
"""

import argparse
import json
import sys

import numpy as np

from crane.models.powered_simplest import PoweredSimplestParams, make_powered_simplest
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_walk, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--push-off", type=float, default=0.08)
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.01)
    args = parser.parse_args()

    p = PoweredSimplestParams(gamma=0.0, push_off=args.push_off)
    model = make_powered_simplest(p)
    run_dir = new_run_dir(f"powered_P{args.push_off:g}")

    fp = find_limit_cycle(model, np.array([0.18, -0.18]))
    if not fp.converged:
        for i, (y, r) in enumerate(fp.history):
            print(f"  newton[{i}] y={y} residual={r:.3e}")
        sys.exit("ERROR: 平地動力サイクルが見つからない（--push-off を調整）")
    lam = np.abs(fp.eigenvalues)
    print(f"converged: y*={fp.y}  |lambda|={lam}  (gamma=0 平地)")
    if lam.max() > 1.0:
        print(f"WARNING: UNSTABLE (max|lambda|={lam.max():.3f})")

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
    animate_walk(strides, gamma=0.0, out=run_dir / "walk.mp4")

    meta = {
        "gamma": 0.0,
        "push_off": args.push_off,
        "pushoff_work_per_step": 0.5 * args.push_off**2,
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": lam.tolist(),
        "stable": bool(lam.max() < 1.0),
        "n_strides_completed": len(strides),
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

> 注: `animate_walk` の引数は既存の simplest デモ `scripts/walk_simplest.py` の呼び出しに合わせること。
> もし `animate_walk(strides, gamma, out)` の正確なシグネチャが異なる場合は walk_simplest.py を確認して合わせる。

- [ ] **Step 2: 平地歩行をスモーク実行**

Run: `uv run python scripts/walk_powered.py --push-off 0.08 --strides 20`
Expected: `converged: y*=...` が出て、平地（gamma=0）で複数歩 deviation が減少（収束）し、
`walk.mp4` / `phase_portrait.png` / `meta.json` が生成。FELL が続く場合は --push-off を 0.06〜0.12 で調整。
実測の y*・push-off 値を報告。

- [ ] **Step 3: README に Phase 5a セクションを追加**

`README.md` の Phase 4a.1 の後に「Phase 5a の結果（動力付き simplest walker, Kuo 2002）」を追加。
**Step 2 の実測値で埋める**: 平地 γ=0 で歩く push-off 値・不動点・max|λ|、push-off→0 退化ゲート
（Phase 1 一致）、エネルギー収支（push-off 仕事 P²/2 = 衝突損失）、Kuo 2002 照合の可否。
使い方セクションに walk_powered.py の行を追加。data/runs/ は gitignore のローカル成果物と明記。

- [ ] **Step 4: GOALS に Phase 5a セクションを追加**

`GOALS.md` の Phase 4a.1 の後に「Phase 5a: 動力付き simplest walker（Kuo 2002）」を追加。ゲート:
```markdown
- [x] powered_simplest.py（push-off 合成衝突写像、P=0 で受動退化）
- [x] **push-off→0 退化ゲート**: P=0 で Phase 1 検証済み受動サイクル (0.2003109, −0.1998325) に全 stride 一致
- [x] 平地 γ=0 リミットサイクル発見 + 安定性 max|λ|<1（実測値を記録）
- [x] エネルギー収支: push-off 仕事 P²/2 = 一歩衝突損失
- [x] Kuo 2002 provenance 記録（取得可否を honest に）
- [ ] **ぽんぽこ殿の目視判定**: 平地 walk.mp4 が「push-off で前進する歩行」に見える
```
能動歩行アークの第一段（平地歩行＝受動には不可能だった能力）達成、として記録。
最後の目視判定は UNCHECKED で残す（コントローラが依頼）。続き Phase 5b（能動 rocker_compass 仮説検証）を今後の課題に。

- [ ] **Step 5: 全テスト + ruff + commit**

```bash
uv run pytest -q
uv run ruff format --check . && uv run ruff check .
git add README.md GOALS.md
git commit -m "$(cat <<'EOF'
docs: record phase 5a powered simplest walker results and judgment gate

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6: 目視判定を依頼**

`walk.mp4` をコントローラ経由でぽんぽこ殿に提示し「平地を push-off で前進する歩行に見えるか」の判定を仰ぐ。
合格後に GOALS のチェックを埋め、PR へ。

---

## Self-Review チェック結果

- **Spec coverage:** 設計 A〜F を網羅。B push-off 物理→Task1、C 平地サイクル→Task4、
  D ゲート（退化→Task2、Kuo→Task3、エネルギー収支→Task4、安定性→Task4、目視→Task5）、
  E モジュール→Task1/3/5、F 命名→Task5。
- **Placeholder scan:** Kuo 2002 の数値のみ「実装時取得」（CLAUDE.md 準拠、取得不能時の honest フォールバック明記、捏造禁止）。他になし。
- **Type consistency:** `PoweredSimplestParams(gamma, push_off)`、`powered_heelstrike_map(x, push_off)`、
  `make_powered_simplest(p)`、push-off 写像 `θ̇⁺=c·θ̇+s·P` / 返り値 `[-θ,-2θ,θ̇⁺,(1-c)θ̇⁺]`、
  KE=(1/2)θ̇²、push-off 仕事=P²/2 が全タスクで一貫。simplest.py は不変（import で再利用のみ）。
- **既知の留意点:** `animate_walk` のシグネチャは walk_simplest.py に合わせる旨を Task5 に明記。
  平地サイクルの seed/push-off は収束しなければ範囲拡大→continuation のフォールバックを Task4 に明記。
