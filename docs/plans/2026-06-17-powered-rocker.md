# Phase 5b: 動力付き rocker_compass（push-off 増強）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** rocker_compass に pre-emptive 撃力 push-off を足し、固定勾配 γ=0.030 で push-off 増強サイクルを発見・検証する（5c の能動 vs 受動 basin 比較の土台）。

**Architecture:** rocker_compass.py の `_build` をミラーした `powered_rocker_compass.py`（rocker_compass.py 不変）。後脚軸方向の撃力モーメント `M_po` を角運動量保存の eq1 にのみ加える（撃力は hip を通るため eq2 は不変）。push_off=0 で rocker_compass に厳密一致。**push-off 公式は実装前 de-risk で検証済み**（push-off→0 で diff 0.0、エネルギー注入符号確定）。

**Tech Stack:** sympy, numpy, scipy, matplotlib（すべて既存）。新依存なし。

**設計ドキュメント:** `docs/2026-06-17-powered-rocker-design.md`

---

## 既存コードの要点

- `src/crane/models/rocker_compass.py`: `RockerCompassParams(m, m_h, c, rho, L, R, gamma, g)`、`_build()`（lru_cache、
  `f_qdd, f_energy, f_kinetic, f_impact` を返す）、`dynamics(t,x,p)`、`kinetic_energy(x,p)`、
  `heelstrike_map(x,p)`（衝突＋脚交換）、`make_rocker_compass(p)`。状態 `x=[θ_st,θ_sw,w_st,w_sw]`、
  断面 `y=(θ_st,w_st,w_sw)`、`lift=[y0,-y0,y1,y2]`、`project=[x0,x2,x3]`、event `x[0]+x[1]`、
  accept `x[0]<0 and x[2]+x[3]<0`。
- `from crane.derive.impact import angular_momentum, cross2d`、`from crane.derive.lagrange import derive_qdd`。
- **Phase 3.5 検証済み**: m=1, m_h=0, c=0.37, rho=0.32, L=1.0, R=0.3, γ=0.030 で
  y*=(0.30844, −1.26256, −0.87914), max|λ|=0.4316。
- `from crane.search import find_limit_cycle`; `from crane.stride import stride, StrideError`。

**de-risk で検証済みの push-off 導出**（そのまま使用）:
`Î_po = push_off·(P_st − hip)/|hip − P_st|`（後脚脚線方向、エネルギー注入符号）、
`M_po = cross2d(P_st − foot_sw_contact, Î_po)` を eq1 にのみ加える。push_off=0 で受動に厳密一致。

## File Structure
- Create: `src/crane/models/powered_rocker_compass.py` — rocker_compass の _build ミラー＋push-off。
- Create: `scripts/walk_powered_rocker.py` — γ=0.030 動力サイクル → 多歩 → 円弧足アニメ。
- Modify: `README.md`, `GOALS.md` — Phase 5b セクション。
- Tests: `test_powered_rocker_impact.py`, `test_powered_rocker_degenerate.py`, `test_powered_rocker_cycle.py`。

---

## Task 1: powered_rocker_compass.py（push-off ＋ factory）

**Files:**
- Create: `src/crane/models/powered_rocker_compass.py`
- Test: `tests/test_powered_rocker_impact.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_powered_rocker_impact.py
import numpy as np

from crane.models.rocker_compass import (
    RockerCompassParams,
    heelstrike_map as passive_heelstrike,
    kinetic_energy,
)
from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
    powered_heelstrike_map,
)

NOMINAL = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)


def test_pushoff_zero_equals_passive():
    pp = PoweredRockerCompassParams(**NOMINAL, push_off=0.0)
    pa = RockerCompassParams(**NOMINAL)
    for x in [
        np.array([0.30844, -0.30844, -1.26256, -0.87914]),
        np.array([0.25, -0.25, -1.1, -0.8]),
    ]:
        np.testing.assert_allclose(
            powered_heelstrike_map(x, pp), passive_heelstrike(x, pa), atol=1e-12
        )


def test_pushoff_injects_energy():
    # push_off>0 で post-collision KE が増加（エネルギー注入、de-risk で確認済み）
    pa = RockerCompassParams(**NOMINAL)
    x = np.array([0.30844, -0.30844, -1.26256, -0.87914])
    pp0 = PoweredRockerCompassParams(**NOMINAL, push_off=0.0)
    ppP = PoweredRockerCompassParams(**NOMINAL, push_off=0.2)
    ke0 = kinetic_energy(powered_heelstrike_map(x, pp0), pa)
    keP = kinetic_energy(powered_heelstrike_map(x, ppP), pa)
    assert keP > ke0


def test_factory_lift_project():
    pp = PoweredRockerCompassParams(**NOMINAL, push_off=0.1)
    model = make_powered_rocker_compass(pp)
    y = np.array([0.3, -1.2, -0.85])
    np.testing.assert_allclose(model.project(model.lift(y)), y, atol=1e-12)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_powered_rocker_impact.py -v`
Expected: FAIL（`ModuleNotFoundError: crane.models.powered_rocker_compass`）

- [ ] **Step 3: Write the implementation**

```python
# src/crane/models/powered_rocker_compass.py
"""動力付き rocker_compass。rocker_compass に pre-emptive 撃力 push-off を足す。

後脚（trailing = old stance）軸方向の撃力 Î_po = push_off·(P_st-hip)/|hip-P_st| を
heelstrike 衝突の直前に注入。撃力は接触点 P_st と hip を結ぶ脚線上にあり hip を通るため、
hip まわりのモーメントはゼロ。よって角運動量保存の eq1（系の新接触点まわり）にのみ
push-off モーメント M_po = cross2d(P_st - foot_sw_contact, Î_po) を加え、eq2（swing 脚の
hip まわり）は不変。push_off=0 で rocker_compass の衝突写像に厳密一致。

push-off 公式は実装前 de-risk で検証済み（push_off→0 で diff 0.0、push_off>0 でエネルギー注入）。
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum, cross2d
from crane.derive.lagrange import derive_qdd
from crane.model import HybridModel, PhaseSpec


@dataclass(frozen=True)
class PoweredRockerCompassParams:
    m: float
    m_h: float
    c: float
    rho: float
    L: float
    R: float
    gamma: float
    push_off: float  # pre-emptive 撃力 push-off の大きさ（0 で受動）
    g: float = 9.81


@lru_cache(maxsize=1)
def _build():
    th_st, th_sw, w_st, w_sw = sp.symbols("th_st th_sw w_st w_sw")
    m, m_h, c, rho, L, R, gamma, g, po = sp.symbols(
        "m m_h c rho L R gamma g po", positive=True
    )
    q, qd = [th_st, th_sw], [w_st, w_sw]

    def down(theta):
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    C_st = sp.Matrix([-R * th_st, R])
    hip = C_st - (L - R) * down(th_st)
    p_st1 = hip + (c - rho) * down(th_st)
    p_st2 = hip + (c + rho) * down(th_st)
    p_sw1 = hip + (c - rho) * down(th_sw)
    p_sw2 = hip + (c + rho) * down(th_sw)
    C_sw = hip + (L - R) * down(th_sw)
    foot_sw_contact = sp.Matrix([C_sw[0], 0])
    P_st = sp.Matrix([-R * th_st, 0])

    bodies = [
        (m_h, hip),
        (m / 2, p_st1),
        (m / 2, p_st2),
        (m / 2, p_sw1),
        (m / 2, p_sw2),
    ]

    def vel(r):
        return r.jacobian(q) * sp.Matrix(qd)

    T = sum(mass * vel(r).dot(vel(r)) for mass, r in bodies) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])
    V = -g * sum(mass * r.dot(g_dir) for mass, r in bodies)

    qdd = derive_qdd(q, qd, T, V, simplify=False)
    params = (m, m_h, c, rho, L, R, gamma, g, po)
    args = (th_st, th_sw, w_st, w_sw, *params)
    f_qdd = sp.lambdify(args, [qdd[0], qdd[1]], "numpy")
    f_kinetic = sp.lambdify(args, T, "numpy")

    st_masses = [(m / 2, p_st1), (m / 2, p_st2)]
    sw_masses = [(m / 2, p_sw1), (m / 2, p_sw2)]
    L_sys_stance = angular_momentum(bodies, q, qd, P_st)
    L_sys_swfoot = angular_momentum(bodies, q, qd, foot_sw_contact)
    L_st_hip = angular_momentum(st_masses, q, qd, hip)
    L_sw_hip = angular_momentum(sw_masses, q, qd, hip)

    # pre-emptive push-off: 後脚脚線方向の撃力（hip を通るため eq2 は不変、eq1 のみ）
    leg = hip - P_st
    legn = sp.sqrt(leg.dot(leg))
    impulse = po * (P_st - hip) / legn  # エネルギー注入符号（de-risk 検証済み）
    M_po = cross2d(P_st - foot_sw_contact, impulse)

    wp_st, wp_sw = sp.symbols("wp_st wp_sw")
    swap = {th_st: th_sw, th_sw: th_st, w_st: wp_st, w_sw: wp_sw}
    eqs = [
        L_sys_stance.subs(swap, simultaneous=True) - L_sys_swfoot - M_po,
        L_sw_hip.subs(swap, simultaneous=True) - L_st_hip,
    ]
    A, rhs = sp.linear_eq_to_matrix(eqs, [wp_st, wp_sw])
    qd_post = A.LUsolve(rhs)
    f_impact = sp.lambdify(args, [qd_post[0], qd_post[1]], "numpy")

    return f_qdd, f_kinetic, f_impact


def _args(x, p: PoweredRockerCompassParams):
    return (*x, p.m, p.m_h, p.c, p.rho, p.L, p.R, p.gamma, p.g, p.push_off)


def dynamics(t: float, x, p: PoweredRockerCompassParams):
    f_qdd = _build()[0]
    qdd = f_qdd(*_args(x, p))
    return [x[2], x[3], qdd[0], qdd[1]]


def powered_heelstrike_map(x: np.ndarray, p: PoweredRockerCompassParams) -> np.ndarray:
    """衝突写像（push-off モーメント付き）＋脚交換。push_off=0 で受動に一致。"""
    wp = _build()[2](*_args(x, p))
    return np.array([x[1], x[0], float(wp[0]), float(wp[1])])


def make_powered_rocker_compass(p: PoweredRockerCompassParams) -> HybridModel:
    return HybridModel(
        phases=(
            PhaseSpec(
                dynamics=lambda t, x: dynamics(t, x, p),
                event_value=lambda x: x[0] + x[1],
                event_accept=lambda x: x[0] < 0.0 and (x[2] + x[3]) < 0.0,
                impact=lambda x: powered_heelstrike_map(x, p),
            ),
        ),
        lift=lambda y: np.array([y[0], -y[0], y[1], y[2]]),
        project=lambda x: np.array([x[0], x[2], x[3]]),
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_powered_rocker_impact.py -v`
Expected: PASS（3 tests）。sympy ビルドで初回数秒。

- [ ] **Step 5: ruff then commit**

```bash
uv run ruff format src/crane/models/powered_rocker_compass.py tests/test_powered_rocker_impact.py
uv run ruff check src/crane/models/powered_rocker_compass.py tests/test_powered_rocker_impact.py
git add src/crane/models/powered_rocker_compass.py tests/test_powered_rocker_impact.py
git commit -m "$(cat <<'EOF'
feat: powered rocker_compass push-off impact (push_off=0 reduces to passive)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: push-off→0 退化ゲート

**Files:**
- Test: `tests/test_powered_rocker_degenerate.py`

- [ ] **Step 1: Write the test**

```python
# tests/test_powered_rocker_degenerate.py
import numpy as np

from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.search import find_limit_cycle
from crane.stride import stride

NOMINAL = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)
GUESS = np.array([0.30, -1.0, -0.40])


def test_single_stride_matches_passive_at_pushoff_zero():
    passive = make_rocker_compass(RockerCompassParams(**NOMINAL))
    powered = make_powered_rocker_compass(PoweredRockerCompassParams(**NOMINAL, push_off=0.0))
    x0 = passive.lift(np.array([0.30844, -1.26256, -0.87914]))
    np.testing.assert_allclose(
        stride(powered, x0).x_end, stride(passive, x0).x_end, atol=1e-10
    )


def test_limit_cycle_matches_phase35_at_pushoff_zero():
    powered = make_powered_rocker_compass(PoweredRockerCompassParams(**NOMINAL, push_off=0.0))
    fp = find_limit_cycle(powered, GUESS)
    assert fp.converged
    np.testing.assert_allclose(fp.y, [0.30844, -1.26256, -0.87914], atol=1e-4)
    assert np.max(np.abs(fp.eigenvalues)) < 1.0  # max|λ|≈0.4316
```

- [ ] **Step 2: Run + verify (no weakening)**

Run: `uv run pytest tests/test_powered_rocker_degenerate.py -v`
Expected: PASS（2 tests）。FAIL 時は push-off 導出のバグ。**緩和せず**原因調査（Phase 3.5 値
(0.30844, −1.26256, −0.87914) は検証済み ground truth）。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format tests/test_powered_rocker_degenerate.py
uv run ruff check tests/test_powered_rocker_degenerate.py
git add tests/test_powered_rocker_degenerate.py
git commit -m "$(cat <<'EOF'
test: powered rocker push-off->0 degenerate gate (matches Phase 3.5 passive)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: 動力サイクル ＋ 安定性（γ=0.030, push-off 掃引）

**Files:**
- Test: `tests/test_powered_rocker_cycle.py`

- [ ] **Step 1: Write the test**

```python
# tests/test_powered_rocker_cycle.py
import numpy as np

from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.search import find_limit_cycle

NOMINAL = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)
GUESS = np.array([0.30844, -1.26256, -0.87914])


def test_powered_cycles_exist_and_stable_across_pushoff():
    """γ=0.030 で push_off∈{0,0.04,0.08} の安定サイクルが受動 seed から直接収束。"""
    guess = GUESS
    results = {}
    for po in [0.0, 0.04, 0.08]:
        model = make_powered_rocker_compass(PoweredRockerCompassParams(**NOMINAL, push_off=po))
        fp = find_limit_cycle(model, guess)
        assert fp.converged, f"push_off={po} not converged"
        assert np.max(np.abs(fp.eigenvalues)) < 1.0, f"push_off={po} unstable"
        results[po] = fp.y
        guess = fp.y
    # push-off で不動点が動く（受動と異なる動力サイクル）
    assert not np.allclose(results[0.0], results[0.08], atol=1e-3)
```

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/test_powered_rocker_cycle.py -v`
Expected: PASS。push_off=0,0.04,0.08 すべて安定サイクル収束（de-risk で確認済み: max|λ| 0.43→0.46 程度）。
**安定性 assert は緩めない**。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format tests/test_powered_rocker_cycle.py
uv run ruff check tests/test_powered_rocker_cycle.py
git add tests/test_powered_rocker_cycle.py
git commit -m "$(cat <<'EOF'
test: powered rocker limit cycles stable across push-off at fixed slope

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: walk_powered_rocker.py ＋ 記録 ＋ 歩容判定

**Files:**
- Create: `scripts/walk_powered_rocker.py`
- Modify: `README.md`, `GOALS.md`

- [ ] **Step 1: walk_powered_rocker.py を書く**

`scripts/walk_rocker.py`（Phase 3.5）を参考に、円弧足アニメの呼び出しを合わせる。`walk_rocker.py` が
使う animation 関数（`animate_rocker` 等）を確認し、同じ関数・引数で呼ぶこと。骨子:

```python
# scripts/walk_powered_rocker.py
"""Phase 5b デモ: 動力付き rocker_compass（γ=0.030 + push-off 増強）の動力サイクル → 多歩 → アニメ。

Usage: uv run python scripts/walk_powered_rocker.py [--push-off 0.08] [--strides 30] [--perturb 0.005]
"""

import argparse
import json
import sys

import numpy as np

from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_rocker, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--push-off", type=float, default=0.08)
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.005)
    args = parser.parse_args()

    p = PoweredRockerCompassParams(
        m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030,
        push_off=args.push_off, g=9.81,
    )
    model = make_powered_rocker_compass(p)
    fp = find_limit_cycle(model, np.array([0.30844, -1.26256, -0.87914]))
    if not fp.converged:
        sys.exit("ERROR: 動力サイクルが見つからない（--push-off を調整）")
    lam = np.abs(fp.eigenvalues)
    print(f"converged (gamma=0.030, push_off={args.push_off}): y*={fp.y} |lambda|={lam}")

    run_dir = new_run_dir(f"powered_rocker_P{args.push_off:g}")
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
    animate_rocker(strides, p.L, p.R, p.gamma, run_dir / "walk.mp4")

    meta = {
        "gamma": 0.030, "push_off": args.push_off, "R": p.R,
        "fixed_point": fp.y.tolist(), "eigenvalues_abs": lam.tolist(),
        "stable": bool(lam.max() < 1.0), "n_strides_completed": len(strides),
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

注: `animate_rocker(strides, L, R, gamma, out)` は walk_rocker.py（Phase 3.5）と同一シグネチャ。

- [ ] **Step 2: スモーク実行**

Run: `uv run python scripts/walk_powered_rocker.py --push-off 0.08 --strides 20`
Expected: `converged ...` が出て、複数歩 deviation 減少、`walk.mp4`/`phase_portrait.png`/`meta.json` 生成。
実測 y*・max|λ| を報告。

- [ ] **Step 3: ruff + commit（スクリプト）**

```bash
uv run ruff format scripts/walk_powered_rocker.py && uv run ruff check scripts/walk_powered_rocker.py
git add scripts/walk_powered_rocker.py
git commit -m "$(cat <<'EOF'
feat: walk_powered_rocker.py demo (gamma=0.030 push-off augmentation)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 4: README + GOALS に Phase 5b セクション追加（実測値で埋める）**

README に「## Phase 5b の結果（動力付き rocker_compass, push-off 増強）」、GOALS に
「## Phase 5b: 動力付き rocker_compass（push-off 増強）— 完了 (2026-06-17)」を追加。ゲート:
```markdown
- [x] powered_rocker_compass.py（push-off モーメントを eq1 に追加、push_off=0 で受動退化）
- [x] **push-off→0 退化ゲート**: push_off=0 で Phase 3.5 検証済み (0.30844, −1.26256, −0.87914) に全 stride 一致
- [x] 動力サイクル γ=0.030・push_off∈[0,0.08] で安定（実測値を記録）
- [x] push-off がエネルギー注入（post-collision KE 増加、walkable slope 拡張）
- [ ] **ぽんぽこ殿の目視判定**: walk.mp4 が「円弧足で転がりつつ push-off で歩く」に見える
```
知見: rocker は push-off でも完全な平地 γ=0 には届かず fold で消失（R=0.3）。4a.1 の車輪極限 fold と
響き合う。続き Phase 5c（R 掃引で能動 vs 受動 basin 比較＝本命仮説）を今後の課題に。目視判定は UNCHECKED。

- [ ] **Step 5: 全テスト + ruff + commit（docs）**

```bash
uv run pytest -q
uv run ruff format --check . && uv run ruff check .
git add README.md GOALS.md
git commit -m "$(cat <<'EOF'
docs: record phase 5b powered rocker results and judgment gate

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6: 目視判定を依頼**

walk.mp4 をコントローラ経由でぽんぽこ殿に提示し判定を仰ぐ。合格後に GOALS のチェックを埋め PR へ。

---

## Self-Review チェック結果
- **Spec coverage:** 設計 A〜F を網羅。B push-off→Task1（de-risk 検証済み公式）、C 固定 γ サイクル→Task3、
  D ゲート（退化→Task2、サイクル+安定→Task3、エネルギー注入→Task1、目視→Task4）、E モジュール→Task1/4、F 命名→Task4。
- **Placeholder scan:** Task4 の円弧足アニメ呼び出しのみ「walk_rocker.py に合わせる」とした（既存 viz 再利用で確実、
  関数名は実装者が walk_rocker.py で確認）。他になし。
- **Type consistency:** `PoweredRockerCompassParams(m,m_h,c,rho,L,R,gamma,push_off,g=9.81)`、
  `powered_heelstrike_map(x,p)`、`make_powered_rocker_compass(p)`、_build が `(f_qdd,f_kinetic,f_impact)` を返す、
  `_args` が push_off を末尾に含む、が全タスクで一貫。rocker_compass.py は不変（mirror）。
- **de-risk 反映:** push-off 公式（eq1 のみ M_po、エネルギー注入符号）は scratch で検証済み（push_off→0 diff 0.0）。
