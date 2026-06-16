# Phase 4a.1: 受動歩行→車輪の連続軸（R-continuation）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 固定勾配 γ=0.030 で rocker_compass の円弧足半径 R を R∈(0,L=1) で掃引し、相対損失 δ(R)・安定性 max|λ|(R)・basin_fraction(R) を測って「コンパスは車輪への漸近」を検証する。

**Architecture:** 既存の `find_limit_cycle`（continuation seed）・`stride`（x_strike/x_end）・`rocker_compass.kinetic_energy`・`basin_slice`（Phase 4a）を再利用。新規は (a) model 非依存の効率指標モジュール `efficiency.py`、(b) 二段掃引スクリプト `wheel_limit.py`、(c) rimless-wheel の provenance。

**Tech Stack:** numpy, scipy, matplotlib, multiprocessing（すべて既存）。新依存なし。

**設計ドキュメント:** `docs/2026-06-16-wheel-limit-design.md`

---

## 既存コードの要点（実装者向け）

- `rocker_compass`: 状態 `x = [θ_st, θ_sw, w_st, w_sw]`（4D）。`make_rocker_compass(p)`、
  `kinetic_energy(x, p) -> float`、`energy(x, p) -> float`。
  project: `x -> [x0, x2, x3]`（断面 y=(θ_st, w_st, w_sw)）。
- `RockerCompassParams(m, m_h, c, rho, L, R, gamma, g=9.81)`（frozen dataclass）。
- references: `from crane import references_mcgeer as ref` →
  `ref.M_LEG=1.0, ref.M_HIP=0.0, ref.C_HIP_TO_COM=0.37, ref.RHO_GYR=0.32, ref.L_LEG=1.0,
  ref.R_FOOT=0.3, ref.GAMMA_GAIT=0.030, ref.G=9.81, ref.SECTION_GUESS=(0.30, -1.0, -0.40)`。
- `find_limit_cycle(model, y_guess) -> FixedPoint(.y, .eigenvalues, .converged)`（search.py）。
- `stride(model, x) -> StrideResult(.x_end, .x_strike, .t_step, .t, .x)`、転倒時 `StrideError`。
- `basin_slice(make_fn, params, fixed_point, *, axes, half_widths, resolution, model_name,
  max_strides=20, converge_tol=1e-3, n_workers=None) -> BasinResult(.grid, .basin_fraction, ...)`。
- **Phase 3.5 検証済みアンカー**: R=0.3, γ=0.030 で y\* = (0.30844, −1.26256, −0.87914)、max|λ|=0.4316。

## File Structure

- Create: `src/crane/efficiency.py` — model 非依存の効率指標（相対損失 δ、衝突損失、機械的 COT）。
- Create: `scripts/wheel_limit.py` — 二段 R 掃引（continuation + δ/COT + basin）、曲線・JSON・montage 出力。
- Modify: `src/crane/references_mcgeer.py` — rimless-wheel provenance（McGeer 1990a）。
- Create: `tests/test_efficiency.py` — 効率指標の単体。
- Create: `tests/test_wheel_limit_anchor.py` — R=0.3 アンカー + 連続性ゲート。

---

## Task 1: 効率指標モジュール `efficiency.py`

**Files:**
- Create: `src/crane/efficiency.py`
- Test: `tests/test_efficiency.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_efficiency.py
import numpy as np

from crane.efficiency import mechanical_cot, relative_loss, step_collision_loss
from crane.models.rocker_compass import RockerCompassParams, kinetic_energy, make_rocker_compass
from crane.search import find_limit_cycle
from crane.stride import stride
from crane import references_mcgeer as ref


def test_relative_loss_basic():
    # δ = (KE_pre - KE_post)/KE_pre
    assert relative_loss(1.0, 0.75) == 0.25
    assert relative_loss(2.0, 2.0) == 0.0


def test_mechanical_cot_basic():
    # COT = loss / (m g d)
    assert np.isclose(mechanical_cot(0.5, m=1.0, g=10.0, step_length=0.5), 0.1)


def test_step_collision_loss_on_rocker_limit_cycle():
    # 実モデルの 1 歩で 0 < δ < 1、KE_pre > KE_post（衝突は散逸）。
    p = RockerCompassParams(
        m=ref.M_LEG, m_h=ref.M_HIP, c=ref.C_HIP_TO_COM, rho=ref.RHO_GYR,
        L=ref.L_LEG, R=ref.R_FOOT, gamma=ref.GAMMA_GAIT, g=ref.G,
    )
    model = make_rocker_compass(p)
    fp = find_limit_cycle(model, np.array(ref.SECTION_GUESS))
    assert fp.converged
    result = stride(model, model.lift(fp.y))
    loss, ke_pre, ke_post = step_collision_loss(
        result.x_strike, result.x_end, lambda x: kinetic_energy(x, p)
    )
    assert ke_pre > ke_post > 0.0  # 衝突で KE 減少
    assert loss > 0.0
    delta = relative_loss(ke_pre, ke_post)
    assert 0.0 < delta < 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_efficiency.py -v`
Expected: FAIL (`ModuleNotFoundError: crane.efficiency`)

- [ ] **Step 3: Write minimal implementation**

```python
# src/crane/efficiency.py
"""受動歩行の効率指標（model 非依存）。

主指標は一歩あたり相対損失 δ = (KE_pre - KE_post)/KE_pre（車輪極限で →0）。
機械的 COT はリミットサイクル上で勾配に固定される（≈ sin γ）ので副次の内部チェック用。
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def relative_loss(ke_pre: float, ke_post: float) -> float:
    """一歩あたり相対損失 δ = (KE_pre - KE_post)/KE_pre。車輪極限で →0。"""
    return (ke_pre - ke_post) / ke_pre


def step_collision_loss(
    x_strike: np.ndarray,
    x_end: np.ndarray,
    kinetic_energy: Callable[[np.ndarray], float],
) -> tuple[float, float, float]:
    """heel-strike の衝突損失。(loss, ke_pre, ke_post) を返す。

    x_strike: 衝突直前の全状態、x_end: 衝突写像適用後（次断面）の全状態。
    KE は脚ラベル交換で不変なので ke_post は衝突直後 KE に等しい。
    """
    ke_pre = float(kinetic_energy(x_strike))
    ke_post = float(kinetic_energy(x_end))
    return ke_pre - ke_post, ke_pre, ke_post


def mechanical_cot(loss: float, *, m: float, g: float, step_length: float) -> float:
    """機械的 cost of transport = 衝突損失 / (m·g·一歩水平距離)。

    リミットサイクル上ではエネルギー収支より ≈ sin γ になるはず（内部チェック）。
    """
    return loss / (m * g * step_length)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_efficiency.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: ruff then commit**

```bash
uv run ruff format src/crane/efficiency.py tests/test_efficiency.py && uv run ruff check src/crane/efficiency.py tests/test_efficiency.py
git add src/crane/efficiency.py tests/test_efficiency.py
git commit -m "$(cat <<'EOF'
feat: efficiency metrics (relative loss delta, collision loss, mechanical COT)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: R=0.3 アンカー + 連続性ゲート（テストのみ）

**Files:**
- Test: `tests/test_wheel_limit_anchor.py`

R 掃引を Phase 3.5 検証済みサイクルに固定する最強ゲート。新規プロダクションコードは不要
（既存 `find_limit_cycle` + `rocker_compass` を使う）。

- [ ] **Step 1: Write the test**

```python
# tests/test_wheel_limit_anchor.py
import numpy as np

from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.search import find_limit_cycle
from crane import references_mcgeer as ref


def _fp_at_R(R, guess):
    p = RockerCompassParams(
        m=ref.M_LEG, m_h=ref.M_HIP, c=ref.C_HIP_TO_COM, rho=ref.RHO_GYR,
        L=ref.L_LEG, R=R, gamma=ref.GAMMA_GAIT, g=ref.G,
    )
    return find_limit_cycle(make_rocker_compass(p), np.array(guess))


def test_R030_reproduces_phase35_cycle():
    """最強アンカー: R=0.3 で Phase 3.5 検証済みサイクルを再現。"""
    fp = _fp_at_R(0.3, ref.SECTION_GUESS)
    assert fp.converged
    np.testing.assert_allclose(
        fp.y, [0.30844, -1.26256, -0.87914], atol=1e-3
    )
    assert np.max(np.abs(fp.eigenvalues)) < 1.0  # 安定


def test_continuation_connects_neighboring_R():
    """連続性: R=0.3 の解を seed に R=0.25 が近傍へ連続接続する（continuation の健全性）。"""
    fp030 = _fp_at_R(0.3, ref.SECTION_GUESS)
    assert fp030.converged
    fp025 = _fp_at_R(0.25, fp030.y)  # 前の解を seed
    assert fp025.converged
    # 小さな ΔR では不動点も小さく動く（branch hopping していない）
    assert np.linalg.norm(fp025.y - fp030.y) < 0.3
```

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/test_wheel_limit_anchor.py -v`
Expected: PASS (2 tests). これはゲートであり、FAIL する場合は R=0.3 の値が Phase 3.5 と
食い違う（重大）か continuation seed が効いていない。緩和せず原因を調査して報告すること。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format tests/test_wheel_limit_anchor.py && uv run ruff check tests/test_wheel_limit_anchor.py
git add tests/test_wheel_limit_anchor.py
git commit -m "$(cat <<'EOF'
test: R=0.3 anchor to phase 3.5 cycle + continuation continuity gate

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: 二段掃引スクリプト `wheel_limit.py`

**Files:**
- Create: `scripts/wheel_limit.py`

- [ ] **Step 1: Write the script**

```python
# scripts/wheel_limit.py
"""Phase 4a.1: rocker_compass の R を掃引し「車輪への漸近」を検証。

固定 γ=0.030 で R∈(0,L) を continuation 追跡。各 R で δ(R)・max|λ|(R) を測り、
粗い R で basin_fraction(R) を測る。曲線・JSON・basin montage を data/runs/ に出力。
Usage: uv run python scripts/wheel_limit.py [--dr 0.05] [--basin-res 50] [--workers N]
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from crane import references_mcgeer as ref  # noqa: E402
from crane.basin import basin_slice  # noqa: E402
from crane.efficiency import mechanical_cot, relative_loss, step_collision_loss  # noqa: E402
from crane.models.rocker_compass import (  # noqa: E402
    RockerCompassParams,
    kinetic_energy,
    make_rocker_compass,
)
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle  # noqa: E402
from crane.stride import StrideError, stride  # noqa: E402


def _params(R):
    return RockerCompassParams(
        m=ref.M_LEG, m_h=ref.M_HIP, c=ref.C_HIP_TO_COM, rho=ref.RHO_GYR,
        L=ref.L_LEG, R=R, gamma=ref.GAMMA_GAIT, g=ref.G,
    )


def _hip_x(theta_st, R, L):
    """stance 接触フレームでの hip 水平位置: -R*theta - (L-R)*sin(theta)。"""
    return -R * theta_st - (L - R) * np.sin(theta_st)


def _measure(R, guess):
    """R で limit cycle を求め、(fp, delta, cot, lam_max) を返す。失敗時 None。"""
    p = _params(R)
    model = make_rocker_compass(p)
    fp = find_limit_cycle(model, np.array(guess))
    if not fp.converged:
        return None
    lam_max = float(np.max(np.abs(fp.eigenvalues)))
    try:
        result = stride(model, model.lift(fp.y))
    except StrideError:
        return None
    loss, ke_pre, ke_post = step_collision_loss(
        result.x_strike, result.x_end, lambda x: kinetic_energy(x, p)
    )
    delta = relative_loss(ke_pre, ke_post)
    # 一歩水平距離 ≈ stance 接触フレームでの hip 前進（report-only の COT 用）
    theta_start = float(result.x[0, 0])
    theta_strike = float(result.x_strike[0])
    step_len = abs(_hip_x(theta_strike, R, p.L) - _hip_x(theta_start, R, p.L))
    cot = mechanical_cot(loss, m=p.m, g=p.g, step_length=step_len) if step_len > 1e-9 else None
    return fp, delta, cot, lam_max


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dr", type=float, default=0.05)
    parser.add_argument("--basin-res", type=int, default=50)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    # --- 細かい掃引（R=0.3 起点に上下） ---
    records = {}  # R -> dict
    # 0.3 起点
    base = _measure(0.3, ref.SECTION_GUESS)
    if base is None:
        raise SystemExit("ERROR: R=0.3 limit cycle not found")
    records[0.3] = base

    def sweep(direction):
        guess = records[0.3][0].y
        R = 0.3
        while True:
            R = round(R + direction * args.dr, 4)
            if R <= 0.0 or R >= ref.L_LEG:
                break
            m = _measure(R, guess)
            if m is None or m[3] >= 1.0:  # 収束失敗 or 不安定 = 歩行消失点
                records[R] = None
                break
            records[R] = m
            guess = m[0].y

    sweep(+1)  # 車輪側
    sweep(-1)  # 点足側

    Rs = sorted(r for r, v in records.items() if v is not None)
    delta_curve = [records[r][1] for r in Rs]
    cot_curve = [records[r][2] for r in Rs]
    lam_curve = [records[r][3] for r in Rs]
    print("R sweep:", [(r, round(records[r][1], 4), round(records[r][3], 3)) for r in Rs])

    # --- 粗い掃引で basin ---
    basin_Rs = [r for r in (0.05, 0.2, 0.4, 0.6, 0.8) if r in records and records[r] is not None]
    if Rs and Rs[-1] not in basin_Rs:
        basin_Rs.append(Rs[-1])  # 消失直前も
    basin_fracs = {}
    basin_grids = {}
    for r in basin_Rs:
        fp = records[r][0]
        res = basin_slice(
            make_rocker_compass, _params(r), fp.y,
            axes=(0, 1), half_widths=(0.16, 0.65), resolution=args.basin_res,
            model_name=f"rocker_compass_R{r}", n_workers=args.workers,
        )
        basin_fracs[r] = res.basin_fraction
        basin_grids[r] = res
        print(f"basin R={r}: fraction={res.basin_fraction:.3f}")

    # --- 曲線プロット ---
    run_dir = new_run_dir("wheel_limit")
    fig, axs = plt.subplots(3, 1, figsize=(7, 10), sharex=True)
    axs[0].plot(Rs, delta_curve, "o-")
    axs[0].set_ylabel("delta (relative loss)")
    axs[0].axhline(0.0, color="gray", lw=0.5)
    axs[1].plot(Rs, lam_curve, "o-")
    axs[1].axhline(1.0, color="red", lw=0.5)
    axs[1].set_ylabel("max|lambda|")
    if basin_fracs:
        axs[2].plot(list(basin_fracs), list(basin_fracs.values()), "s-")
    axs[2].set_ylabel("basin_fraction")
    axs[2].set_xlabel("R (foot radius); R->L=1 is the wheel limit")
    fig.suptitle("Approaching the wheel: rocker_compass R-continuation (gamma=0.030)")
    fig.savefig(run_dir / "R_sweep.png", dpi=150)
    plt.close(fig)

    # --- basin montage ---
    if basin_grids:
        cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])
        order = sorted(basin_grids)
        fig2, axs2 = plt.subplots(1, len(order), figsize=(4 * len(order), 4))
        if len(order) == 1:
            axs2 = [axs2]
        for ax, r in zip(axs2, order):
            g = basin_grids[r]
            ax.imshow(
                g.grid, origin="lower",
                extent=[g.ax0_vals[0], g.ax0_vals[-1], g.ax1_vals[0], g.ax1_vals[-1]],
                aspect="auto", cmap=cmap, vmin=0, vmax=2, interpolation="nearest",
            )
            ax.plot(g.fixed_point[0], g.fixed_point[1], "k*", markersize=10)
            ax.set_title(f"R={r}\nfrac={g.basin_fraction:.3f}")
        fig2.tight_layout()
        fig2.savefig(run_dir / "basin_R_montage.png", dpi=150)
        plt.close(fig2)

    data = {
        "gamma": ref.GAMMA_GAIT,
        "sin_gamma": float(np.sin(ref.GAMMA_GAIT)),
        "R": Rs,
        "delta": delta_curve,
        "cot": cot_curve,
        "max_lambda": lam_curve,
        "fixed_points": [records[r][0].y.tolist() for r in Rs],
        "basin_fraction": basin_fracs,
        "walking_vanishes_at": (
            min((r for r, v in records.items() if v is None), default=None)
        ),
    }
    (run_dir / "R_sweep.json").write_text(json.dumps(data, indent=2, default=str))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 低コストのスモーク実行**

Run: `uv run python scripts/wheel_limit.py --dr 0.1 --basin-res 16 --workers 4`
Expected: `R sweep: [...]` が表示され、δ・max|λ| のリストが出る。`basin R=...` が数点表示。
`R_sweep.png` / `basin_R_montage.png` / `R_sweep.json` が `data/runs/.../wheel_limit/` に生成。
report: 実測の δ(R) と basin_fraction(R) の傾向（R↑ で δ↓ か、basin↑ か）と、歩行消失点
（`walking_vanishes_at`）を**正直に**。仮説に合わせ込まない。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format scripts/wheel_limit.py && uv run ruff check scripts/wheel_limit.py
git add scripts/wheel_limit.py
git commit -m "$(cat <<'EOF'
feat: wheel_limit.py R-continuation sweep (delta, stability, basin curves)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: rimless wheel の provenance（McGeer 1990a）

**Files:**
- Modify: `src/crane/references_mcgeer.py`

設計ゲート 4: 「rimless wheel = 最簡受動歩行モデル」を一次資料で裏取りし provenance 記録。

- [ ] **Step 1: 一次資料を確認して記録**

WebSearch/WebFetch で **McGeer 1990a "Passive Dynamic Walking" (IJRR 9(2):62-82)** を当たり、
論文が straight-leg walker の導入として **rimless wheel（"synthetic wheel" / spoked wheel）** を
出発点に使っているか確認する。McGeer は §2 付近で rimless wheel を passive walking の
最簡モデルとして導入していることが多い。実際に読めた記述を provenance 付きで
`src/crane/references_mcgeer.py` の末尾に追記:

```python
# --- rimless wheel = 最簡受動歩行モデル（issue #19 / Phase 4a.1 の概念的支柱）---
# McGeer, T. (1990) "Passive Dynamic Walking", IJRR 9(2):62-82。
# URL: <実際に参照した URL>  取得日: 2026-06-16
# 記録する事実（原典の記述に基づく。読めた内容のみ。読めなければ None + その旨）:
RIMLESS_WHEEL_PROVENANCE: str | None = "<原典で確認した記述の要約、または取得失敗の明記>"
```

**記憶からの断定はしない。** McGeer 1990a の本文が読めなければ、その旨を文字列に明記し
（例: "McGeer 1990a 本文未取得; rimless wheel の位置づけは要再確認"）、`RIMLESS_WHEEL_PROVENANCE`
を honest に書く。捏造禁止。可能なら Coleman & Ruina / Garcia 1998 等の二次資料で補強し
出典を明記する。

- [ ] **Step 2: ruff then commit**

```bash
uv run ruff format src/crane/references_mcgeer.py && uv run ruff check src/crane/references_mcgeer.py
git add src/crane/references_mcgeer.py
git commit -m "$(cat <<'EOF'
docs: record rimless-wheel provenance for wheel-limit framing (McGeer 1990a)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: 本番実行・記録・目視判定ゲート

**Files:**
- Modify: `README.md`
- Modify: `GOALS.md`

- [ ] **Step 1: 本番解像度で掃引を実行**

Run: `uv run python scripts/wheel_limit.py --dr 0.05 --basin-res 50 --workers <cores>`
（数分〜30分）。`R_sweep.png` / `basin_R_montage.png` / `R_sweep.json` を得る。
basin 窓 (0.16,0.65) が全 R で basin を切り落としていないか確認（極端な R で収束領域が
窓端に接していたら half_widths を広げて再実行し、その旨を記録）。

- [ ] **Step 2: 全テスト + エネルギー収支チェックを確認**

Run: `uv run pytest -q && uv run ruff format --check . && uv run ruff check .`
Expected: 全 green、ruff クリーン。
加えて `R_sweep.json` の `cot` 各値が `sin_gamma`（≈0.0300）に近いか目視確認
（機械的 COT ≈ sin γ の内部整合チェック。大きくずれる場合は step_length 幾何の近似誤差
として README に明記）。

- [ ] **Step 3: README に Phase 4a.1 セクションを追加**

`README.md` の Phase 4a の後に「Phase 4a.1 の結果」を追加。**Step 1 の実測値で埋める**:
δ(R) の傾向（R↑ で減少か）、basin_fraction(R) の傾向、歩行消失点 R、max|λ|(R)、
機械的 COT≈sin γ の確認結果、「コンパスは車輪への漸近」仮説 H1/H2/H3 の正直な判定。
成果物パス（gitignore のローカル）も明記。

- [ ] **Step 4: GOALS に Phase 4a.1 セクションを追加**

`GOALS.md` の Phase 4a の後に「Phase 4a.1: 車輪への漸近（R 連続軸）」を追加。ゲート:
```markdown
- [x] efficiency.py（δ・衝突損失・COT）+ R=0.3 アンカー + 連続性ゲート（テスト pass）
- [x] R 掃引: δ(R)・max|λ|(R)・basin_fraction(R) を測定（実測値を記録）
- [x] エネルギー収支チェック: 機械的 COT(R) ≈ sin γ
- [x] rimless wheel provenance（McGeer 1990a、取得可否を honest に記録）
- [ ] **ぽんぽこ殿の目視判定**: R_sweep.png が「車輪への漸近」（δ↓・basin↑）を支持して見える
```
H1/H2/H3 の判定を正直に記録（仮説に反する結果も finding として明記）。
最後の目視判定は UNCHECKED で残す（コントローラが依頼）。

- [ ] **Step 5: Commit**

```bash
git add README.md GOALS.md
git commit -m "$(cat <<'EOF'
docs: record phase 4a.1 wheel-limit results and visual judgment gate

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6: 目視判定を依頼**

`R_sweep.png` をコントローラ経由でぽんぽこ殿に提示し、「δ↓・basin↑ が車輪への漸近を
支持して見えるか」の判定を仰ぐ。合格後に GOALS のチェックを埋め、PR へ。

---

## Self-Review チェック結果

- **Spec coverage:** 設計 A〜F を網羅。A スコープ→全タスク、B 二段掃引→Task3、C 指標→Task1+Task3、
  D ゲート（R=0.3 アンカー→Task2、R→0 連続性→Task2、COT≈sinγ→Task5 Step2、provenance→Task4、
  目視→Task5）、E モジュール/出力→Task1/Task3、F 命名→Task5。
- **Placeholder scan:** McGeer 1990a の provenance のみ「実装時に取得」とした（CLAUDE.md 準拠、
  取得不能時の honest フォールバックを明記、捏造禁止）。他にプレースホルダなし。
- **Type consistency:** `relative_loss(ke_pre, ke_post)`、`step_collision_loss(x_strike, x_end,
  kinetic_energy) -> (loss, ke_pre, ke_post)`、`mechanical_cot(loss, *, m, g, step_length)`、
  `_params(R)` / `_measure(R, guess)` / `_hip_x(theta, R, L)`、basin_slice の引数、
  rocker_compass の (project=[x0,x2,x3], 4D state) が全タスクで一貫。
- **既知の近似（明示済み）:** 機械的 COT の step_length は stance 接触フレームの hip 前進近似で
  report-only。脚ラベル交換のフレーム遷移を厳密には追わないため、COT は内部チェック（≈sinγ を
  期待）に留め、ハードゲートにしない。主指標は δ(R)。
