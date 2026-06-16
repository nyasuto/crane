# Phase 5c: 能動 vs 受動 basin 比較 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** γ=0.030 で R×push_off 格子の能動 basin_fraction(R,P) を測り、受動（P=0）と比較して「制御は受動順位を保存するか／形態差を平坦化するか」を検証する。

**Architecture:** 新規物理なし。`basin_slice`（Phase 4a）に `make_powered_rocker_compass`（Phase 5b）+ `PoweredRockerCompassParams` を渡して掃引するだけ。受動不動点は R-continuation（Phase 4a.1 と同手法）で各 R に供給し、能動不動点はそれを seed に直接収束。

**Tech Stack:** numpy, scipy, matplotlib, multiprocessing（すべて既存）。新依存なし。

**設計ドキュメント:** `docs/2026-06-17-active-basin-design.md`

---

## 既存コードの要点
- `from crane.basin import basin_slice` → `basin_slice(make_fn, params, fixed_point, *, axes=(0,1), half_widths, resolution, model_name, max_strides=20, converge_tol=1e-3, n_workers=None) -> BasinResult(.grid, .basin_fraction, .ax0_vals, .ax1_vals, .fixed_point)`。
- `from crane.models.powered_rocker_compass import PoweredRockerCompassParams, make_powered_rocker_compass`（Phase 5b、push_off=0 で受動 rocker に厳密退化）。
- `from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass`。
- `from crane.search import find_limit_cycle`。
- nominal: m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, γ=0.030。Phase 3.5 受動不動点 @ R=0.3: (0.30844, −1.26256, −0.87914)。
- **Phase 4a.1 受動 basin(R) @ γ=0.030, 窓(0.16,0.65), res50**: R=0.05→0.037, 0.20→0.072, 0.40→0.200, 0.60→0.430。

## File Structure
- Create: `scripts/active_basin_sweep.py` — (R,P) 掃引、basin_slice、曲線＋montage＋JSON。
- Create: `tests/test_active_basin_gate.py` — push_off→0 で能動 basin が受動 basin に一致（内部ゲート）。
- Modify: `README.md`, `GOALS.md` — Phase 5c セクション。

---

## Task 1: push_off→0 内部ゲート

**Files:**
- Test: `tests/test_active_basin_gate.py`

powered(P=0) の basin が受動 rocker_compass の basin に一致することを低解像度で確認。検証済み 4a.1 に錨。

- [ ] **Step 1: Write the test**

```python
# tests/test_active_basin_gate.py
import numpy as np

from crane.basin import basin_slice
from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.search import find_limit_cycle

NOM = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)
GUESS = np.array([0.30, -1.0, -0.40])


def test_powered_pushoff_zero_basin_matches_passive():
    """P=0 の能動 basin が受動 rocker_compass の basin に一致（内部ゲート、低解像度）。"""
    passive_params = RockerCompassParams(**NOM)
    powered_params = PoweredRockerCompassParams(**NOM, push_off=0.0)

    fp_passive = find_limit_cycle(make_rocker_compass(passive_params), GUESS)
    fp_powered = find_limit_cycle(make_powered_rocker_compass(powered_params), GUESS)
    assert fp_passive.converged and fp_powered.converged
    np.testing.assert_allclose(fp_powered.y, fp_passive.y, atol=1e-8)

    kw = dict(axes=(0, 1), half_widths=(0.16, 0.65), resolution=15, n_workers=1)
    b_passive = basin_slice(
        make_rocker_compass, passive_params, fp_passive.y, model_name="passive", **kw
    )
    b_powered = basin_slice(
        make_powered_rocker_compass, powered_params, fp_powered.y, model_name="powered", **kw
    )
    assert np.array_equal(b_powered.grid, b_passive.grid)
    assert b_powered.basin_fraction == b_passive.basin_fraction
```

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/test_active_basin_gate.py -v`
Expected: PASS。powered(P=0) は受動モデルと厳密同一なので basin grid も一致。
FAIL 時は push-off 退化のバグ。**緩和せず**調査。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format tests/test_active_basin_gate.py
uv run ruff check tests/test_active_basin_gate.py
git add tests/test_active_basin_gate.py
git commit -m "$(cat <<'EOF'
test: active basin push-off->0 gate (matches passive rocker basin)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: active_basin_sweep.py（(R,P) 掃引）

**Files:**
- Create: `scripts/active_basin_sweep.py`

- [ ] **Step 1: Write the script**

```python
# scripts/active_basin_sweep.py
"""Phase 5c: γ=0.030 で R×push_off 格子の能動 basin を測り受動と比較。

受動不動点は R=0.3 から R-continuation で各 R に供給し、能動不動点はそれを seed に直接収束。
basin_fraction(R,P) の曲線・montage・JSON（順位・変動係数つき）を data/runs/ に出力。
Usage: uv run python scripts/active_basin_sweep.py [--basin-res 50] [--workers N]
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from crane.basin import basin_slice  # noqa: E402
from crane.models.powered_rocker_compass import (  # noqa: E402
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass  # noqa: E402
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle  # noqa: E402

NOM = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, gamma=0.030, g=9.81)
R_VALUES = [0.05, 0.20, 0.40, 0.60]
P_VALUES = [0.0, 0.04, 0.08]


def passive_fp(R, n_steps=12):
    """R=0.3 から R-continuation で R の受動不動点を返す。失敗時 None。"""
    guess = np.array([0.30844, -1.26256, -0.87914])
    for r in np.linspace(0.3, R, n_steps):
        fp = find_limit_cycle(make_rocker_compass(RockerCompassParams(**NOM, R=float(r))), guess)
        if not fp.converged:
            return None
        guess = fp.y
    return guess


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basin-res", type=int, default=50)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    fracs = {}  # (R,P) -> basin_fraction
    grids = {}  # (R,P) -> BasinResult
    for R in R_VALUES:
        seed = passive_fp(R)
        if seed is None:
            print(f"R={R}: passive fp not found, skip")
            continue
        for P in P_VALUES:
            params = PoweredRockerCompassParams(**NOM, R=R, push_off=P)
            fp = find_limit_cycle(make_powered_rocker_compass(params), seed)
            if not fp.converged:
                print(f"R={R} P={P}: cycle not found, skip")
                continue
            res = basin_slice(
                make_powered_rocker_compass, params, fp.y,
                axes=(0, 1), half_widths=(0.16, 0.65), resolution=args.basin_res,
                model_name=f"R{R}_P{P}", n_workers=args.workers,
            )
            fracs[(R, P)] = res.basin_fraction
            grids[(R, P)] = res
            seed = fp.y  # 次の P を前の能動不動点で seed
            print(f"R={R} P={P}: basin_fraction={res.basin_fraction:.3f}")

    run_dir = new_run_dir("active_basin")

    # 曲線: basin vs R, P ごとに 1 本
    fig, ax = plt.subplots(figsize=(7, 5))
    for P in P_VALUES:
        xs = [R for R in R_VALUES if (R, P) in fracs]
        ys = [fracs[(R, P)] for R in xs]
        ax.plot(xs, ys, "o-", label=f"push_off={P}")
    ax.set_xlabel("R (foot radius)")
    ax.set_ylabel("basin_fraction")
    ax.set_title("Active vs passive basin: basin(R) for each push-off (gamma=0.030)")
    ax.legend()
    fig.savefig(run_dir / "active_basin_curves.png", dpi=150)
    plt.close(fig)

    # montage: 各 (R,P) の basin
    cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])
    keys = [(R, P) for R in R_VALUES for P in P_VALUES if (R, P) in grids]
    if keys:
        nP = len(P_VALUES)
        nR = len([R for R in R_VALUES if any((R, P) in grids for P in P_VALUES)])
        fig2, axs = plt.subplots(nR, nP, figsize=(3 * nP, 3 * nR))
        axs = np.atleast_2d(axs)
        rows = [R for R in R_VALUES if any((R, P) in grids for P in P_VALUES)]
        for i, R in enumerate(rows):
            for j, P in enumerate(P_VALUES):
                a = axs[i, j]
                if (R, P) not in grids:
                    a.axis("off")
                    continue
                g = grids[(R, P)]
                a.imshow(
                    g.grid, origin="lower",
                    extent=[g.ax0_vals[0], g.ax0_vals[-1], g.ax1_vals[0], g.ax1_vals[-1]],
                    aspect="auto", cmap=cmap, vmin=0, vmax=2, interpolation="nearest",
                )
                a.set_title(f"R={R} P={P}\nfrac={g.basin_fraction:.3f}", fontsize=8)
        fig2.tight_layout()
        fig2.savefig(run_dir / "active_basin_montage.png", dpi=150)
        plt.close(fig2)

    # 仮説指標: 順位（各 P で argmax_R）と変動係数（R 間の std/mean）
    summary = {}
    for P in P_VALUES:
        ys = {R: fracs[(R, P)] for R in R_VALUES if (R, P) in fracs}
        if not ys:
            continue
        vals = np.array(list(ys.values()))
        argmax_R = max(ys, key=ys.get)
        cv = float(np.std(vals) / np.mean(vals)) if np.mean(vals) > 0 else None
        summary[P] = {"argmax_R": argmax_R, "cv_across_R": cv, "fractions": ys}

    data = {
        "gamma": 0.030,
        "R_values": R_VALUES,
        "P_values": P_VALUES,
        "basin_fraction": {f"R{R}_P{P}": fracs[(R, P)] for (R, P) in fracs},
        "summary_by_pushoff": {str(P): summary[P] for P in summary},
    }
    (run_dir / "active_basin.json").write_text(json.dumps(data, indent=2, default=str))
    print("summary:", {P: (summary[P]["argmax_R"], round(summary[P]["cv_across_R"], 3))
                        for P in summary if summary[P]["cv_across_R"] is not None})
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 低解像度スモーク実行**

Run: `uv run python scripts/active_basin_sweep.py --basin-res 16 --workers 4`
Expected: 各 (R,P) の `basin_fraction` が表示され、`active_basin_curves.png`/`active_basin_montage.png`/
`active_basin.json` が生成。`summary`（各 P の argmax_R と変動係数）が出る。実測の傾向を報告。
内部整合チェック: P=0 行の basin_fraction が 4a.1 受動値（R=0.05→0.037, 0.2→0.072, 0.4→0.200, 0.6→0.430）
の近く（低解像度なのでズレは許容、傾向一致を見る）。

- [ ] **Step 3: ruff then commit**

```bash
uv run ruff format scripts/active_basin_sweep.py && uv run ruff check scripts/active_basin_sweep.py
git add scripts/active_basin_sweep.py
git commit -m "$(cat <<'EOF'
feat: active_basin_sweep.py (basin(R,P) active vs passive at gamma=0.030)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: 本番実行・記録・歩容判定

**Files:**
- Modify: `README.md`, `GOALS.md`

- [ ] **Step 1: 本番解像度で実行**

Run: `uv run python scripts/active_basin_sweep.py --basin-res 50 --workers <cores>`
（数十分）。`active_basin_curves.png`/`active_basin_montage.png`/`active_basin.json` を得る。
P=0 行が 4a.1 受動 basin(R)（0.037/0.072/0.200/0.430）を再現するか内部整合確認（窓・解像度同一）。
基準と大きくズレたら原因（窓・seed・解像度）を調査。

- [ ] **Step 2: 全テスト + ruff**

Run: `uv run pytest -q && uv run ruff format --check . && uv run ruff check .`
Expected: 全 green、ruff クリーン。

- [ ] **Step 3: README に Phase 5c セクションを追加**

`README.md` の Phase 5b の後に「Phase 5c の結果（能動 vs 受動 basin 比較）」を追加。**実測値で埋める**:
basin_fraction(R,P) 表、各 P の argmax_R（順位保存か）、R 間変動係数の P 依存（平坦化か）、
仮説の正直な判定（受動順位は保存されるか／制御は形態差を消すか）。成果物パス・gitignore を明記。

- [ ] **Step 4: GOALS に Phase 5c セクションを追加**

`GOALS.md` の Phase 5b の後に「Phase 5c: 能動 vs 受動 basin 比較 — 完了 (2026-06-17)」を追加。ゲート:
```markdown
- [x] push_off→0 内部ゲート: P=0 の能動 basin が受動 rocker_compass の basin に一致
- [x] (R,P) 格子で能動 basin_fraction を測定（実測表を記録）
- [x] 仮説判定: 順位保存／平坦化を正直に報告（実測値）
- [ ] **ぽんぽこ殿の目視判定**: active_basin_curves.png / montage が主張を支持して見える
```
H（順位保存）/（平坦化）の判定を実測通り記録。能動歩行アーク（5a/5b/5c）の締めくくりと位置づけ。
目視判定は UNCHECKED で残す。

- [ ] **Step 5: Commit**

```bash
git add README.md GOALS.md
git commit -m "$(cat <<'EOF'
docs: record phase 5c active vs passive basin results and judgment gate

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6: 目視判定を依頼**

`active_basin_curves.png` / `active_basin_montage.png` をコントローラ経由でぽんぽこ殿に提示し判定を仰ぐ。
合格後に GOALS のチェックを埋め PR へ。

---

## Self-Review チェック結果
- **Spec coverage:** 設計 A〜F を網羅。B 方法→Task2、C 指標（順位・変動係数）→Task2、
  D ゲート（push_off→0 内部→Task1、正直報告→Task3、目視→Task3）、E ファイル→Task1/2/3、F 命名→Task3。
- **Placeholder scan:** なし。実測値で埋める箇所は Task3 で明示。
- **Type consistency:** `basin_slice(make_fn, params, fixed_point, *, axes, half_widths, resolution, model_name, n_workers)`、
  `PoweredRockerCompassParams(**NOM, R, push_off)`、`RockerCompassParams(**NOM, R)`、`passive_fp(R)`、
  `fracs[(R,P)]` / `grids[(R,P)]` が一貫。NOM に R は含めず各所で R を渡す（NOM は gamma 等のみ）。
- **既知の留意点:** 受動 fp は R-continuation（R=0.3 起点）で供給。能動 fp は受動 fp を seed に直接収束
  （push-off は摂動）。非収束 (R,P) は skip して記録（finding）。
