# Crane

受動歩行機械のリミットサイクルを **hybrid dynamics + Poincaré shooting** で
数学的に発見する研究プロジェクト。Heron の後継。

物理エンジン（Genesis / MuJoCo）でのランダム探索では basin が薄すぎて
リミットサイクルを狙い撃ちできないことが Heron で判明した。
Crane はアプローチを逆転させ、先に EOM・衝突写像を解析的に立てて不動点を求め、
その後でシミュレーションに渡す。

## Phase 1 の結果 (Garcia 1998 Simplest Walker, γ = 0.009)

| 種別 | θ* | θ̇* | max\|λ\| | 安定性 |
|------|-----|-----|---------|-------|
| long-period | 0.2003109 | −0.1998325 | 0.589 | 安定 |
| short-period | 0.1939374 | −0.2038669 | 4.00 | 不安定 |

文献 O(γ) 近似値 (0.199529, −0.198983) との差: 成分あたり ~8.5e-4 < 許容 1.5e-3（O(γ^(5/3)) 打ち切り誤差の範囲内）

**θ* ∝ γ^(1/3) スケーリング則**: γ ∈ {0.004, 0.006, 0.009, 0.012} で log-log slope = 0.3287（理論値 1/3、誤差 0.005）

デモ実行 (1% 摂動から 30 歩、収束):

```
stride 0:  deviation 1.8e-3
stride 29: deviation 6.6e-10   ← リミットサイクル収束
```

## Phase 2 の結果 (Goswami 1996 Compass Gait, γ = 3°)

| 種別 | θ_st* | θ̇_st* | θ̇_sw* | max\|λ\| | 安定性 |
|------|--------|--------|--------|---------|-------|
| period-1 (γ=3°) | 0.27103 | −1.09238 | −0.37737 | 0.580 | 安定 |

文献値 (Goswami INRIA RR-2996) との一致: T=0.7343s（文献 0.735s、0.1%）、strike 角 0.2710（文献 0.271）

**固有値**: 印刷値 0.332 は2歩合成写像の multiplier。1歩写像では 0.580²=0.336（1.3% 一致）

**Period-doubling カスケード再現**:
- period-1 不安定化: γ=4.40°（文献 4.37°、分解能 0.05° 以内）
- period-2 branch: 4.40°–5.95° 追跡、period-2 不安定化 ~4.95°（文献の第2分岐 4.9°）
- flip bifurcation（実固有値が −1 通過）確認

デモ実行 (0.5% 摂動から 30 歩、収束):

```
stride 0:  deviation 6.7e-3
stride 29: deviation 5.3e-10   ← リミットサイクル収束
```

## Phase 3 の結果 (Hsu Chen 2007 点足 Kneed Walker, γ = 0.0504 rad)

4相 hybrid（unlocked 3-DOF → knee-strike → locked 2-DOF → heel-strike）。
断面は heel-strike 直後の 3D `y = (θ_st, θ̇_st, θ̇_sw)`。

| 量 | 本実装 | 文献 (MIT thesis Ch.3) | 差 |
|------|--------|------------------------|-----|
| q1 (stance 角) | 0.1882 | 0.1877 | 0.27% |
| q̇1 | −1.1096 | −1.1014 | 0.7% |
| step period | 0.5643 s | 0.559（図読取 ±2%）| 1.0% |
| heel-strike 直前 stance 角 | −0.2386 | −0.2380 | 0.25% |
| 安定性 max\|λ\| | **0.655 < 1（安定）** | locally stable と明言 | — |

Newton shooting が残差 **5e-14**（4反復）で不動点 y*=(0.23859, −1.10959, −0.05715) に収束。
固有値は複素対 0.353±0.552i と実 0.144。膝屈曲は毎 stride ~30°。

**compass 退化ゲート**: shank 質量 → 0 で locked 力学・heel-strike・全 stride 不動点・上位固有値が
Phase 2 検証済み compass (0.27103, −1.09238, −0.37737) に一致 — 記号導出レイヤーの最強検証。

デモ実行 (0.5% 摂動から 12 歩、収束):

```
stride 0:  deviation 9.7e-3  knee_flexion_max 30.0deg
stride 11: deviation 1.6e-4  knee_flexion_max 29.9deg   ← リミットサイクル収束
```

## セットアップ

```bash
uv sync
```

Python 3.12 固定。

## 使い方

```bash
# Phase 1 デモ歩行: Simplest Walker（walk.mp4 + phase_portrait.png + meta.json を data/runs/ に出力）
uv run python scripts/walk_simplest.py [--gamma 0.009] [--strides 30] [--perturb 0.01]

# Phase 2 デモ歩行: Compass Gait
uv run python scripts/walk_compass.py [--gamma-deg 3.0] [--strides 30] [--perturb 0.005]

# Phase 2 分岐図: slope continuation + period-doubling（bifurcation.csv + bifurcation.png）
uv run python scripts/bifurcation_compass.py [--gamma-max-deg 6.0] [--step-deg 0.05]

# Phase 3 デモ歩行: 点足 Kneed Walker（4セグメント walk.mp4 + phase_portrait.png + meta.json）
uv run python scripts/walk_kneed.py [--gamma-deg <published>] [--strides 30] [--perturb 0.005]

# テスト（63 tests）
uv run pytest
```

## ドキュメント

- 設計: `docs/2026-06-12-crane-design.md`
- フェーズゲート: `GOALS.md`

## Heron との関係

Heron (Genesis/MuJoCo ポート) で確認されたこと: 物理エンジン上でのランダム探索は
basin が薄すぎてリミットサイクル候補に到達できない。
Crane はその教訓から出発し、解析的にサイクルを求めてからシミュレーションに渡す逆順戦略を採る。
