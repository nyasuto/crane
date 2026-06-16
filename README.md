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

## Phase 3.5 の結果 (McGeer 1990a Rocker-foot Compass, γ = 0.030, R = 0.3)

円弧足（半径 R の round foot）2D compass。剛体脚を 2点質量（m/2 を hip 距離 c±ρ に配置）で
表現し、既存の点質量 derive レイヤーをそのまま再利用。stance 円弧足は転がり接触
（曲率中心 C=[−R·θ_st, R]）。R→0 で点足 compass に退化。

| 量 | 本実装 | 文献 (McGeer IJRR 9(2), §5) | 差 |
|------|--------|-----------------------------|-----|
| q* = (θ_st, θ̇_st, θ̇_sw) | (0.30844, −1.26256, −0.87914) | — | — |
| step period | 0.8417 s | 2.5·√(l/g) = 0.7982 s | 5.45%（±10% 内）|
| strike 時 half-interleg 角 \|θ_st\| | 0.308 | α₀ = 0.30 | 傍証（ゲート外）|
| 安定性 max\|λ\| | **0.4316 < 1（安定）** | Fig.9（R=0.3 で安定化）と整合 | — |

固有値は 0.4316 と複素対 −0.1529±0.2986j（大きさ 0.4316, 0.3355, 0.3355）。
McGeer §5 straight-leg モデルには印刷数値固有値がない（Table 1 は kneed テスト機械）ため、
固有値ゲートは reference-only、ゲートは存在性・安定性・step period・R→0 退化。

**R→0 退化ゲート（最強検証）**: R→1e-9, ρ→1e-9, c=b で全 stride 不動点・固有値が
Phase 2 検証済み点足 compass (0.27103, −1.09238, −0.37737) に一致（不動点 ~1e-10 / 固有値 ~1e-8）。

**R-continuation**: gait family を R=0.30→0.24→0.18→0.12 まで追跡、全て収束・安定。

デモ実行 (摂動から 30 歩、収束):

```
stride 0:  deviation 3.95e-3
stride 29: deviation 8.5e-14   ← リミットサイクル収束
```

## Phase 3.6 の結果 (McGeer 1990b Rocker-foot Kneed Walker)

Phase 3 点足 kneed（kneed.py、4相 hybrid）に Phase 3.5 の円弧足（半径 R の round foot
転がり接触）を合流させた McGeer 1990b 原機械の**点質量再現**。`models/rocker_kneed.py` は
kneed.py の _build() に 5 箇所のみの外科的変更を施したもの（derive レイヤー / kneed.py は不変）。

| 量 | 本実装 | 文献 (McGeer ICRA 1990b) | 差 |
|------|--------|--------------------------|-----|
| y* = (θ_st, θ̇_st, θ̇_sw) | (0.35661, −1.20906, −0.54601) | — | — |
| step period | 0.7945 s | 2.7·√(l/g) = 0.862 s | 7.84%（±15% 内）|
| 膝屈曲 | ≈ 27.6° | — | — |
| 固有値（大きさ）| 0.616, 0.616, 0.046 | 印刷値 0.447（4D 断面）| reference-only |
| 安定性 max\|λ\| | **0.616 < 1（安定）** | 安定サイクルと明言 | — |

不動点 y* は references SECTION_GUESS から Newton 直接収束。固有値は複素対
0.39527±0.47248i と実 0.04579。固有値の差は分布質量（McGeer 実機）vs 点質量（本実装）+
断面次元差（4D vs 3D）に由来する documented な reference check で、ゲートは
存在性・安定性・step period・R→0 退化。

**R→0 退化ゲート（最強検証）**: R→1e-9 で全 stride 不動点・固有値が
Phase 3 検証済み点足 kneed (0.23859, −1.10959, −0.05715) に一致（不動点 ~1e-10 / 固有値 ~1e-8）。

デモ実行 (摂動から 30 歩、収束):

```
stride 0:  deviation 5.93e-3
stride 29: deviation 4.87e-9   ← リミットサイクル収束
```

## Phase 4a の結果 (Basin of attraction 可視化)

各 Poincaré 断面点から stride 写像を最大 20 反復し「収束 / 転倒 / 未決」を分類して
basin（吸引域）を可視化する。統制ペア（点足 vs 円弧足）で **「円弧足が basin を
広げるか」** を `basin_fraction`（窓内の収束セル率＝面積の代理スカラー）で比較した。

**較正した窓・解像度**: 4つの 3D 断面モデルは公平比較のため**同一窓 half_widths=(0.16, 0.65)**
（θ_st 軸 ±0.16、θ̇_st 軸 ±0.65、各不動点中心）を共有。res=41 較正で各モデルの収束領域は
概ね窓内に収まり、残る漏れは +θ_st / -θ̇_st 方向の細い舌状部の先端 0〜3 セルのみ
（完全包含は可視構造を潰すため不採用）。simplest は文献ゲートが固定する 2D 窓 (0.08, 0.08)。
本番は **解像度 R=60、14 並列、所要 ≈35 分**（R=80 は ≈17 分/モデルで 10 分/モデル基準超のため R=60 採用）。

| モデル | 足形状 | basin_fraction (R=60) |
|--------|--------|----------------------|
| simplest | 点足 | 0.029 |
| compass | 点足 | 0.061 |
| rocker_compass | 円弧足 | **0.117** |
| kneed | 点足 | **0.108** |
| rocker_kneed | 円弧足 | 0.061 |

3D 窓はペア内で同一面積（(2·0.16)·(2·0.65)=0.416）なので fraction が面積比較そのもの。

**円弧足 vs 点足の判定（正直に）**:
- **Compass ペア**: rocker_compass 0.117 > compass 0.061（相対 +92%）→ **円弧足が basin を約 1.9 倍に拡大**。
  仮説を支持、robust（差は 15% 基準を大きく超える）。
- **Kneed ペア**: rocker_kneed 0.061 < kneed 0.108（相対 −43%）→ **円弧足が basin を約 0.57 倍に縮小**。
  仮説に**反する** robust な結果。較正時に点足 kneed の方が窓端漏れが大きく（真の basin はさらに大）、
  この反転は clipping のアーティファクトではなく実在。**failure ではなく要調査の finding** として記録する
  （膝つき機構と円弧足の相互作用が compass 系と異なる可能性）。

**simplest のフラクタル的 basin**: basin は窓を斜めに横切る非常に薄い連結スリバーで、
fraction 0.029。Schwab & Wisse (2001) の "very small and thin, fractal-like" 記述と定性一致
（薄すぎて物理エンジン上のランダム探索ではリミットサイクルに到達できない、という Heron 教訓の根拠）。

これは Phase 4a の成果物（basin 可視化）。次は **Phase 4b（物理エンジンへの seed 受け渡し）** で
解析的に求めた不動点を初期値として MuJoCo 等に渡す逆順戦略を検証する。

montage PNG は `data/runs/` 配下（**gitignore 対象＝ローカル成果物**）:
`data/runs/20260616_144724_basin_compare/basin_compare.png`（+ `basin_fractions.json`）。

## Phase 4a.1 の結果 (受動歩行→車輪の連続軸、issue #19)

「コンパスは車輪への漸近」という概念枠組み（rimless wheel が最簡受動歩行モデル）を、
rocker_compass の円弧足半径 R を固定勾配 γ=0.030 で R∈(0, L=1) に掃引して検証した。
`hip = C − (L−R)·down(θ)` より R→L で hip が足中心に乗る＝半径 L の車輪の極限。
各 R で相対損失 δ(R)=(KE直前−KE直後)/KE直前（主指標、車輪極限で →0）、安定性 max|λ|(R)、
basin_fraction(R) を測定。R=0.3 起点に continuation で両方向へ追跡。

| R | δ（効率損失） | max\|λ\|（安定性） | basin_fraction |
|---|--------------|-------------------|----------------|
| 0.05 | 0.303 | 0.437 | 0.037 |
| 0.20 | 0.292 | **0.373**（最安定） | 0.072 |
| 0.40 | 0.272 | 0.517 | 0.200 |
| 0.60 | 0.241 | 0.703 | **0.430**（最大） |
| 0.70 | **0.210**（最効率） | 0.867 | 0.026（崩壊） |
| 0.75 | — | — | 歩行消失 |

**3指標が別々の R で最適化される三つ巴のトレードオフ**（正直な判定）:
- **効率（δ）**: R↑ で単調減少（0.303→0.210）。車輪に近いほど効率↑ → **H3 支持**。
- **安定性（max|λ|）**: R≈0.2 で最小（最も強く安定）、車輪に近づくと急速に不安定化（0.867 @ R=0.70）。
- **basin**: R≈0.6 でピーク（0.430）後、R=0.70 で崩壊（0.026）→ **H1（単調増）は不成立**。
- **車輪極限**: R=0.75 でリミットサイクル消失＝**純粋な車輪には到達できない** → **H2 支持**。
- **エネルギー収支チェック**: 機械的 COT(R) = 衝突損失/(総質量·g·一歩水平距離) が全 R で sin γ=0.030 に
  **7桁一致**（リミットサイクル上で重力入力＝衝突損失の厳密な健全性確認。副次的に step_length 幾何の正しさも裏付け）。

結論: 円弧足を大きくするほど効率は車輪へ向かって単調改善するが、安定性と頑健性（basin）は
犠牲になり、純粋な車輪に達する手前（R≈0.75）で歩行が消滅する。「車輪が最も効率的だが、
歩行を維持するには車輪極限から離れる必要がある」という効率・安定性・頑健性のトレードオフ。

成果物（gitignore 対象）: `data/runs/20260616_225010_wheel_limit/`
（`R_sweep.png` 3段曲線・`basin_R_montage.png`・`R_sweep.json`）。

## Phase 5a の結果 (動力付き simplest walker, Kuo 2002)

受動歩行は勾配の重力でしか駆動できず、平地 γ=0 では歩けない。Kuo 2002 の simplest walker に
pre-emptive **push-off**（後脚軸方向の撃力 P）を足し、**受動では不可能だった平地歩行を実現**した
（能動歩行アークの第一段）。連続相は simplest と同一（γ 経由）で、heel-strike だけを
「push-off 撃力 P → 角運動量保存衝突 → 脚交換」の合成写像に差し替える（P=0 で受動に厳密退化）。

平地 γ=0 では θ 部分系が push-off と分離し、素朴な shoot は静止解に落ちる。そこで受動 γ=0.009
サイクルから γ→0・push_off→target へ **continuation** して平地サイクルを発見した。

- **平地リミットサイクル（γ=0, push_off=0.115）**: y*=(0.32659862, −0.33950461)、max|λ|=0.650（安定）。
  不動点近傍からの 25 歩シミュで deviation が単調に減衰（4.0e-3 → 5.9e-7）し、平地サイクルへ収束。
- **push-off→0 退化ゲート**: P=0 で Phase 1 検証済み受動サイクル (0.2003109, −0.1998325) に全 stride 一致。
- **エネルギー収支**: push-off 仕事 P²/2 = 一歩衝突損失（~1e-13 で一致）。
- **Kuo 2002 provenance**: abstract を逐語確認（"toe-off impulse is four times less costly" + power law、
  `src/crane/references_kuo.py`）。本文は paywall のため照合用の単一印刷数値なし＝数値ゲートは
  エネルギー収支（P²/2）に委譲。

成果物（gitignore 対象）: `data/runs/<timestamp>_powered_P0.115/`（`walk.mp4`・`phase_portrait.png`・`meta.json`）。

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

# Phase 3.5 デモ歩行: 円弧足 Rocker-foot Compass（walk.mp4 + phase_portrait.png + meta.json）
uv run python scripts/walk_rocker.py [--gamma-deg <published>] [--strides 30] [--perturb 0.005]

# Phase 3.6 デモ歩行: 円弧足 Rocker-foot Kneed Walker（4セグメント walk.mp4 + phase_portrait.png + meta.json）
uv run python scripts/walk_rocker_kneed.py [--gamma-deg <published>] [--strides 30] [--perturb 0.005]

# Phase 4a basin 比較: 5モデルの吸引域 montage + basin_fraction 表（data/runs/ に出力）
uv run python scripts/basin_compare.py [--resolution 60] [--workers N]

# Phase 4a.1 車輪への漸近: rocker_compass の R 掃引で δ(R)・安定性・basin 曲線（data/runs/ に出力）
uv run python scripts/wheel_limit.py [--dr 0.05] [--basin-res 50] [--workers N]

# Phase 5a デモ歩行: 動力付き simplest walker の平地 γ=0 歩行（continuation で平地サイクル発見）
uv run python scripts/walk_powered.py [--push-off 0.115] [--strides 30] [--perturb 0.01]

# テスト
uv run pytest
```

## ドキュメント

- 設計: `docs/2026-06-12-crane-design.md`
- フェーズゲート: `GOALS.md`

## Heron との関係

Heron (Genesis/MuJoCo ポート) で確認されたこと: 物理エンジン上でのランダム探索は
basin が薄すぎてリミットサイクル候補に到達できない。
Crane はその教訓から出発し、解析的にサイクルを求めてからシミュレーションに渡す逆順戦略を採る。
