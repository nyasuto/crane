# Phase 4a.1 設計ドキュメント: 受動歩行→車輪の連続軸（rocker-foot R 掃引）

**作成日:** 2026-06-16
**ステータス:** 承認済み（ぽんぽこ殿レビュー 2026-06-16）
**対応 issue:** [#19](https://github.com/nyasuto/crane/issues/19)
**前段:** Phase 4a（basin 可視化、PR #18 マージ済み）
**関連:** Phase 3.5（rocker_compass、R=0.3 の検証済みサイクル）

---

## 背景と動機

Phase 4a で、統制比較により **円弧足は basin を広げるが普遍則ではない**ことが判明した
（compass ペアで ×1.9、kneed ペアで ×0.57）。一次資料照合（McGeer 1990b / Hsu Chen 2007）と
合わせ、basin 拡大は膝でなく**円弧足**由来という論旨に整合した。

ここからぽんぽこ殿との考察で次の概念枠組みが生まれた：

```
点足コンパス ──(R↑)── 円弧足コンパス ──(R→脚長 L)── 完全な車輪
   basin 小                basin 大                 歩行が消える（転がるだけ）
```

rimless wheel（スポークだけの車輪）が最簡の受動歩行モデルであり、コンパス2脚は
「2本スポークの rimless wheel」、円弧足は足裏を車輪リムにして転がり接触させる＝最も車輪的。
heel-strike の衝突損失を円弧足の転がり接触が緩和する＝車輪（損失ゼロ極限）に近づく＝
摂動許容が広がる。Phase 4a で rocker_compass が最大 basin だったのは「最も車輪的なモデルが
最も頑健」の実測的帰結だった。

rocker_compass の幾何 `hip = C − (L−R)·down(θ)`（C=[−R·θ, R]）より、**R→L で hip が
足中心高さ L に一致＝半径 L の車輪の極限**そのもの。よって R∈(0, L) の掃引が
幾何学的にクリーンな「車輪への漸近」軸になる。

## ゴール

固定勾配 γ=0.030（Phase 3.5 と同一）で rocker_compass の R を R∈(0, L=1) で掃引し、
各 R のリミットサイクル・安定性・効率・basin を測定して **「コンパスは車輪への漸近」**
という概念枠組みを曲線で検証する。

## スコープ

### やること
- R 掃引でリミットサイクルを continuation 追跡（前の解を seed）
- 各 R で効率指標 δ(R)、安定性 max|λ|(R)、basin_fraction(R) を測定
- 歩行が消失/不安定化する R（車輪極限側の境界）を特定
- 曲線（δ↓・basin↑）として可視化し「車輪への漸近」を検証

### やらないこと（YAGNI）
- rocker_kneed の R 掃引（膝の R 依存は issue #6 後に別途）
- γ_min(R)（各 R での最小勾配）探索 — 魅力的だが高コスト、次フェーズ候補
- 物理エンジン連携（Phase 4b）

## 効率指標の物理（重要）

リミットサイクル上では1歩のエネルギー収支が閉じる：**重力が1歩で与えるエネルギー ＝
heel-strike 衝突で失うエネルギー**（揺動相は保存系）。したがって
機械的 COT = 衝突損失/(m·g·水平距離) = sin γ となり、**R によらず勾配で固定**される。
ゆえに機械的 COT は「車輪への漸近」を判別できない（が、エネルギー収支の内部チェックには使える）。

判別力のある指標は**一歩あたり相対損失** δ(R)：

```
δ(R) = (KE_直前 − KE_直後) / KE_直前
```

無次元で、車輪極限で δ→0。同じ絶対損失でも円弧足ほど速く歩き KE が大きい＝δ が小さい。
R とともに変化し、車輪への漸近を直接捉える。**主指標は δ(R)**、機械的 COT は副次の
内部整合チェック（≈ sin γ を確認）とする。

## アーキテクチャ

### 新規モジュール `src/crane/efficiency.py`
model 非依存の効率指標。`kinetic_energy` を callable として受け取る。

```python
def relative_loss(ke_pre: float, ke_post: float) -> float:
    """δ = (KE_pre - KE_post)/KE_pre。車輪極限で →0。"""

def step_collision_loss(x_strike, x_end, kinetic_energy) -> tuple[float, float, float]:
    """(loss, ke_pre, ke_post) を返す。loss = KE(x_strike) - KE(x_end)。"""

def mechanical_cot(loss: float, m: float, g: float, step_length: float) -> float:
    """衝突損失 / (m·g·一歩水平距離)。リミットサイクル上では ≈ sin γ になるはず。"""
```

一歩水平距離は rocker_compass の hip_x = −R·θ_st − (L−R)·sin θ_st の幾何から算出
（接触フレームは heel-strike でリセット）。幾何が煩雑なら mechanical_cot は secondary
扱いとし、δ(R) を主軸に据える。

### 新規スクリプト `scripts/wheel_limit.py`
二段掃引を実行:
1. **効率・安定性・不動点の細かい掃引**: R=0.02 から ΔR≈0.05 で L 近くまで。各 R で
   前の解を seed に `find_limit_cycle`。収束失敗 or 不安定化（max|λ|≥1）で歩行消失点を特定。
   記録: y\*(R), max|λ|(R), δ(R), 機械的 COT(R)。
2. **basin の粗い掃引**: R∈{0.05, 0.2, 0.4, 0.6, 0.8}＋消失直前、各 R で `basin_slice`
   （Phase 4a 再利用）を**全 R 共通の固定窓**・res≈50 で。fraction と絶対面積を記録。

出力（`data/runs/<timestamp>/`、gitignore 対象）:
- `R_sweep.png`: δ(R)・max|λ|(R)・basin_fraction(R) の曲線（多段プロット）
- `R_sweep.json`: 全測定データ
- `basin_R_montage.png`: 各 R の basin 並置

### 既存資産の再利用
- `find_limit_cycle`（search.py）: continuation の各 R で seed 付き呼び出し
- `stride`（x_strike, x_end）/ `rocker_compass.kinetic_energy` / `energy`
- `basin_slice`（basin.py、Phase 4a）
- references_mcgeer（m=1, m_h=0, c=0.37, ρ=0.32, L=1.0, γ=0.030）

## データフロー
1. references_mcgeer から rocker_compass パラメータ（γ=0.030 固定）。
2. R=0.3 から開始（Phase 3.5 検証済み seed）、R を上下に continuation で広げる。
3. 各 R: find_limit_cycle → stride で δ・COT 算出 → 安定性。
4. 粗い R 集合で basin_slice。
5. 曲線・montage・JSON を出力。

## 検証ゲート

単一の公表曲線は存在しないので、以下で締める:

1. **R=0.3 アンカー（最強）**: 掃引が R=0.3 で Phase 3.5 検証済みサイクル
   (0.30844, −1.26256, −0.87914) を再現（既存 verified ground truth に R 軸を固定）。
2. **R→0 退化**: 小 R で点足 compass 近傍に連続接続（Phase 3.5 の R→0 退化と整合）。
3. **エネルギー収支**: 機械的 COT(R) ≈ sin γ（全 R）。δ(R) と独立な健全性チェック。
4. **rimless wheel provenance**: 「rimless wheel = 最簡受動歩行モデル」を McGeer 1990a で
   裏取りし references に provenance 記録（**記憶からの断定はしない**。取得できなければ
   その旨を明記し主張を控える）。
5. **ぽんぽこ殿の目視判定**: δ(R)↓・basin(R)↑ の曲線が「車輪への漸近」を支持して見えるか。

H1（basin↑単調）/ H2（消失点存在）/ H3（δ↓単調）は仮説。結果は仮説に合わせ込まず正直に報告。

## リスクと対策

| リスク | 対策 |
|---|---|
| R↑ で continuation が分岐/発散 | 前の解を seed に小 ΔR。失敗 R を消失点として記録（H2 の観測そのもの） |
| basin 窓が R で basin を切り落とす | 全 R を含む広めの固定窓。fraction と絶対面積の両方を記録 |
| 機械的 COT の幾何（step length）が煩雑 | δ(R) を主軸に。COT は secondary、算出困難なら energy() 収支で代替チェック |
| R→L 近傍で数値が病的（脚長 0 近傍） | L−R が小さい領域は慎重に。消失点付近は ΔR を細かく |
| 「車輪が最簡受動歩行」の provenance 不在 | McGeer 1990a で裏取り、取れなければ主張を controlled に留める |

## フェーズ命名
GOALS に **Phase 4a.1: 車輪への漸近（R 連続軸）** を追加（Phase 4b 物理エンジンとは独立の小フェーズ）。
