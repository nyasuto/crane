# Phase 5a 設計ドキュメント: 動力付き simplest walker（Kuo 2002）

**作成日:** 2026-06-17
**ステータス:** 承認済み（ぽんぽこ殿が設計一任・実装移行を承認 2026-06-16）
**対応:** 能動歩行アークの第一段（受動歩行アーク Phase 0〜4a.1 完了後）
**続き:** Phase 5b（能動 rocker_compass で「受動で良い個体は能動でも強いか」を検証）

---

## 背景と動機

受動歩行アークで「効率・局所安定性・basin の三つ巴トレードオフ」「basin こそ有限外乱への頑健性」を
得た。次の問い: **能動制御を足すと、basin/効率の良い受動個体は能動でも強いのか**（ぽんぽこ殿の仮説）。
これを検証する前に、Crane の流儀（まず公表モデルに錨を打つ）に従い、能動の最小要素を
**公表値で検証できる動力付き simplest walker（Kuo 2002）**で確立する。

受動歩行は斜面の重力でしか歩けない。**push-off（撃力踏み出し）を足すと平地（γ=0）を歩ける**——
これが受動には不可能だった質的に新しい能力であり、本フェーズの headline。

## ゴール

simplest walker に pre-emptive 撃力 push-off を足し、**平地 γ=0 を push-off 駆動で歩く**リミットサイクルを
発見・安定性証明し、Kuo 2002 の公表関係＋ push-off→0 退化ゲートで検証する。

## スコープ

### やること
- simplest walker をミラーした `powered_simplest.py`（simplest.py は不変）に push-off 撃力写像を実装
- 平地リミットサイクル発見、push-off P 掃引で速度・効率・安定性を測定
- Kuo 2002 公表関係＋ push-off→0 退化ゲート＋エネルギー収支で検証
- 平地歩行アニメ

### やらないこと（YAGNI）
- basin 測定 → Phase 5b（受動 vs 能動の順位比較）に集約
- 連続トルク制御・hip トルク（push-off 撃力のみ）
- rocker/kneed の能動化 → 5b 以降
- 斜面 × push-off（平地に集中）

## 新規物理: pre-emptive 撃力 push-off

Kuo 2002 の核心は**衝突直前の撃力 push-off**（後脚軸方向の撃力で COM 速度を衝突前に再配向し、
後続の衝突損失を最小化。衝突後押し出しや hip トルクより同速度を 1/4 コストで出せる）。

heel-strike 条件（φ−2θ=0, θ<0）の瞬間に合成写像:
1. **後脚（trailing 脚）軸方向の撃力 P で速度ジャンプ**（KE 注入）
2. 前脚の角運動量保存衝突（既存の散逸写像）
3. 脚交換

### 導出（reduced 座標 y=(θ, θ̇)、M=1, l=g=1）
受動の hip 速度は stance 脚に垂直（大きさ |θ̇|）。push-off は trailing 脚軸方向の撃力 P を hip 点質量に
与える（push-off 前の脚軸方向速度は 0）。衝突は新 stance 脚（leading 脚、interleg 角 2θ）軸方向成分を除去:

- 受動の垂直成分 → 新脚垂直への射影 = cos(2θ)（＝既存の c·θ̇）
- push-off（trailing 脚軸方向）→ 新脚垂直への射影 = sin(2θ)

よって `c = cos(2θ)`, `s = sin(2θ)` として

```
θ̇⁺ = c·θ̇ + s·P                  # P=0 で c·θ̇（受動）に厳密一致
powered_heelstrike_map(x, P):
    return [-θ, -2θ, θ̇⁺, (1-c)·θ̇⁺]
```

P=0 で `[-θ, -2θ, c·θ̇, (1-c)c·θ̇]` ＝ 既存 `heelstrike_map` に一致（退化ゲートが構造的に成立）。
符号: θ̇<0（前進）, θ<0 ⇒ s=sin(2θ)<0, よって P>0 で θ̇⁺ がより負＝加速（前進にエネルギー注入）。✓

### エネルギー（検証アンカー）
push-off 前の脚軸方向速度 0、後 P（M=1）⇒ **push-off 仕事 = P²/2**。これは Kuo の「push-off 仕事 ∝
push-off²」と一致し、強い整合チェックになる。平地リミットサイクル上では **push-off 仕事 P²/2 = 一歩衝突損失**。

## 平地リミットサイクル
γ=0、連続相は既存 dynamics（γ=0 で theta_dd=sin θ）。push-off P>0 で `find_limit_cycle`（既存 Newton
shooting）で y*=(θ*, θ̇*) を求める。重力でなく push-off がエネルギー供給。断面・lift/project は受動と同形。

## 検証ゲート

1. **push-off→0 退化ゲート（最強）**: powered model（γ を param に持つ）が P=0 で全 stride・不動点・固有値が
   Phase 1 検証済み受動 simplest walker と一致（γ=0.009 で (0.2003109, −0.1998325)）。検証済み ground truth に錨。
2. **Kuo 2002 効率/速度ゲート**: push-off ↔ 速度 ↔ 力学コストの公表関係を照合。**Kuo 2002 本文を実装時に
   取り寄せ provenance 記録、記憶からの数値は不使用**。取得不能なら定性＋エネルギー収支ゲートに委譲。
3. **エネルギー収支**: 平地サイクルで push-off 仕事 P²/2 = 一歩衝突損失（独立計算で照合）。
4. **安定性**: 動力サイクルで max|λ|<1。
5. **ぽんぽこ殿の目視判定**: 平地 walk.mp4 が「push-off で前進する歩行」に見える。

## アーキテクチャ / ファイル

- `src/crane/models/powered_simplest.py`（新規）: `PoweredSimplestParams(gamma, push_off)`、
  `powered_heelstrike_map`、`make_powered_simplest`。simplest.py をミラーし**改変しない**
  （rocker_kneed が kneed をミラーした流儀）。連続 dynamics は simplest と同じ（γ 経由）。
  P=0 で simplest と数値一致。
- `src/crane/references_kuo.py`（新規）: Kuo 2002 provenance ＋公表数値（push-off↔速度↔コスト）。
- `scripts/walk_powered.py`（新規）: 平地動力サイクル発見→多歩シミュ→平地アニメ＋meta.json。
- テスト: `test_references_kuo`、`test_powered_simplest_impact`（push-off 写像・エネルギー P²/2・P=0 退化）、
  `test_powered_simplest_degenerate`（P→0 で Phase 1 受動に全 stride 一致）、`test_powered_simplest_cycle`
  （平地サイクル発見・安定性・push-off 仕事=衝突損失・Kuo 照合）、`test_viz_powered`。

## データフロー
1. `PoweredSimplestParams(gamma=0, push_off=P)` → `make_powered_simplest`。
2. `find_limit_cycle` で平地不動点。
3. `stride` の x_strike/x_end と push-off 仕事 P²/2 でエネルギー収支照合。
4. P 掃引で速度・効率・max|λ|。
5. walk_powered.py で平地アニメ＋meta。

## リスクと対策

| リスク | 対策 |
|---|---|
| push-off 写像の符号/係数ミス | P=0 退化ゲート（受動と全 stride 一致）＋ push-off 仕事=P²/2 ＋ 平地で仕事=損失、の三重ガード |
| Kuo 2002 本文が paywall で数値取得不可 | 定性関係＋エネルギー収支＋退化ゲートに委譲。取得可否を honest に記録（捏造禁止） |
| 平地サイクルが Newton で収束しない | 受動 γ=0.009 サイクルを連続変形（γ↓0 と P↑ を同時に）して seed を供給 |
| simplest.py を壊す | powered を別ファイルにミラー、simplest.py は byte 不変（テストで担保） |

## フェーズ命名
**Phase 5a: 動力付き simplest walker（Kuo 2002）**。続く 5b で能動 rocker_compass の仮説検証。
