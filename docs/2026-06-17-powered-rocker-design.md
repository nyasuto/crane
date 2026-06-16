# Phase 5b 設計ドキュメント: 動力付き rocker_compass（push-off 増強）

**作成日:** 2026-06-17
**ステータス:** 承認済み（ぽんぽこ殿が設計一任・実装移行を承認 2026-06-16）+ 実装前 de-risk で regime を訂正
**前段:** Phase 5a（動力付き simplest walker）、Phase 4a.1（rocker_compass R 掃引・受動 basin）
**続き:** Phase 5c（R 掃引で能動 vs 受動 basin を比較＝本命仮説）

---

## 背景と動機

ぽんぽこ殿の仮説「受動で basin/効率の良い個体は能動でも強いか／制御は形態差を消すか」を検証する
ための土台。Phase 5a で能動の最小要素（pre-emptive push-off）を simplest walker で確立した。
本フェーズは rocker_compass（2-DOF・転がり足）に push-off を足し、5c の比較の道具を揃える。

## 実装前 de-risk の知見（regime 訂正）

scratch で push-off 定式化を検証した結果:
- **push-off 定式化は正しい**（push-off→0 で受動 rocker に厳密一致、diff 0.0）。
- **R=0.3 では平地 γ=0 サイクルは存在しない**: push-off を上げても θ\* が縮み γ→0 手前で fold（消失）。
  「平地歩行」（5a 由来）は rocker では成立しない。これは 4a.1 の車輪極限 fold と響き合う知見。
- **固定勾配 γ=0.030 + push-off 増強は機能**: po∈[0,0.10] で安定サイクルが直接収束。push-off は
  walkable slope を γ=0.030→~0.002 へ拡張（エネルギー注入の証拠）。

→ regime を「**受動と同一の勾配 γ=0.030 で push-off 増強**」に変更。5c の比較（同一 γ で受動 vs 能動 basin）に
最適で、制御の効果を形態から分離できる。

## ゴール

rocker_compass に pre-emptive 撃力 push-off を足した動力付きモデルを確立し、固定勾配 γ=0.030 で
push-off 増強サイクルを発見・安定性証明する。push-off→0 退化ゲートで Phase 3.5 検証済み実装に錨を打つ。

## スコープ

### やること
- rocker_compass をミラーした `powered_rocker_compass.py`（rocker_compass.py 不変）に push-off を実装
- 固定 γ=0.030 で push-off 増強サイクルを発見、安定性・push-off 依存を測定
- push-off→0 退化ゲート＋push-off がエネルギー注入することの確認

### やらないこと（YAGNI）
- R 掃引の能動 basin 比較 → Phase 5c（本命仮説）
- 平地 γ=0 歩行（rocker は fold で到達不可、上記 de-risk）
- rocker_kneed の能動化

## 新規物理: push-off の記号導出（de-risk で検証済み）

pre-emptive push-off ＝ 後脚（old stance = trailing）軸方向の撃力を heelstrike 衝突の直前に注入。

**撃力**: `Î_po = push_off · (P_st − hip) / |hip − P_st|`（後脚の脚線方向、エネルギー注入符号）。
**鍵**: 撃力は接触点 P_st と hip を結ぶ脚線上にあり **hip を通る** → hip まわりのモーメントはゼロ。
よって既存の角運動量保存2式のうち:
- **eq1（系の新接触点まわり角運動量）にのみ** push-off モーメント `M_po = cross2d(P_st − foot_sw_contact, Î_po)` を加える。
- **eq2（swing 脚の hip まわり角運動量）は不変**。

```
M_po = cross2d(P_st - foot_sw_contact, Î_po)
eq1' = L_sys_stance.subs(swap) - L_sys_swfoot - M_po = 0
eq2' = L_sw_hip.subs(swap)   - L_st_hip            = 0   # 不変
```

push_off=0 で M_po=0 → 既存 rocker_compass の衝突写像に厳密一致（de-risk で diff 0.0 確認）。
push_off>0 で post-collision エネルギーが増加（de-risk で確認、符号確定）。

## 固定勾配サイクル
γ=0.030 固定、push_off>0 で `find_limit_cycle`（既存）を受動不動点 (0.30844, −1.26256, −0.87914) を
seed に直接収束（de-risk で po∈[0,0.10] 確認、continuation 不要）。

## 検証ゲート
1. **push-off→0 退化ゲート（最強）**: push_off=0 で全 stride・不動点・固有値が Phase 3.5 検証済み
   受動 rocker サイクル (0.30844, −1.26256, −0.87914) に一致。
2. **動力サイクル存在＋安定性**: γ=0.030・push_off∈[0,0.10] で安定サイクル（max|λ|<1）が収束。
3. **エネルギー注入の確認**: push_off を上げると walkable slope が γ=0.030 より下（平地方向）へ拡張する
   （continuation で γ↓ しつつ po↑ がサイクルを維持＝push-off がエネルギーを供給している証拠）。
4. **ぽんぽこ殿の目視判定**: push-off 増強 walk.mp4 が「円弧足で転がりつつ push-off で歩く」に見える。

## アーキテクチャ／ファイル
- `src/crane/models/powered_rocker_compass.py`（新規）: rocker_compass の `_build` をミラーし、
  push-off モーメント `M_po` を eq1 に加える（rocker_compass.py 不変）。`PoweredRockerCompassParams(..., push_off)`、
  `make_powered_rocker_compass`。push_off=0 で rocker_compass に一致。
- `scripts/walk_powered_rocker.py`（新規）: γ=0.030・push_off で動力サイクル → 多歩 → 円弧足アニメ（既存 viz 再利用）。
- テスト: push-off→0 退化（受動一致）・動力サイクル＋安定性・エネルギー注入（walkable slope 拡張）・viz。

## リスクと対策
| リスク | 対策 |
|---|---|
| push-off 記号導出のミス | de-risk で公式検証済み（push-off→0 で diff 0.0、符号確定）。退化ゲートで常時ガード |
| rocker_compass.py を壊す | powered を別ファイルにミラー、rocker_compass.py は byte 不変（テストで担保） |
| 平地に届かない | regime を γ=0.030 固定の push-off 増強に変更（de-risk 反映）。fold は知見として記録 |

## フェーズ命名
**Phase 5b: 動力付き rocker_compass（push-off 増強）**。続く **Phase 5c** で R 掃引の能動 vs 受動 basin 比較（本命仮説）。
