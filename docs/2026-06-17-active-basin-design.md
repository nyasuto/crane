# Phase 5c 設計ドキュメント: 能動 vs 受動 basin 比較（本命仮説）

**作成日:** 2026-06-17
**ステータス:** 承認済み（ぽんぽこ殿が設計一任・実装移行を承認 2026-06-17）
**前段:** Phase 4a.1（受動 basin(R)）、Phase 5b（動力付き rocker_compass）
**位置づけ:** 能動歩行アークの締めくくり。ぽんぽこ殿の本命仮説の検証。

---

## 背景と動機

ぽんぽこ殿の仮説:
> 受動で basin/効率の良い個体は、能動でも強いのか？ それとも制御は形態差を消すのか？

Phase 4a.1 で受動 rocker_compass の basin(R) を得た（γ=0.030、basin 最大は R≈0.6）。
Phase 5b で能動化の道具（push-off 増強モデル）を確立した。本フェーズは両者を合流させ、
**同一勾配 γ=0.030 で各 R に push-off を載せ、能動 basin(R) を受動 basin(R) と比較**して仮説を検証する。

## ゴール

γ=0.030 で R を掃引し、push_off P∈{0, 0.04, 0.08} の能動 basin_fraction(R, P) を測定。
受動（P=0）と比較し、(1) 制御は受動の順位を保存するか、(2) 制御は R 間の形態差を平坦化するか、を答える。

## スコープ

### やること
- (R, P) 格子で能動リミットサイクルを求め `basin_slice`（4a.1 と同一窓）で basin_fraction を測定
- 受動（P=0）vs 能動（P>0）の basin(R) 曲線・順位・ばらつきを比較
- push_off→0 で受動 basin を再現する内部ゲート

### やらないこと（YAGNI）
- 新規物理（basin_slice ＝ Phase 4a、powered_rocker_compass ＝ Phase 5b を再利用）
- push-off の R ごとチューニング（固定 P で骨格を出す。何を揃えるかは別研究問題）
- rocker_kneed / simplest の能動 basin

## 方法／データフロー

対象: R ∈ {0.05, 0.20, 0.40, 0.60}（4a.1 の basin 測定点。R=0.7 は受動 basin が崩壊済みで除外）、
push_off P ∈ {0.0, 0.04, 0.08}。窓: 4a.1 と同一 half_widths=(0.16, 0.65)、各サイクルの不動点中心、
解像度 res≈50（重ければ 40）。固定 γ=0.030。

各 (R, P):
1. `make_powered_rocker_compass(PoweredRockerCompassParams(m=1,m_h=0,c=0.37,rho=0.32,L=1,R,gamma=0.030,push_off=P))`。
2. `find_limit_cycle` を受動不動点近傍 seed で収束（push-off は受動サイクルの摂動なので直接収束）。
   収束しなければ push_off を 0 から continuation。
3. `basin_slice(make_powered_rocker_compass, params, fp.y, axes=(0,1), half_widths=(0.16,0.65),
   resolution=res, model_name=f"R{R}_P{P}", n_workers=cores)` で basin_fraction を測定。
4. basin_fraction(R, P) を表・曲線・montage に集約。

## 仮説の指標
- **basin_fraction(R, P) 曲線**: P ごとに 1 本。受動（P=0）が基準。
- **順位保存?**: argmax_R basin が P>0 でも R≈0.6 のままか。
- **平坦化?**: R 間の basin のばらつき（変動係数 = std/mean across R）が P 増加で縮むか
  ＝制御が形態差を消すか。
- 結果は仮説に合わせ込まず正直に報告（保存/非保存、平坦化/非平坦化のいずれでも finding）。

## 検証ゲート
1. **push_off→0 内部ゲート（最強）**: P=0 で能動 basin(R) が Phase 4a.1 の受動 basin(R) を再現
   （powered モデルの push_off=0 ＝ 受動 rocker_compass なので basin も一致）。検証済み結果に錨。
   テストは 1 R 点・低解像度で「powered(P=0) basin ＝ rocker_compass basin」を確認。
2. **仮説の正直な報告**: 順位・平坦化を実測通り記録。
3. **ぽんぽこ殿の目視判定**: basin(R,P) 比較図（曲線＋montage）が読み取れ主張を支持して見えるか。

## アーキテクチャ／ファイル
- `scripts/active_basin_sweep.py`（新規）: (R,P) 掃引、各 `basin_slice`、basin_fraction(R,P) の
  曲線＋ basin montage ＋ JSON を `data/runs/` に出力。順位・変動係数も算出・記録。
- テスト `tests/test_active_basin_gate.py`: 低解像度で powered(P=0) の basin が
  rocker_compass（受動）の basin に一致することを確認（内部ゲート）。
- README/GOALS の Phase 5c セクション。

## リスクと対策
| リスク | 対策 |
|---|---|
| 計算量（12 basin） | res40 へ落とす／並列。曲線が主眼なので中解像度で十分 |
| 一部 (R,P) でサイクル非収束 | push_off=0 から continuation で seed 供給。非収束は記録（finding） |
| 窓が basin を切り落とす | 4a.1 と同一窓で比較整合。収束領域が窓端に接したら記録 |
| 内部ゲートが重い | 1 R 点・低解像度（res~20）に限定 |

## フェーズ命名
**Phase 5c: 能動 vs 受動 basin 比較（本命仮説）**。能動歩行アークの締めくくり。
