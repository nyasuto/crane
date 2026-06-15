# Crane 設計ドキュメント: Rocker-foot Kneed Walker（McGeer 1990b 原機械, math-first）

**作成日:** 2026-06-16
**ステータス:** ぽんぽこ殿レビュー待ち
**親設計:** `docs/2026-06-12-crane-design.md`
**前提フェーズ:** Phase 3（点足 kneed, Hsu Chen 2007）+ Phase 3.5（円弧足 compass, McGeer 1990a）
**位置づけ:** Phase 3.6（Phase 3 と Phase 3.5 の合流。物理エンジン Phase 4 はスコープ外のまま）

---

## 背景と動機

親設計の核心的観察: 「McGeer のオリジナル機械が成立した条件（rocker foot + 質点集中 +
knee latch + 質量分布）は**全部揃って初めて歩く**。部分的な再現では歩かない」。

Crane はこれまで要素を一つずつ隔離して検証してきた:
- Phase 3: 膝付き4相 hybrid（点足）— Hsu Chen 2007 と照合、完了
- Phase 3.5: 円弧足の転がり接触（膝なし compass）— McGeer 1990a と照合、完了

本フェーズは両者を**合流**させ、McGeer 1990b "Passive Walking with Knees" の
**原機械（膝 + 円弧足）を完全再現**する。受動歩行の王道機械が数学的に成立することを示す。

## ゴールと完了ゲート

**円弧足・膝付き2脚のリミットサイクルを数学的に発見し、安定性を固有値で証明する。**

完了ゲート（Phase 1-3.5 と同じ流儀）:

1. Newton shooting で不動点を発見（収束履歴をログ）
2. 固有値で漸近安定性（max|λ| < 1）を確認
3. **R→0 退化ゲート（主軸・最強の検証）**: Phase 3 検証済み点足 kneed (Hsu Chen) の
   不動点 (0.23859, −1.10959, −0.05715) と固有値に一致
4. McGeer 1990b の公表数値との照合（実装の文献タスクで精読。6 ページの会議論文なので
   取れる数値だけを provenance 付きで照合。図読取値は ±X% ラベル）
5. ぽんぽこ殿の歩容判定: walk.mp4 が「円弧足で転がりつつ膝を使って歩く」ように見える

## 検証アンカー（文献）

主軸は **R→0 退化ゲート**（自前の Phase 3 検証済み実装への一致）。
文献は補助:
- **McGeer, T. (1990) "Passive Walking with Knees", Proc. IEEE ICRA**（Cornell/Ruina ラボ
  公開 PDF、6 ページ、2026-06-16 取得済み）。原機械の rocker-foot kneed。
- 退化先: **Hsu Chen 2007 点足 kneed**（Phase 3、`references_kneed.py`、検証済み）。

McGeer 1990b は短い会議論文のため印刷数値が薄い可能性がある。Phase 3.5 で McGeer 1990a の
固有値が test machine 由来で直脚モデルに使えなかったのと同様、**取れる数値だけをゲートにし、
無ければ「存在 + 安定性 + R→0 退化 + 歩容判定」を主ゲートに縮退**する（provenance 規律）。
数値はプランに直接書かない。

## 物理モデル

Phase 3 kneed に「stance 足の転がり」だけを足す。状態・相機械・質量モデル・knee-strike は不変。

- **状態**: Phase 3 kneed と同一 6D `x = [θ_st, θ_th, θ_sh, θ̇_st, θ̇_th, θ̇_sh]`（slope 法線
  基準の絶対角）。断面 `y = (θ_st, θ̇_st, θ̇_sw)`。lift/project は kneed と同一。
- **4相機械**: unlocked 3-DOF swing → knee-strike（θ_sh−θ_th=0）→ locked 2-DOF →
  heel-strike（θ_st+θ_sw=0）。Phase 3 と同一。
- **質量モデル**: Phase 3 kneed の Hsu Chen 点質量（hip m_h、各脚に thigh m_t・shank m_s の
  点質量）+ **足半径 R**。R→0 で点足 kneed に厳密縮退。
- **stance 足の転がり（唯一の追加）**: stance 円弧足の曲率中心 C は高さ R、転がりで
  `C_st = [−R·θ_st, R]`。`hip = C_st − (L−R)·down(θ_st)`（rocker_compass と同一、R→0 で
  点足 kneed の hip 式に一致）。stance 脚の全質点は hip 経由でこの転がりを継承。
  **unlocked / locked の両相の力学に適用する**（stance 脚は両相とも単一剛体リンク長 L）。
- **knee-strike は不変**: swing 脚内部の衝突（脛が大腿と整列）。足形状と無関係なので
  Phase 3 の定式化（系全体を stance 接触回り + swing 脚を hip 回りで角運動量保存）をそのまま。
  ただし「stance 接触」は転がり接触点 `P_st=[−R·θ_st, 0]` に置き換わる。
- **heel-strike**: swing 円弧足の転がり接触 `P_sw=[C_sw_x, 0]` 回りで角運動量保存 + 脚交換
  （rocker_compass と同型 + kneed の rigid-through-collision θ̇_sh⁺=θ̇_th⁺）。両足同半径 R の
  幾何対称性から strike 面は `θ_st+θ_sw=0`（kneed/compass と同一、厳密）。
- **R→0** で Phase 3 点足 kneed に厳密縮退（退化ゲートの土台）。

## コード構成（既存資産を最大再利用）

| パス | 責務 |
|---|---|
| `src/crane/models/rocker_kneed.py` | 新規。kneed.py の構造をミラーしつつ hip を転がりに。自前 `_build()`（derive 層は無改変）、`RockerKneedParams`（kneed 質量パラメータ + R）、`make_rocker_kneed()` |
| `src/crane/references_mcgeer_knees.py` | 新規。McGeer 1990b 公表値（provenance、図読取ラベル） |
| `src/crane/viz.py` | `animate_rocker_kneed` 追加（4セグメント + 円弧足の転がり = Phase 3 viz + Phase 3.5 viz の合成） |
| `scripts/walk_rocker_kneed.py` | 新規。Phase 3.6 デモ |
| 再利用 | `stride.py`・`search.py`・`derive/lagrange.py`・`derive/impact.py`・`references_kneed.py`（R→0 退化ゲートの照合先） |

derive レイヤー（lagrange/impact）・kneed.py・rocker_compass.py は**無改変**。

## 検証戦略（三段防御、Phase 3/3.5 と同じ）

1. **不変量テスト**: 平衡・unlocked/locked 両相のエネルギー保存（drift < 1e-7）・
   knee-strike / heel-strike の KE 散逸＆脚交換幾何・零速度不動点・R→0 で両相力学＆両衝突が
   Phase 3 点足 kneed に一致
2. **退化ゲート（主軸・最強）**: R→1e-9 で全 stride 不動点・固有値が Phase 3 検証済み点足 kneed
   (0.23859, −1.10959, −0.05715) に一致
3. **文献ゲート**: McGeer 1990b パラメータでリミットサイクル存在、取れる公表数値（不動点・
   step period・固有値）を照合（無ければ存在＋安定性＋R→0 退化に縮退）
4. **歩容判定**: walk.mp4（円弧足の転がり + 膝屈曲）をぽんぽこ殿が判定

## リスクと対策

| リスク | 対策 |
|---|---|
| 転がりが unlocked / locked の両相に正しく入っているか | R→0 退化ゲートを両相力学で個別検証 + 全 stride で検証。符号/相のミスはここで捕捉 |
| knee-strike の stance 接触が転がり点に変わる影響 | R→0 退化で Phase 3 knee-strike に一致を確認。m_s/m_t=0.1 で影響小だが検証で固定 |
| McGeer 1990b（6 ページ）の数値が薄い | R→0 退化ゲートを backbone に。McGeer 値は取れた分だけ（Phase 3.5 と同じ縮退規律） |
| basin が点足 kneed より細い | Armijo（導入済み）+ R=0（Phase 3 不動点）からの R-continuation で seed 供給 |
| McGeer 1990b と Hsu Chen の質量/足パラメータ規約差 | 文献タスクで独立検算（Phase 3/3.5 と同じ）。退化ゲートは Hsu Chen 規約で行う |

## スコープ外（やらないこと）

- 物理エンジン（MuJoCo / Genesis）再現 — 親設計の Phase 4 のまま（記録のみ）
- rocker foot の3D化・能動制御
- period-doubling / カオス領域の追跡（別 issue）
