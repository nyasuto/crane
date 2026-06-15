# Crane 設計ドキュメント: Rocker-foot Compass（McGeer 1990a, math-first）

**作成日:** 2026-06-16
**ステータス:** ぽんぽこ殿レビュー待ち
**親設計:** `docs/2026-06-12-crane-design.md`
**位置づけ:** Phase 3.5（完了済み Phase 3 と、スコープ外の物理エンジン Phase 4 の間の新マイルストーン）

---

## 背景と動機

Crane は受動歩行のリミットサイクルを hybrid dynamics + Poincaré shooting で数学的に発見する
プロジェクト。Phase 1（simplest walker）・Phase 2（点足 compass）・Phase 3（点足 kneed walker）
で、点足/ピン足モデルの3段階を完了した。

McGeer のオリジナル機械が成立した条件のひとつが **rocker foot（円弧足）** である
（親設計 背景節）。点足では支点が固定だが、円弧足は接触点が転がって移動し、連続相の力学と
heel-strike 幾何の双方が変わる。本フェーズは、その rocker foot 機構を**最も単純な 2リンク
compass で隔離して**導入・検証する。kneed への適用（McGeer 1990b 原機械）は将来の後続とする。

## ゴールと完了ゲート

**円弧足 2D compass のリミットサイクルを数学的に発見し、安定性を固有値で証明する。**

完了ゲート（Phase 1-3 と同じ流儀）:

1. Newton shooting で不動点を発見（収束履歴をログ）
2. 固有値で漸近安定性（max|λ| < 1）を確認
3. McGeer 1990a の公表数値との一致（provenance 付き。図読取値は ±X% の正直ラベル）
4. **R→0 退化ゲート**: Phase 2 検証済み点足 compass の不動点
   (0.27103, −1.09238, −0.37737) と固有値に一致
5. ぽんぽこ殿の歩容判定: walk.mp4 が「円弧足で転がりながら歩く」ように見える

## 検証アンカー（文献）

**McGeer, T. (1990) "Passive Dynamic Walking", IJRR 9(2):62–82.**
（Cornell / Ruina ラボ公開 PDF、2026-06-16 取得。スキャン PDF のためページ画像から読取）

§5「Steady Walking of a General 2D Biped」の**直脚・円弧足モデル**が本フェーズの対象。
本セッションで読み取った数値（実装時の文献タスクで再検証・provenance 固定する）:

- 足形状 = 半径 R の円弧足（"round feet"）。`R/L` が鍵。R→0 で点足に縮退
- パラメータ: 足半径 R、脚の回転半径 r_gyr（分布質量）、脚質量中心位置、hip 質量、脚長 L、斜面 γ
- §8 / Table 1: test machine の S-to-S 固有値 ≈ 0.70（speed）, −0.65（swing）, −0.43（totter）、
  斜面 2.9%（≈ 1.66°）。step period / pendulum period 比 ≈ 0.36–0.39
- 注意: Table 1 は test machine（膝あり機）の値の可能性があり、§5 直脚モデルの印刷値と
  区別が要る。実装時の文献タスクで「直脚円弧足モデルの不動点・固有値」を特定し、
  図のみの量は digitize ±X% でラベルする（Phase 3 と同じ規律）

数値はプランに直接書かない。Claude の記憶からの数値をゲートに使わない（親設計の流儀）。

## 物理モデル

- **状態**: Phase 2 compass と同じ断面座標 `y = (θ_st, θ̇_st, θ̇_sw)`、2リンク（stance / swing）。
- **脚**: 一般剛体 — 脚質量 m、hip→CoM 距離、慣性 I（または回転半径 r_gyr）、hip 質量 m_h、
  脚長 L、**足半径 R**、斜面 γ。Goswami の質点脚（足から a、hip から b）はこの特例。
- **連続相（単相）**: stance 脚が円弧足で転がる。曲率中心 C は地面から高さ R にあり、脚角 θ_st に
  応じて接触点 P が移動する（転がり: C は水平に −R·θ_st 方向へ）。この転がり拘束を**位置運動学に
  組み込んで** derive レイヤー（Euler-Lagrange）で T, V を導出する。
- **heel-strike**: swing 脚の円弧足が接地した瞬間に瞬間衝突（新接触点回りの角運動量保存）＋脚交換。
  R が接地角・歩幅（step 幾何）を変える。受理条件は swing 足高さ 0 の下降接地（compass と同型、
  ただし円弧足の接地条件に補正）。
- **R→0** で点足 compass に厳密縮退（退化ゲートの土台）。

## コード構成（既存資産を最大再利用）

| パス | 責務 |
|---|---|
| `src/crane/models/rocker_compass.py` | 新規。`_build()`（lru_cache で記号導出 + lambdify）、`RockerCompassParams`、`make_rocker_compass()` → 単相 `HybridModel`（roll 力学 + heel-strike + lift/project）。kneed.py と同じ作法 |
| `src/crane/references_mcgeer.py` | 新規。McGeer 公表値を provenance 付き記録（パラメータ変換根拠 + 固有値 + step period、図読取はラベル） |
| `src/crane/viz.py` | `animate_rocker` 追加（円弧足を弧で描画。stance 足の転がりを接触点移動で表現） |
| `scripts/walk_rocker.py` | 新規。Phase 3.5 デモ（walk_compass.py 骨格 + 不安定警告 + meta.json） |
| 再利用 | `stride.py`（相機械）、`search.py`（Newton+Armijo+固有値）、`derive/lagrange.py`・`derive/impact.py` |

## 検証戦略（三段防御、Phase 3 と同じ）

1. **不変量テスト**: 平衡（揃った静止で加速度ゼロ）・転がり相のエネルギー保存（drift < 1e-7）・
   heel-strike の KE 散逸 / 脚交換幾何 / 零速度不動点・R→0 で力学＆衝突が compass に一致
2. **退化ゲート（最強）**: R→1e-9 + 質点等価パラメータ → 全 stride 不動点・固有値が Phase 2
   検証済み compass (0.27103, −1.09238, −0.37737) に一致
3. **文献ゲート**: McGeer パラメータ（変換・provenance）でリミットサイクル存在、固有値が
   §5/Table 1 値と図読取許容内で一致、step period 比照合
4. **歩容判定**: walk.mp4（円弧足の転がり込み）をぽんぽこ殿が判定

## リスクと対策

| リスク | 対策 |
|---|---|
| 転がり接触の力学導出（新規の肝）のバグ | R→0 退化ゲート + 転がり相のエネルギー保存テストで二重検証 |
| McGeer パラメータ変換（r_gyr/CoM 規約）の誤り | Phase 3 と同じく図・式（慣性行列・重力項）の独立検算で provenance 固定。勝手に直さない |
| basin が点足より細く Newton が発散 | 先日入れた Armijo 条件が効く。R=0（= Phase 2 不動点）から R を増やす continuation で初期推測を供給 |
| Table 1 が膝あり機の値で直脚モデルと不一致 | 実装文献タスクで直脚円弧足モデルの印刷値を特定。無ければ「存在 + 安定性 + R→0 退化」を主ゲートにし、固有値は参考比較に縮退（判断は文献タスクの結果で下す） |

## スコープ外（やらないこと）

- kneed + rocker（McGeer 1990b 原機械）— 本フェーズ完了後の後続候補
- 物理エンジン（MuJoCo / Genesis）再現 — 親設計の Phase 4 のまま（スコープ外・記録のみ）
- rocker foot の3D化・能動制御
