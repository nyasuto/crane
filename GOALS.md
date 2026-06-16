# ゴール設定 (Crane)

設計は docs/2026-06-12-crane-design.md。各 Phase のゲートを満たすまで次へ進まない。

## Phase 0: 雛形と数値基盤 — 完了 (2026-06-12)

- [x] uv / ruff / pytest / src layout (Python 3.12 固定)
- [x] solve_ivp イベント検出が解析解 π/2 と 1e-8 で一致

## Phase 1: Simplest Walker (Garcia 1998) — 完了 (2026-06-13)

- [x] EOM・衝突写像の実装と原論文 PDF との照合（完全一致）
- [x] 文献値の provenance 付き記録 (references.py)
- [x] Poincaré shooting で long-period 不動点発見 (0.2003109, -0.1998325)、文献 O(γ) 導出値と 1.5e-3 で一致
- [x] 固有値で漸近安定性を確認 (max|λ| = 0.589 < 1)
- [x] short-period 不動点も発見 (0.1939374, -0.2038669)、不安定 (max|λ| = 4.0) — 論文の記述と整合
- [x] θ\* ∝ γ^(1/3) スケーリング則を再現 (log-log slope 0.3287、理論 1/3)
- [x] **ぽんぽこ殿の歩容判定**: data/runs の walk.mp4 が「歩行」に見える（2026-06-13 合格判定）

### Phase 1 で得た知見

- heel-strike の正しい受理条件は「θ < 0 かつ ġ > 0（swing 足が下降して接地）」。
  交差点では ẏ_sw = sin(θ)·ġ が成り立ち、θ<0 だけでは scuff-dip 出口の上昇交差を
  誤って受理する（short-period 歩容と long-period basin の 15/49 軌道を破壊していた）
- 文献値ゲートの許容誤差は O(γ^(5/3)) 打ち切り誤差を考慮した絶対 1.5e-3 が適正
- continuation は単調チェーン必須（距離順だと branch hopping、γ=0.012 で不安定枝に乗り移る実例）
- basin 境界近くでは Poincaré 写像の固有値が +1 を跨ぎ Newton が発散する（±2% 摂動で実測）

### 今後の課題（issue 化済み）

- [#1](https://github.com/nyasuto/crane/issues/1) Newton バックトラッキングに Armijo 条件（残差増加ステップの排除）
- [#2](https://github.com/nyasuto/crane/issues/2) walk_simplest.py: 不安定サイクル収束時の警告 (max|λ|>1、γ=0.02 で λ=5.4 の実例)
- [#3](https://github.com/nyasuto/crane/issues/3) period-doubling カスケード探索 (γ → 0.0151+、Garcia Sec. 5.2-5.3)

## Phase 2: sympy 導出レイヤー + Compass Gait (Goswami 点足系) — 完了 (2026-06-13)

- [x] HybridModel プロトコル refactor: stride/search を model 非依存化（挙動変更なし）
- [x] derive/lagrange.py: Euler-Lagrange 自動導出。二重振り子の教科書 EOM と 1e-15 一致、Noether 保存 6e-15（tests/test_derive_lagrange.py に恒久化）
- [x] derive/impact.py: 角運動量・衝突写像導出。pivot 乗り移り解析解 ψ̇=−0.5 再現
- [x] models/compass.py: Goswami 点足 compass gait を記号導出で実装。不変量テスト（平衡・swing 相エネルギー保存・衝突 KE 減少・ラベル交換）全 pass
- [x] **Garcia 退化ゲート**: compass(a=0, m=1e-9, m_h=1, l=1, g=1) が Phase 1 検証済み実装と一致 — dynamics 5e-9 / impact 8e-11 / 不動点 4e-10 / 固有値 0.58911 / slave 方向 multiplier 1e-10
- [x] Goswami 文献値 provenance 付き記録 (references_goswami.py): INRIA RR-2996 (HAL) から取得。EOM・衝突定式化は sympy 記号比較で厳密一致（差ゼロ）
- [x] **Compass リミットサイクル発見** (γ=3°): y\*=(0.27103, −1.09238, −0.37737)、T=0.7343s（文献 0.735、0.1%）、strike 角 0.2710（文献 0.271）、max|λ|=0.580 安定
- [x] **固有値の1歩/2歩決着**: 印刷値 0.332 は2歩合成写像の multiplier（0.580²=0.336、1.3% 一致で数値確定）
- [x] **period-doubling カスケード再現**: period-1 不安定化 4.40°（文献 4.37°、分解能 0.05° 内）、period-2 branch 4.40°–5.95° 追跡、period-2 不安定化 ~4.95°（文献の第2分岐 4.9°）。flip bifurcation（実固有値が −1 通過）確認
- [x] viz model 非依存化 + walk_compass.py: 30歩で deviation 6.7e-3 → 5.3e-10
- [x] テスト 41 本全 green (commit ac66731)
- [x] **ぽんぽこ殿の歩容判定**: walk.mp4 が「歩行」に見える（2026-06-13 合格判定）

### Phase 2 で得た知見

- 衝突写像の正しい定式化: pre-pinned / post-pinned パラメータ化の二重評価 + swap 代入。
  有限質量でも保存則が機械精度 (1e-15) で成立する
- 記号導出レイヤーの検証は「検証済み実装への退化ゲート」が最強（Garcia 退化）。
  EOM・衝突それぞれで独立に確認できるため、バグの切り分けが明確
- 文献固有値の写像定義（1歩 vs 2歩）は本文と図で混在しうる → 自前数値で確定する。
  Goswami の印刷値 0.332 は 2 歩写像の multiplier、1 歩写像では 0.580²=0.336 で 1.3% 一致
- period-2 探索の seed は分岐点近傍で 1e-3 摂動では period-1 に回帰、1e-2 が必要だった

### 今後の課題

- viz の脚長 l=1 固定の明記/一般化（compass は l=1m なので現状問題なし）
- search のバックトラッキングに Armijo 条件（issue #1 と同根、compass でも同じ）
- カオス領域 (γ>5°) の Lyapunov 指数 / period-4, 8 の追跡（Phase 2.5 候補）

## Phase 3: Kneed Walker (点足, Hsu Chen 2007) — 完了 (2026-06-15)

足形状 = 点足（ぽんぽこ殿決定 2026-06-13。rocker foot は Phase 4 候補）。
4相 hybrid: unlocked swing (3 DOF) → knee-strike → locked swing (2 DOF) → heel-strike。
状態は全相 6D 統一 `x = [θ_st, θ_th, θ_sh, θ̇_st, θ̇_th, θ̇_sh]`、断面は heel-strike 直後の 3D
`y = (θ_st, θ̇_st, θ̇_sw)`。

- [x] 相機械への refactor: PhaseSpec 列 + stride 一般化（simplest/compass は挙動不変、全 green）
- [x] 文献値の provenance 付き記録 (references_kneed.py): Hsu Chen 2007 MIT 修論 Ch.3 から取得。
      パラメータ・座標規約・衝突定式化を記号比較で照合（Table 3.1 と全パラメータ一致、heel-strike は
      原典 eq 3.5 = rigid-through-collision と一致）
- [x] derive レイヤーで両相力学を記号導出 (kneed.py): 平衡・unlocked/locked 相のエネルギー保存
      (drift < 1e-7)・locked 相の θ_th≡θ_sh 維持を不変量テストで確認
- [x] 衝突写像 (knee-strike / heel-strike): 位置不変・速度ロック・KE 散逸・脚交換幾何・剛体保持・
      零速度不動点を不変量テストで確認
- [x] **compass 退化ゲート**: m_s→1e-9 で locked 力学・heel-strike・全 stride 不動点・上位2固有値が
      Phase 2 検証済み compass (0.27103, −1.09238, −0.37737) に一致 — 記号導出レイヤーの最強検証
- [x] **Kneed リミットサイクル発見** (γ=0.0504 rad): y\*=(0.23859, −1.10959, −0.05715)、
      Newton 残差 5e-14（4反復）。q1=0.1882（文献 0.1877、0.27%）、q̇1=−1.1096（文献 −1.1014、0.7%）、
      step period 0.5643s（文献図読取 0.559、1.0%）、heel-strike 直前 stance 角 −0.2386（文献 −0.2380、0.25%）
- [x] **安定性**: Jacobian multiplier 0.655, 0.655, 0.144（複素対 0.353±0.552i + 実 0.144）、
      max|λ|=0.655 < 1 で漸近安定。文献 Sec.3.3.1 の locally stable 記述と整合
- [x] gait family continuation: 公表 slope の 0.8 倍まで追跡
- [x] viz 4セグメント + walk_kneed.py: 12歩で deviation 9.7e-3 → 1.5e-4、膝屈曲 ~30°/stride 観察可能
- [x] テスト 63 本全 green
- [x] **ぽんぽこ殿の歩容判定**: walk.mp4 が「歩行」に見える（膝の屈曲込み）。
      「歩行が素晴らしく自然で、膝を使っている」と合格判定（2026-06-15）

### Phase 3 で得た知見

- 文献固有値が「Jacobian 有限差分」でなく「摂動軌道の最小二乗推定」のことがある（Hsu Chen は
  25 軌道×10s の least-squares で 0.405 を得た）。本実装の Jacobian 値 0.655 とは差があるが、
  原典が「参考値」と位置づけており両者とも安定。固有値ゲートは自前数値で確定する流儀が正しい
- 零近傍の速度成分（θ̇_sw: 本実装 −0.057 vs 文献 −0.040）は固定刻み前進積分（原典 dt=1e-3）の
  バイアスが出やすい。確信を持って印刷された量（q1, q̇1, interleg 角, step period）は全て <1% 一致
- heel-strike の 3 本目（新 swing 脛速度）は rigid-through-collision (θ̇_sh⁺=θ̇_th⁺) が McGeer の
  機械式ロック物理と整合。当初案の「脛を knee 回りで独立保存」は不採用（references_kneed.py 参照）
- 退化ゲート (massless shank) は記号導出レイヤーの検証として極めて有効。knee-strike 撃力→0、
  脛が slave 振り子になり、Phase 2 検証済み compass に厳密一致する

### 今後の課題

- rocker foot 化（McGeer 1990 原機械、Phase 4 候補）
- search のバックトラッキングに Armijo 条件（issue #1 と同根）
- basin 可視化 → Phase 4a で実施（完了）。物理エンジン (MuJoCo) seed 受け渡し → Phase 4b（未着手）

## Phase 3.5: Rocker-foot Compass (McGeer 1990a) — 完了 (2026-06-16)

円弧足（半径 R の round foot）2D compass を記号導出レイヤーで実装。剛体脚を
2点質量（m/2 を hip 距離 c±ρ に配置）で表現することで質量・CoM・慣性を再現し、
既存の点質量 derive レイヤーを**そのまま再利用**。stance 円弧足は転がり接触
（曲率中心 C=[−R·θ_st, R]）。R→0 で点足 compass に退化する。

- [x] 文献値の provenance 付き記録 (references_mcgeer.py): McGeer, T. (1990)
      "Passive Dynamic Walking", IJRR 9(2):62–82, §5 straight-leg round-foot model。
      パラメータ m=1.0, m_h=0.0, c=0.37, ρ=0.32, L=1.0, R=0.3, γ=0.030 rad（図読取 ±10%）
- [x] derive レイヤーで round-foot 力学を記号導出 (dynamics)・衝突写像 (impact)・
      4D stride を実装。不変量テスト（エネルギー保存・転がり接触・脚交換）全 pass
- [x] **R→0 退化ゲート（最強検証）**: R→1e-9, ρ→1e-9, c=b で全 stride 不動点・固有値が
      Phase 2 検証済み点足 compass (0.27103, −1.09238, −0.37737) に一致（不動点 ~1e-10 / 固有値 ~1e-8）
- [x] **Rocker リミットサイクル発見** (γ=0.030, R=0.3): y\*=(0.30844, −1.26256, −0.87914)、
      references SECTION_GUESS から Newton 収束。step period 0.8417s（McGeer 2.5·√(l/g)=0.7982s、5.45%、±10% 内）
- [x] **安定性**: 固有値 0.4316 と複素対 −0.1529±0.2986j（大きさ 0.4316, 0.3355, 0.3355）、
      max|λ|=0.4316 < 1 で漸近安定。McGeer Fig.9 と整合（公称質量配置は点足では不安定、R=0.3 の円弧足で安定化）
- [x] **独立傍証（ゲート外）**: strike 時の half-interleg 角 |θ_st|=0.308 ≈ McGeer 印刷値 α₀=0.30
- [x] R-continuation: gait family を R=0.30→0.24→0.18→0.12 まで追跡、全て収束・安定
- [x] walk_rocker.py: 30歩で deviation 3.95e-3 → 8.5e-14（リミットサイクル収束）
- [x] テスト 77 本全 green
- [x] **ぽんぽこ殿の歩容判定**: walk.mp4 が「円弧足で転がりながら歩く」ように見える（2026-06-16 合格判定）

### Phase 3.5 で得た知見

- 円弧足の静的平衡角 ≠ γ（接触点 P_st=[−R·θ,0] が θ とともに動くため）。R→0 で初めて γ に近づく。
  点足モデル由来の素朴な「θ=γ が不動点」チェックは転がり足では誤り
- 剛体脚を 2 点質量（m/2 を CoM±ρ に配置）で表すと質量・CoM・慣性を再現でき、
  既存の点質量 derive レイヤーを**改変なし**で流用できる。回転慣性を加えるクリーンな手法
- McGeer の印刷固有値（Table 1: 0.70/−0.05/−0.83）は kneed テスト機械のもので、§5 straight-leg
  モデルの値ではない（後者は Fig.9 の locus のみで数値固有値の印刷なし）。
  ゆえに固有値ゲートは reference-only (POINCARE_EIGENVALUES=None)、ゲートは存在性・安定性・step period・R→0 退化

### 今後の課題

- kneed + rocker = McGeer 1990b 原機械の完全再現（Phase 3 + Phase 3.5 の合流、次候補）
- search のバックトラッキングに Armijo 条件（issue #1 と同根）
- basin 可視化 → Phase 4a で実施（完了）。物理エンジン (MuJoCo) seed 受け渡し → Phase 4b（未着手）

## Phase 3.6: Rocker-foot Kneed Walker (McGeer 1990b 原機械) — 完了 (2026-06-16)

Phase 3 点足 kneed walker（kneed.py、4相 hybrid）に Phase 3.5 の円弧足（半径 R の
round foot 転がり接触）を合流させ、McGeer 1990b "Passive Walking with Knees" の原機械を
点質量で再現。`models/rocker_kneed.py` = kneed.py の \_build() に **5 箇所のみ**の
外科的変更（転がり hip + 転がり接触の衝突 pivot）を施したもの。derive レイヤー /
kneed.py / rocker_compass.py は一切改変せず。

- [x] 文献値の provenance 付き記録 (references_mcgeer_knees.py): McGeer, T. (1990)
      "Passive Walking with Knees", Proc. IEEE ICRA（6ページ会議論文）。
      パラメータ m_h=0.6885（hip 質量 68%）, m_t=0.1, m_s=0.062, l_t=0.46, l_s=0.54,
      b_t=0.20, b_s=0.24, R=0.20, γ=0.0456 rad
- [x] kneed.py の \_build() に 5 箇所の変更で rocker 化（転がり hip C_st=[−R·θ_st, R];
      stance 接触 P_st; swing 円弧足接触 C_sw=knee_sw+(l_s−R)·down(θ_sh);
      knee-strike & heel-strike の系 pivot [0,0]→P_st と foot_sw→P_sw; パラメータ R）。
      不変量テスト（references・dynamics・impact・stride・cycle・viz）全 pass
- [x] **R→0 退化ゲート（最強検証）**: R→1e-9 で全 stride 不動点・固有値が
      Phase 3 検証済み点足 kneed (0.23859, −1.10959, −0.05715) に一致（不動点 ~1e-10 / 固有値 ~1e-8）。
      unlocked/locked 力学・knee-strike・heel-strike の 4 経路すべてで一発成立
- [x] **Rocker-kneed リミットサイクル発見** (McGeer config): y\*=(0.35661, −1.20906, −0.54601)、
      references SECTION_GUESS から Newton 直接収束。膝屈曲 ≈ 27.6°
- [x] **安定性**: 固有値 0.39527±0.47248i と実 0.04579（大きさ 0.616, 0.616, 0.046）、
      max|λ|=0.616 < 1 で漸近安定。McGeer の安定サイクルと整合
- [x] **step period**: t_step=0.7945s（McGeer 2.7·√(l/g)=0.862s、7.84%、±15%（GAIT_TOLERANCE）内）。
      McGeer 値は無次元（√(l/g)）なので換算
- [x] walk_rocker_kneed.py: 30歩で deviation 5.93e-3 → 4.87e-9（リミットサイクル収束）
- [x] テスト 89 本全 green（Phase 3.6 分: references・dynamics・impact・stride・R→0 退化ゲート・cycle・viz）
- [x] **ぽんぽこ殿の歩容判定**: walk.mp4 が「円弧足で転がりつつ膝を使って歩く」ように見える（2026-06-16 合格判定）

### Phase 3.6 で得た知見

- 検証済みの 2 サブシステム（Phase 3 kneed + Phase 3.5 円弧足）の合流は、kneed.py の
  \_build() への **5 箇所**の外科的変更のみで済んだ（転がり hip C_st=[−R·θ_st, R]; stance 接触 P_st;
  swing 円弧足接触 C_sw=knee_sw+(l_s−R)·down(θ_sh); knee-strike/heel-strike の系 pivot を
  接触点へ; パラメータ R）。R→0 退化が 4 経路すべてで一発成立したのは合流設計の健全性の証左
- 本実装の kneed は点質量モデル（Phase 3 Hsu Chen と同様）だが、McGeer 1990b の実機械は
  分布質量脚（回転半径 ≠ 0）。ゆえにこれは McGeer 原機械の**点質量再現**であり、
  存在性・安定性・step period・half-interleg は一致するが、印刷固有値は一致しない（reference-only）
- McGeer 1990b は（1990a と違い）固有値を印刷しているが、4D 断面のもので本実装の 3D 断面とは
  次元が異なる。よって配列の直接比較は不可。我々 max|λ|=0.616 vs McGeer 印刷値 0.447 の差は
  分布質量 vs 点質量 + 断面次元差に由来する documented な reference check（ゲート外）

### 今後の課題

- 分布質量脚（回転半径 ρ ≠ 0）で McGeer 1990b の印刷固有値（0.447）を厳密一致させる
  （Phase 3.5 の 2 点質量手法を kneed の各脚セグメントに適用すれば可能）
- search のバックトラッキングに Armijo 条件（issue #1 と同根）
- basin 可視化は Phase 4a で実施（下記）。物理エンジン (MuJoCo) への seed 受け渡しは Phase 4b（未着手）

## Phase 4a: Basin of attraction 可視化 — 完了 (2026-06-16)

各 Poincaré 断面点から stride 写像を最大 20 反復し「収束 / 転倒 / 未決」を分類して
basin（吸引域）を地図化（`basin.py`）。統制ペア（点足 vs 円弧足）で **円弧足が basin を
広げるか** を `basin_fraction`（窓内収束セル率＝面積の代理）で検証。窓はペア内で同一に
適用し、仮説が偽でも正直に報告する。

**較正した窓**: 4つの 3D 断面モデル（compass / rocker_compass / kneed / rocker_kneed）は
公平比較のため**同一窓 half_widths=(0.16, 0.65)**（θ_st 軸 ±0.16、θ̇_st 軸 ±0.65、不動点中心）を
共有。res=41 較正で収束領域は概ね窓内に収まり、残る漏れは +θ_st / -θ̇_st 方向の細い舌状部の
先端 0〜3 セルのみ（完全包含は可視構造を潰すため不採用。漏れは basin が大きい側ほど大きく、
両判定を弱めず強める向き）。simplest は文献ゲートが固定する 2D 窓 (0.08, 0.08)。
本番は **解像度 R=60、14 並列、≈35 分**（R=80 は ≈17 分/モデルで基準超のため R=60）。

**実測 basin_fraction (R=60)**:

| モデル | 足形状 | basin_fraction |
|--------|--------|----------------|
| simplest | 点足 | 0.029 |
| compass | 点足 | 0.061 |
| rocker_compass | 円弧足 | 0.117 |
| kneed | 点足 | 0.108 |
| rocker_kneed | 円弧足 | 0.061 |

**円弧足 vs 点足の判定（正直に）**:
- **Compass ペア**: rocker 0.117 > point 0.061（相対 +92%）→ 円弧足が basin を約 1.9 倍に拡大。robust、仮説支持。
- **Kneed ペア**: rocker 0.061 < point 0.108（相対 −43%）→ 円弧足が basin を約 0.57 倍に縮小。robust に**仮説と反する**。
  較正時に点足 kneed の方が窓端漏れが大きく、clipping アーティファクトではない。failure ではなく
  **要調査の finding**（膝つき機構と円弧足の相互作用が compass 系と異なる可能性）。

### ゲート
- [x] basin.py 分類 + 並列スライス + plot_basin（テスト pass）
- [x] 内部整合ゲート（不動点近傍球が全 CONVERGED ＝ max|λ|<1 と整合）
- [x] simplest 文献ゲート（Schwab & Wisse 2001 provenance、定性一致：薄くフラクタル的）
- [x] 5モデル統制比較（calibrated 窓 (0.16,0.65)・解像度 R=60）：実測 basin_fraction を記録
      （simplest 0.029 / compass 0.061 / rocker_compass 0.117 / kneed 0.108 / rocker_kneed 0.061;
      compass ペアは rocker>point で仮説支持、kneed ペアは rocker<point で仮説に反する＝要調査 finding）
- [x] **ぽんぽこ殿の目視判定**: basin_compare.png が実測と整合（rocker_compass は compass より明らかに広い、
      rocker_kneed は kneed より狭い、simplest は最も狭く薄い）と確認（2026-06-16 合格判定）

### Phase 4a で得た知見

- 「円弧足が basin を広げる」は普遍則ではない。compass ペアでは成立（×1.9）するが、
  kneed ペアでは反転（×0.57）。膝つきハイブリッド（4相）と円弧足転がり接触の相互作用が
  compass 系（2相）と質的に異なる可能性。Phase 4b 以降の調査対象
- **文献照合（一次資料、2026-06-16 取得）**: この結果は文献と整合する。McGeer 1990b
  ("Passive Walking with Knees", ICRA) は膝を必須ではなくクリアランス装置と位置づけ
  （遊脚の "toe stubbing" を受動的に解決；直脚機は 3 J/step のリトラクションモータが必要だった）、
  安定性向上は "in some cases more stable" と控えめ。basin の拡大は**膝ではなく円弧足**に帰属
  （Hsu Chen 2007 MIT thesis §3.3.2:「curved feet ... increase their basins of attraction」）。
  さらに Hsu Chen は「足擦り＝転倒を課すと kneed の basin は縮む」と示す。kneed vs 直脚の
  basin **サイズ**を定量比較した査読論文は見当たらず、本結果は文献ギャップに位置する
- **重要 caveat（issue #6 連動）**: Hsu Chen の「kneed basin 縮小」は「曲がったまま接地＝転倒」を
  課した場合の話。本実装の basin 計算がこの bent-knee 接地を `StrideError` で転倒検出できているか
  未確認（= 未実装の issue #6）。未検出なら kneed の basin_fraction は過大評価で、真の物理 basin は
  さらに小さい可能性。issue #6 実装後に kneed basin を測り直すのが次の一手
- simplest の basin は窓を斜めに横切る非常に薄いフラクタル的スリバー（fraction 0.029）。
  Schwab & Wisse 2001 の "very small and thin, fractal-like" と定性一致。物理エンジン上の
  ランダム探索が届かない（Heron 教訓）理由を可視化で裏付け
- basin_compare.png / basin_fractions.json は `data/runs/` 配下＝**gitignore 対象のローカル成果物**:
  `data/runs/20260616_144724_basin_compare/basin_compare.png`

### 今後の課題（Phase 4b 以降）

- **Phase 4b（未着手）**: 解析的に求めた不動点を初期値として物理エンジン (MuJoCo 等) に
  受け渡す逆順戦略の検証（Heron が失敗したランダム探索を解析 seed で置き換える）
- kneed ペアで円弧足が basin を縮小する機構の解明（4相 hybrid × 円弧足の相互作用）
- [#19](https://github.com/nyasuto/crane/issues/19) 受動歩行→車輪の連続軸 → **Phase 4a.1 で実施（下記）**

## Phase 4a.1: 受動歩行→車輪の連続軸（R 連続軸、issue #19）— 完了 (2026-06-16)

「コンパスは車輪への漸近」という概念枠組み（rimless wheel が最簡受動歩行モデル）を、
rocker_compass の円弧足半径 R を固定勾配 γ=0.030 で R∈(0, L=1) に掃引して検証。
`hip = C − (L−R)·down(θ)` より R→L で hip が足中心に乗る＝半径 L の車輪の極限。
R=0.3 起点に continuation で両方向追跡。主指標は相対損失 δ(R)=(KE直前−KE直後)/KE直前。

**実測（ΔR=0.05、basin res=50、12 並列）**:

| R | δ（効率損失） | max\|λ\| | basin_fraction |
|---|--------------|----------|----------------|
| 0.05 | 0.303 | 0.437 | 0.037 |
| 0.20 | 0.292 | 0.373（最安定） | 0.072 |
| 0.40 | 0.272 | 0.517 | 0.200 |
| 0.60 | 0.241 | 0.703 | 0.430（最大） |
| 0.70 | 0.210（最効率） | 0.867 | 0.026（崩壊） |
| 0.75 | — | — | 歩行消失 |

### ゲート
- [x] efficiency.py（δ・衝突損失・COT）+ R=0.3 アンカー + 連続性ゲート（テスト pass）
- [x] R=0.3 アンカー: 掃引が Phase 3.5 検証済みサイクル (0.30844, −1.26256, −0.87914) を 1e-3 で再現
- [x] R 掃引: δ(R)・max|λ|(R)・basin_fraction(R) を測定（上表）
- [x] エネルギー収支チェック: 機械的 COT(R) = 衝突損失/(総質量·g·一歩水平距離) が全 R で sin γ=0.030 に **7桁一致**
- [x] rimless wheel provenance（McGeer 1990a を Tedrake "Underactuated Robotics" Ch.4 で確認、references_mcgeer に記録）
- [x] **ぽんぽこ殿の目視判定**: R_sweep.png が「車輪への漸近」（δ↓）と効率/安定性/頑健性トレードオフを支持して見える（2026-06-16 合格判定）

### Phase 4a.1 で得た知見
- **3指標が別々の R で最適化される三つ巴のトレードオフ**: 効率 δ は R↑ で単調改善（H3 支持）、
  安定性 max|λ| は R≈0.2 で最小、basin（頑健性）は R≈0.6 でピーク後 R=0.70 で崩壊、
  R=0.75 で歩行消失（**純粋な車輪には到達できない**、H2 支持）。H1（basin 単調増）は**不成立**。
- 「車輪が最も効率的だが、歩行を維持するには車輪極限から離れる必要がある」。ぽんぽこ殿の
  「半径を上げれば良いが極端は取れない」という直感が 3 指標の分離として定量化された。
- 機械的 COT がリミットサイクル上で sin γ に固定される（重力入力＝衝突損失）ことを 7 桁で実証。
  COT は効率の判別力を持たず、判別には相対損失 δ が必要——という当初の物理的見立ても確認。
- **歩行消失は fold（鞍点-ノード）分岐**: 安定枝の支配的固有値は実で正、R↑ で単調増加し
  R=0.725 で +0.960、R≈0.75 で安定枝が消失（0.30→0.432, 0.60→0.703, 0.70→0.867, 0.725→0.960）。
  実固有値が +1 を通過＝安定サイクルが共存する不安定サイクル（鞍点）と対消滅する fold の指紋。
  実際 R=0.6 付近に max|λ|≈5.77 の不安定サイクルが共存（continuation の誤シードで偶然観測）。
- **安定性と basin が一致しない物理的理由（局所 vs 大域）**: max|λ| は不動点での写像の傾き＝
  局所回復速度（井戸の急峻さ、R≈0.2 で最良）。basin の大きさは共存不安定鞍点の安定多様体＝
  basin 境界の位置で決まる大域量（井戸の広さ、R≈0.6 で最大）。支配する対象が異なるため同じ R で
  最適化されない。R→fold で鞍点が安定サイクルへ接近＝basin 境界が掃き寄せられ、局所固有値が 1 に
  達する前に basin が急崩壊（0.43→0.026）。McGeer の「固有値は回復速度を与えるが許容摂動の大きさは
  語らない」（井戸の深さ≠広さ）を定量化・機構同定したもの。Heron 教訓（ランダム探索は大域 basin に
  着地せねばならず局所固有値だけでは不十分）の数理的裏付けでもある。
- 成果物（gitignore 対象）: `data/runs/20260616_225010_wheel_limit/`（R_sweep.png・basin_R_montage.png・R_sweep.json）。
- **新規性評価（2026-06-16、`docs/2026-06-16-phase4a1-novelty-assessment.md`）**: 一次資料サーベイの結果、
  本 Phase の概念的中身はほぼ既出。効率↑（Adamczyk/Collins/Kuo 2006）、安定性≠basin（**Schwab & Wisse 2001**
  ——我々が simplest basin で引用済みの同論文が「basin と線形安定性に直接の関係は無い／最も頑健な設計は
  basin 最大の方」と明言）、円弧足の fold＋共存鞍点＋basin 崩壊（**MMT 2025/26**）。決定的テストで GSN 2007 /
  MMT 2025/26 が**足半径を制御軸として既に使用**と判明し、「足半径という新軸」の擁護線も弱い。残りうる狭い
  貢献は (1) 効率×固有値安定性×basin の3軸統合提示、(2) basin の内部最大が単調傾向と対立しうる点のみ（要
  paywall 全文確認）。F（basin＝淘汰の頑健性）は文献に無く我々の仮説。**firm な新規性主張は現状の証拠では
  支持されない**。位置づけは「既知結果の円弧足 compass での再現・統合」。

## Phase 5a: 動力付き simplest walker（Kuo 2002）— 完了 (2026-06-17)

simplest walker（Garcia 1998）に pre-emptive push-off（後脚軸方向の撃力 P）を足し、
**受動では不可能だった平地 γ=0 歩行**を実現した（能動歩行アークの第一段）。連続相は
simplest と同一（γ 経由）で、heel-strike だけを「push-off 撃力 P → 角運動量保存衝突 →
脚交換」の合成写像に差し替える（P=0 で受動に厳密退化）。平地 γ=0 では θ 部分系が push-off と
分離し素朴な shoot が静止解に落ちるため、受動 γ=0.009 から γ→0・push_off→target へ
continuation して平地サイクルを発見。

### ゲート
- [x] powered_simplest.py（push-off 合成衝突写像、P=0 で受動退化）
- [x] **push-off→0 退化ゲート**: P=0 で Phase 1 検証済み受動サイクル (0.2003109, −0.1998325) に全 stride 一致
- [x] 平地 γ=0 リミットサイクル発見 + 安定性 max|λ|<1（continuation で探索; 実測 y*=(0.3266,−0.3395), P=0.115, max|λ|=0.650）
- [x] エネルギー収支: push-off 仕事 P²/2 = 一歩衝突損失（~1e-13）
- [x] Kuo 2002 provenance 記録（abstract 逐語確認 "toe-off は 1/4 コスト" + power law、本文 paywall を明記）
- [x] **ぽんぽこ殿の目視判定**: 平地 walk.mp4 が「push-off で前進する歩行」に見える（2026-06-17 合格判定）

### Phase 5a で得た知見
- **受動には不可能だった平地歩行を pre-emptive push-off で実現**（能動歩行アークの第一段）。
  受動歩行は勾配の重力でしか駆動できないが、heel-strike 前の push-off 撃力で COM を redirection し、
  平地 γ=0 でも安定リミットサイクル（max|λ|=0.650）を成立させた。
- 平地 γ=0 では θ 部分系が push-off と分離するため、素朴な shoot は静止解に落ちる。受動サイクルから
  γ→0・push_off→target へ **continuation 必須**。25 歩シミュで deviation が全体として減衰（4.0e-3→5.9e-7、
  初期は振動）して平地サイクルへ収束することを確認。
- 成果物（gitignore 対象）: `data/runs/<timestamp>_powered_P0.115/`（walk.mp4・phase_portrait.png・meta.json）。
- 続き **Phase 5b（能動 rocker_compass で「受動で良い個体は能動でも強いか」検証）** を今後の課題に。

## Phase 5b: 動力付き rocker_compass（push-off 増強）— 完了 (2026-06-17)

円弧足 rocker_compass（Phase 3.5）に pre-emptive push-off（後脚軸方向の撃力）を足し、
固定勾配 γ=0.030 で push-off がエネルギーを注入する様子を観察した。連続相は rocker_compass と
同一で、heel-strike を「push-off 撃力 → 角運動量保存衝突 → 脚交換」の合成写像に差し替える
（push_off=0 で受動に厳密退化）。

### ゲート
- [x] powered_rocker_compass.py（push-off モーメントを eq1 に追加、push_off=0 で受動退化）
- [x] **push-off→0 退化ゲート**: push_off=0 で Phase 3.5 検証済み (0.30844, −1.26256, −0.87914) に全 stride 一致
- [x] 動力サイクル γ=0.030・push_off∈[0,0.08] で安定（max|λ| 0.432→0.462、実測値記録）
- [x] push-off がエネルギー注入（post-collision KE 増加＝テストでゲート化。walkable slope の平地側拡張は de-risk 継続で確認）
- [x] **ぽんぽこ殿の目視判定**: walk.mp4 が「円弧足で転がりつつ push-off で歩く」に見える（2026-06-17 合格判定）

### Phase 5b で得た知見
- rocker は push-off でも完全な平地 γ=0 には届かず fold でサイクルが消失（R=0.3）。Phase 4a.1 の
  車輪極限 fold と響き合う（平地歩行のデモではなく記録された知見）。
- push-off 公式は後脚軸撃力が hip を通るため eq1 のみにモーメントを追加（実装前 de-risk で検証）。
- 動力サイクル γ=0.030, push_off=0.08 は y*=(0.29247,−1.20217,−0.85407), max|λ|=0.462 で安定。
  不動点近傍からの 20 歩シミュで deviation が単調減衰（3.6e-3→1.3e-9）して収束。
- 成果物（gitignore 対象）: `data/runs/<timestamp>_powered_rocker_P0.08/`（walk.mp4・phase_portrait.png・meta.json）。
- 続き **Phase 5c（R 掃引で能動 vs 受動 basin 比較＝本命仮説）** を今後の課題に。

