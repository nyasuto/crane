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
- [x] θ* ∝ γ^(1/3) スケーリング則を再現 (log-log slope 0.3287、理論 1/3)
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
- [x] **Compass リミットサイクル発見** (γ=3°): y*=(0.27103, −1.09238, −0.37737)、T=0.7343s（文献 0.735、0.1%）、strike 角 0.2710（文献 0.271）、max|λ|=0.580 安定
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
- [x] **Kneed リミットサイクル発見** (γ=0.0504 rad): y*=(0.23859, −1.10959, −0.05715)、
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
- 物理エンジン (MuJoCo) 再現と basin 可視化（Phase 4、設計上スコープ外・記録のみ）
