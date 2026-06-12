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

## Phase 2: sympy 導出レイヤー + Compass Gait (McGeer 1990)
Phase 1 ゲート通過後に計画策定。

## Phase 3: Kneed Walker (McGeer 1990)
Phase 2 完了後。
