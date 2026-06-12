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

文献 O(γ) 近似値 (0.199529, −0.198983) との差: 1.5e-3（O(γ^(5/3)) 打ち切り誤差の範囲内）

**θ* ∝ γ^(1/3) スケーリング則**: γ ∈ {0.004, 0.006, 0.009, 0.012} で log-log slope = 0.3287（理論値 1/3、誤差 0.005）

デモ実行 (1% 摂動から 30 歩、収束):

```
stride 0:  deviation 1.8e-3
stride 29: deviation 6.6e-10   ← リミットサイクル収束
```

## セットアップ

```bash
uv sync
```

Python 3.12 固定。

## 使い方

```bash
# デモ歩行（walk.mp4 + phase_portrait.png + meta.json を data/runs/ に出力）
uv run python scripts/walk_simplest.py [--gamma 0.009] [--strides 30] [--perturb 0.01]

# テスト（19 tests）
uv run pytest
```

## ドキュメント

- 設計: `docs/2026-06-12-crane-design.md`
- フェーズゲート: `GOALS.md`

## Heron との関係

Heron (Genesis/MuJoCo ポート) で確認されたこと: 物理エンジン上でのランダム探索は
basin が薄すぎてリミットサイクル候補に到達できない。
Crane はその教訓から出発し、解析的にサイクルを求めてからシミュレーションに渡す逆順戦略を採る。
