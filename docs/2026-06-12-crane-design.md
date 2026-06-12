# Crane 設計ドキュメント

**作成日:** 2026-06-12
**ステータス:** 承認済み（ぽんぽこ殿レビュー待ち）
**前身プロジェクト:** Heron (`~/git/heron`)

---

## 背景と動機

Heron（2026-05）では、受動歩行機械を Genesis / MuJoCo のフルフィジカルエンジン上で
MAP-Elites により探索したが、**安定リミットサイクルへの到達には失敗した**。
最良個体でも flips=5（5歩で転倒）が上限だった。

Heron の核心的教訓:

1. **受動歩行の basin of attraction は「very small and thin, fractal-like」**（文献
   表現どおり）。フルフィジカルエンジン + 確率的探索で「偶然 basin に乗る」のは
   構造的に困難。
2. **数値だけでは「歩行」と「バネ前進」を区別できない**（Heron issue #13）。
   最終判定は人間の動画観察が必要。
3. McGeer のオリジナル機械が成立した条件（rocker foot + 質点集中 + knee latch +
   質量分布）は**全部揃って初めて歩く**。部分的な再現では歩かない。

Crane はアプローチを反転する。エンジンで探すのではなく、**文献で実際にリミット
サイクルが発見されてきた王道——hybrid dynamics の直接積分 + Poincaré 写像の
不動点探索——を実装する**。basin がどれほど細くても、ニュートン法は不動点に
直接到達できる。

## ゴール

**受動歩行のリミットサイクルを数学的に発見し、その安定性を固有値で証明する。**

物理エンジン（MuJoCo 等）での再現は次段階（Phase 4、本プロジェクトのスコープ外、
記録のみ）。

## フェーズ構成

| Phase | 内容 | ゲート（完了条件） |
|---|---|---|
| 0 | プロジェクト雛形（uv / ruff / src layout）+ solve_ivp イベント検出基盤 | 単振り子のイベント検出デモが動く |
| 1 | Simplest walker（手導出）+ shooting / 固有値解析 | 文献値一致 + 歩容判定 |
| 2 | sympy 導出レイヤー + Compass gait + continuation | simplest 再導出一致 + 文献値一致 + 歩容判定 |
| 3 | Kneed walker（4相 hybrid dynamics） | リミットサイクル発見 + 安定性 + 歩容判定 |
| 4 | （スコープ外・記録のみ）物理エンジン再現、basin 可視化 | — |

## 対象モデルの階段

単純→複雑の3段階。各段階で「リミットサイクル発見 + 文献値照合 + 歩容判定」の
ゲートを通過してから次へ進む。

1. **Simplest walker** (Garcia, Chatterjee, Ruina, Coleman 1998
   "The Simplest Walking Model")
   — 2 DOF、点質量 hip + 微小足質量。公表された不動点・固有値があり、
   数値機構の検証に最適。
2. **Compass gait** (McGeer 1990 "Passive Dynamic Walking")
   — 有限質量分布の膝なし2脚。slope continuation で gait family と
   period-doubling カスケードも観察。
3. **Kneed walker** (McGeer 1990 "Passive Walking with Knees")
   — ぽんぽこ殿の本来の研究対象。4相 hybrid dynamics:
   unlocked swing (3 DOF) → knee-strike (衝突、膝ロック) →
   locked swing (2 DOF) → heel-strike (衝突、脚交換)。

## 実装アプローチ: ハイブリッド（手導出 → sympy 自動導出）

- **Phase 1 (simplest walker) は論文から手導出**。コードが読みやすく、公表数値で
  「積分器・イベント検出・shooting・固有値解析」という数値機構そのものを検証できる。
- **Phase 2 以降は sympy 記号導出レイヤー**。Lagrangian → 運動方程式、
  角運動量保存 → 衝突写像を自動導出して lambdify。代数の手作業ミスが
  膝付きで最大化するのを防ぐ。
- **レイヤーの検証**: sympy 層で simplest walker を再導出し、Phase 1 の
  手書き実装と数値一致することを確認してから compass に進む。

## 中核の数値手法

1. **Hybrid dynamics 直接積分**
   - 連続相: `scipy.integrate.solve_ivp`（高精度 RK、`rtol/atol` 厳しめ）、
     event 関数で heel-strike / knee-strike を検出
   - 離散相: 衝突写像（heel-strike は角運動量保存、knee-strike も同様）
2. **Stride 写像 S(x)**: 「積分 → 衝突 → 脚交換」の合成。副作用なしの純関数
3. **Poincaré shooting**: S(x*) = x* をニュートン法で解く（Jacobian は数値微分）。
   収束履歴は必ずログに残す
4. **安定性解析**: DS の固有値。max |λ| < 1 で漸近安定。文献の固有値と照合
5. **Continuation**: 発見した不動点を初期推定にして slope 等のパラメータを
   少しずつ動かし、gait family を追跡

## 検証ゲート（テスト戦略）

テストカバレッジではなく**文献の公表数値との一致**を各段階のゲートとする。

| Phase | 照合対象 |
|---|---|
| 1 | Garcia 1998 の不動点（stance angle 等）と固有値 |
| 2 | McGeer 1990 / Goswami らの compass gait 数値 |
| 3 | McGeer 1990 (knees) のパラメータ表と gait 数値 |

照合に使う具体的数値は**実装時に原論文から取得する**（Claude の記憶からは
書かない。Heron 教訓 #2 の防波堤）。

加えて各 Phase 末に**歩容アニメーションのぽんぽこ殿判定**（「これは歩行か」）を置く。

## ディレクトリ構成

```
crane/
├── CLAUDE.md / README.md / GOALS.md / DESIGN.md
├── pyproject.toml          # numpy, scipy, matplotlib（sympy は Phase 2 で追加）
├── docs/                   # 設計・知見（本ドキュメント含む）
├── src/crane/
│   ├── models/             # simplest.py, compass.py, kneed.py
│   ├── stride.py           # stride 写像（積分 + イベント + 衝突の合成）
│   ├── search.py           # Newton shooting / 固有値解析 / continuation
│   ├── derive/             # sympy 導出レイヤー（Phase 2〜）
│   └── viz.py              # phase portrait / stick-figure アニメ mp4
├── scripts/                # 各 Phase の実行スクリプト
└── data/runs/              # 出力（git 管理外、timestamp ディレクトリ方式）
```

## 主要抽象

- **`Model`** (protocol): `dynamics(t, x) -> xdot`, イベント関数群,
  `impact(x) -> x'`。パラメータは frozen dataclass
- **`stride(model, x0) -> StrideResult`**: 軌跡付き1歩写像
- **`find_limit_cycle(model, x_guess) -> FixedPoint`**: ニュートン法、
  収束履歴ログ付き

## 役割分担

- **Claude**: 文献調査、導出、実装、数値検証、可視化生成 — 自律推進
- **ぽんぽこ殿**: 歩容アニメーションの観察判定、パラメータ感度の直感照合、
  研究的方向の最終判断

## Heron から引き継ぐ流儀

- uv パッケージ管理、ruff format/check、Python 3.12+ src layout
- dataclass + 型ヒント積極活用、ログを惜しまない
- 1タスク1イテレーション1コミット、Conventional Commits（英語）
- 知見・課題は GitHub issue 化（`nyasuto/crane` リポジトリ想定）
- data/runs/<timestamp>/ への出力、上書き禁止

## スコープ外（やらないこと）

- 物理エンジン（MuJoCo / Genesis）での再現 — Phase 4 候補として記録のみ
- MAP-Elites / QD 探索 — リミットサイクルが手に入ってから次段階で再訪
- 3D 化、能動制御、強化学習
- Web ダッシュボード等の重い可視化（matplotlib で足りる範囲に留める）

## リスクと対策

| リスク | 対策 |
|---|---|
| 論文の式の転記ミス | Phase 1 は公表数値との一致で機構ごと検証 |
| sympy 導出のバグ | simplest を再導出して Phase 1 実装と突き合わせ |
| kneed の4相イベント順序の複雑さ | 相機械を明示的に書き、相遷移を全部ログ |
| ニュートン法の初期推定が悪く収束しない | 文献の既知解近傍から開始 + continuation で広げる |
| 「数値上の不動点」が歩行に見えない | 各 Phase 末にぽんぽこ殿の動画判定ゲート |
