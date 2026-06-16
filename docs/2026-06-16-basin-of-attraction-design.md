# Phase 4a 設計ドキュメント: Basin of Attraction の可視化

**作成日:** 2026-06-16
**ステータス:** 承認済み（ぽんぽこ殿レビュー済み 2026-06-16）
**前段:** Phase 0〜3.6（5モデルのリミットサイクル発見済み）
**後続:** Phase 4b（物理エンジンへの種戻し、本ドキュメントのスコープ外）

---

## 背景と動機

本プロジェクトの原点は Heron の失敗にある。Genesis/MuJoCo のフル物理エンジン上で
MAP-Elites（確率的探索）を回したが、安定リミットサイクルに到達できなかった
（最良でも 5 歩で転倒）。理由は受動歩行の basin of attraction が
"very small and thin, fractal-like" で、ランダム探索では偶然 basin に乗れないこと。

Crane はアプローチを反転し、**解析的に不動点を求め固有値で安定性を証明してから
シミュレーションに渡す**戦略を採った。Phase 0〜3.6 でこの「種」——不動点の厳密値と
安定性固有値——を 5 モデル（simplest / compass / kneed / rocker-compass /
rocker-kneed）について手に入れた。

Phase 4a は原点に立ち戻る二段構えの第一段である。物理エンジンに行く前に、
**Crane の解析モデル上で basin そのものを地図にする**。これにより:

1. **ナラティブ**: 「basin は薄くフラクタル的」を実際の地図で定量的に裏付ける
   （＝Heron の教訓を数値で締め、Crane の解析アプローチの正当性を示す）。
2. **実用情報（Phase 4b への橋）**: 「不動点からどれだけズレた初期条件まで収束するか」の
   許容半径を測り、物理エンジンへ渡す初期条件の誤差バジェットを見積もる。
3. **研究的主張（主軸）**: **円弧足や膝が basin を広げるか**。McGeer は「rocker foot が
   basin を拡大して歩行を安定化する」と論じた。今 Crane には point-foot/rocker-foot ×
   compass/kneed が揃っており、Phase 3.5/3.6 でやった rocker 化の意味を
   basin 面積比較で直接回収できる。

これら三つは排他ではなく、一枚の basin 地図から全て取れる。優先順位は 主軸3 ＞ 副産物2、
ナラティブ1 も同時回収。

## ゴール

受動歩行リミットサイクルの basin of attraction を Poincaré 断面の 2D スライスとして
可視化し、**円弧足が basin を広げるという主張を5モデル統制比較で示す**。
あわせて basin が安定性固有値（既に証明済み）と整合することを内部ゲートで保証する。

## スコープ

### やること
- Poincaré 断面の各初期条件から stride 写像を前進反復し「収束 / 転倒 / 未決」を分類
- 5モデルを同一スライスで描き、basin 面積（収束率）を比較
- simplest は 2D 断面そのものを全体図として描画（文献のフラクタル basin 図に対応）
- `multiprocessing` でグリッド点を並列分類（数分/モデル）

### やらないこと（YAGNI）
- 物理エンジン（MuJoCo/Genesis）連携 — Phase 4b
- 3D basin のボリュームレンダリング
- matplotlib 以外の可視化依存の追加（並列化は標準ライブラリのみ、新依存なし）

## 対象モデル（統制ペア）

主張がクリーンに出るよう統制ペアで構成する:

| モデル | 断面 | スライス | 役割 |
|---|---|---|---|
| simplest (Garcia 1998) | 2D `(θ, θ̇)` | 全体図 | ナラティブ1・文献照合 |
| point-foot compass (Phase 2) | 3D | `θ̇_sw` 固定, `(θ_st, θ̇_st)` 掃引 | 円弧足効果の基準（膝なし） |
| rocker compass (Phase 3.5) | 3D | 同上 | vs point compass |
| point-foot kneed (Phase 3) | 3D | 同上 | 円弧足効果の基準（膝あり） |
| rocker kneed (Phase 3.6) | 3D | 同上 | vs point kneed |

統制ペアは**同一 half_widths** で掃引するので basin 面積が直接比較できる。

## アーキテクチャ

### 新規モジュール `src/crane/basin.py`

```python
CONVERGED, FELL, UNDECIDED = 0, 1, 2

def classify_ic(model, y0, fixed_point, *, max_strides=20, converge_tol=1e-3) -> int:
    """断面点 y0 から stride を最大 max_strides 反復。
    deviation<tol で CONVERGED / StrideError で FELL / 上限まで未収束で UNDECIDED。"""

@dataclass(frozen=True)
class BasinResult:
    grid: np.ndarray          # (res, res) の分類コード
    ax0_vals: np.ndarray
    ax1_vals: np.ndarray
    axes: tuple[int, int]     # 掃引した断面座標インデックス
    fixed_point: np.ndarray
    basin_fraction: float     # 掃引窓内の収束セル率（比較スカラー）
    model_name: str

def basin_slice(model, fixed_point, *, axes=(0, 1), half_widths, resolution,
                max_strides=20, converge_tol=1e-3, n_workers=None) -> BasinResult:
    """fixed_point 中心、axes 以外の断面座標を不動点値に固定して 2D 掃引。
    multiprocessing.Pool でグリッド点を並列分類。"""
```

### `viz.py` への追加
- `plot_basin(result, path)`: 分類グリッドを imshow で描画（収束/転倒/未決を色分け、
  不動点をマーカー表示、軸ラベルは断面座標）。

### 新規スクリプト
- `scripts/basin_map.py`: 1モデルの basin 地図 ＋ meta.json を出力。
- `scripts/basin_compare.py`: 5モデル一括、比較 montage ＋ basin_fraction 表を出力。

## データフロー

1. `find_limit_cycle(model, guess)` で不動点を取得（既存 search.py を再利用）。
2. スライス定義: simplest は `(θ, θ̇)` 全体。3D断面モデルは `axes=(0,1)`、
   `θ̇_sw` を不動点値に固定して `(θ_st, θ̇_st)` を掃引。範囲は不動点中心 ± half_widths。
3. 各グリッド点: 断面座標 → `model.lift` で全状態へ → `stride` を最大反復 → 分類。
4. グリッド集約、basin_fraction（収束セル/窓内総セル）算出、画像＋スカラー保存。

## 分類ロジック

既存 `stride(model, x)`（転倒時に `StrideError`）を流用:

```
x = model.lift(y0)
for k in range(max_strides):
    try:
        result = stride(model, x)
    except StrideError:
        return FELL
    if norm(project(result.x_end) - fixed_point) < converge_tol:
        return CONVERGED
    x = result.x_end
return UNDECIDED
```

- `max_strides=20`, `converge_tol=1e-3` を既定値とし、実装時に収束/転倒の分離が
  シャープになるよう調整。
- 収束セルは早期に CONVERGED で抜け、転倒セルは早期に FELL で抜けるので、
  1点あたりの平均反復は max_strides より十分小さい。

## 窓幅（自動拡張）

不動点中心に既定 half_widths から開始し、**収束領域が窓端に接していたら**
（basin を切り落としている兆候）窓を一段広げる軽い自動拡張を入れる。
統制ペアでは最終的な half_widths を揃えて basin_fraction を比較する。

## 検証ゲート

これまでの Phase と異なり「狙う単一公表値」が乏しい。誠実なゲートは四段:

1. **内部整合（固有値との接続・自動テスト）**:
   - 不動点自身は CONVERGED（deviation≈0）。
   - 十分遠方は FELL。
   - 不動点の十分小さい近傍球は全 CONVERGED ＝ Phase 1〜3.6 で証明済みの
     `max|λ|<1`（局所安定）と整合する。basin が固有値計算と矛盾しないことを保証。
2. **文献照合（simplest）**:
   - simplest walker の basin は公表がある（**Schwab & Wisse 2001 "Basin of
     Attraction of the Simplest Walking Model"**、Garcia 1998 にも basin 図）。
   - 実装時に一次資料を取り寄せ provenance 付きで `references*.py` に記録し、
     **形状の定性一致 ＋（あれば）basin 面積系の数値**をゲートにする。
   - **記憶からの数値はゲートに使わない（CLAUDE.md 厳守）。** 数値が一次資料から
     取れなければ定性一致と内部整合のみをゲートとする。
3. **モデル間の主張（軸3・方向性ゲート）**:
   - rocker の basin_fraction > point-foot の basin_fraction
     （McGeer の「円弧足が basin を広げ安定化する」主張の方向性）。
4. **ぽんぽこ殿の目視判定**: basin 画像が主張を支持して見えるか（観察フェーズ）。

## 出力

- 各モデル: `data/runs/<timestamp>/basin.png` ＋ `meta.json`
  （basin_fraction・許容半径・params・スライス仕様・分類カウント）。
- 比較: 統制ペアを並べた montage ＋ basin_fraction 表（README/GOALS に記録）。

## フェーズ命名

GOALS の従来 Phase 4「物理エンジン再現＋basin 可視化」を分割:
- **Phase 4a**: basin 可視化（本ドキュメント）。
- **Phase 4b**: 物理エンジンへの種戻し（後続、別 spec）。

## リスクと対策

| リスク | 対策 |
|---|---|
| kneed の basin 計算が重い | 中解像度（100〜150²）＋ multiprocessing 並列。早期 break で平均反復削減 |
| 窓が basin を切り落とす | 収束領域が窓端に接したら自動拡張。統制ペアは最終窓幅を揃える |
| 2D スライスが basin の代表になっていない | 断面の固定座標を不動点値に取る（不動点を必ず含む面）。複数スライス角も検討可 |
| 単一公表値が乏しくゲートが弱い | 固有値整合（自動）＋ Schwab/Wisse 照合 ＋ 方向性主張 ＋ 目視の四段で締める |
| converge_tol/max_strides 依存で分類が揺れる | 収束/転倒の分離がシャープになる値を実装時に較正、meta に記録 |
