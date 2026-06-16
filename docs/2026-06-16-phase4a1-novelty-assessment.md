# Phase 4a.1 新規性評価（受動歩行→車輪の連続軸）

**作成日:** 2026-06-16
**対象:** Phase 4a.1（issue #19、PR #21 マージ済み）の結果と、ぽんぽこ殿の考察の新規性
**方法:** 一次/査読資料の文献サーベイ（3 並行調査エージェント＋決定的全文取得テスト）。
記憶からの断定はせず、取得元 URL と取得形態（全文 / abstract / 検索要約）を明記。
**重要な但し書き:** 「文献が見つからない」ことは新規性の証明ではない。本評価は
「文献で覆われている／覆われていないように見える」境界の暫定報告であり、
firm な新規性主張には下記の paywall 全文2本の確認が必要。

---

## 我々の主張（6点に分解）

- **A**: 円弧足半径 R↑ で効率↑（衝突損失↓、車輪極限で 0）。
- **B**: 円弧足が basin を広げる。
- **C**: 局所安定性（固有値）と basin サイズが**別の R で最適化される**（乖離）。
- **D**: 足半径→車輪極限で歩行枝が fold（鞍点-ノード、支配固有値→+1）で消滅、basin は分岐前に崩壊。
- **E**: 効率×局所安定性×basin の**3軸トレードオフ**（各々別の R で最適）。
- **F**: basin（有限外乱耐性）こそ生物/淘汰に効く頑健性。

## 我々の実測（参照）

R∈(0.05, 0.70) 掃引、γ=0.030: δ 単調減少（0.303→0.210）; max|λ| は R≈0.2 で最小（0.373）、
R=0.725 で 0.960、R≈0.75 で安定枝消失; basin は R≈0.6 で最大（0.430）後 R=0.70 で崩壊（0.026）;
支配固有値は実で正、+1 へ（fold）; 機械的 COT ≈ sin γ を 7 桁一致。

---

## 主張別 新規性評価

| 主張 | 判定 | 最も近い先行研究（取得形態） |
|---|---|---|
| A | **既出・確立** | Adamczyk, Collins & Kuo 2006（衝突仕事 ∝ (1−ρ)²、ρ→1 で 0）; Kuo 2002; McGeer 1990a |
| B | **既出** | McGeer 1990a; GSN 2007（足半径を外乱除去ノブに使用）。一次の「basin 面積増加」明示文は未取得 |
| C | **概念は既出（2001）。足半径軸も GSN が先行** | Schwab & Wisse 2001 Fig 2.8（全文取得・逐語）; GSN 2007（足半径×Floquet-vs-外乱、abstract） |
| D | **機構は既出（2025/26）。足半径も同論文の軸** | MMT 2025/26「Falling dynamics」（SN fold＋共存鞍点＋crisis、foot radius×slope、要約） |
| E | **対は既出、3軸統合は未発見（ギャップ候補）** | Kwan & Hubbard 2007（2軸: 速度/効率 vs 局所安定性、abstract） |
| F | 工学側は既出／**生物・淘汰側は我々の解釈（文献に無し）** | Hobbelen & Wisse Limit Cycle Walking; Byl & Tedrake 2009; Heim & Spröwitz 2018（basin 批判） |

### Claim C の決定的資料（最重要）
**Schwab & Wisse 2001「Basin of Attraction of the Simplest Walking Model」**（= Wisse 学位論文 Ch.2、
全文 PDF を pdftotext で抽出、逐語確認）が、勾配 γ に対し max|λ| と basin 面積を一枚に併記（Fig 2.8）し:

> "Clearly, there is no direct relation between the stability of the cyclic motion and its basin of attraction."
> "the most robust design would probably not be the one with the best linearized stability, but the one with the largest basin of attraction."

→ **Claim C の概念（安定性≠basin）は 2001 年に明言済み**。我々が simplest basin で既に引用している同じ論文。
我々の貢献は「点足×勾配」ではなく「**円弧足×足半径**」での定量的具体例という位置づけに留まる。

### 決定的全文テストの結果（2026-06-16 実施、両本文とも paywall 403 で未取得 → 検索要約/abstract ベース）
- **GSN 2007**（Hobbelen & Wisse, IEEE T-RO 23(6)）: 「外乱除去を slope・**foot radius**・hip spring の変化で検証」。
  既存指標（ZMP・**basin of attraction**・Floquet）を「不満足」とし、Floquet と実外乱除去の相関は限定的、
  GSN は 93% 相関。→ **Claim C の「局所安定性≠外乱耐性」を足半径軸で既に扱っている**（ただし指標は GSN で、
  basin 面積 vs 固有値の二最適分離という我々の framing とは厳密には異なる）。
- **MMT 2025/26**（「Falling dynamics」）: 「**無次元 foot radius と slope angle の相互作用**が周期歩容を支配」
  「転倒確率は **foot radius**・slope・hip spring で変わる」。SN fold が period-doubling を迂回、共存鞍点・crisis。
  → **Claim D の機構も足半径の関与も既にこの論文の射程**にある可能性が高い。「足半径→車輪極限の fold」
  という厳密な framing が含まれるかは全文未取得で未確認。

---

## 正直な結論

- **新規でない**: A（効率）、B（円弧足→basin）、**C の概念**（Schwab & Wisse 2001）、**D の機構**（MMT 2025/26）。
- **決定的テストで新規性がさらに減じた点**: C・D とも、先行研究（GSN 2007 / MMT 2025/26）が**足半径を制御軸**として
  既に扱っており、「足半径という新軸」という当初の擁護線は弱い。
- **なお残りうる狭い貢献（要全文確認）**:
  1. **E の3軸統合**（効率×固有値安定性×basin が3つの別 R で最適）を一枚に提示した点。Kwan & Hubbard は2軸まで。
  2. **basin の内部最大**（R≈0.6）が、一部文献の「半径↑で外乱耐性↑（単調）」傾向と対立しうる点。
     ただし MMT 2025/26 が foot radius×slope の転倒構造を扱う以上、ここも既出の可能性は残る。
- **F の生物・淘汰側**は文献に無く、我々の**解釈/仮説**。主張時は仮説と明示し、Heim & Spröwitz
  「Beyond Basins of Attraction」（basin を限定的指標と批判）と対話すべき。

## 位置づけ（率直に）

本 Phase の科学的価値は「**基礎的新発見**」ではなく「**既知の諸結果を、点質量 derive レイヤーを流用した
円弧足 compass で再現・統合し、固有値・効率・basin・fold を一貫した数値パイプラインで確認した**」点にある。
論文化を本気で検討するなら、(1) 下記 paywall 全文2本の精読で C・D が完全に覆われているか確定、
(2) E（3軸統合）と「内部最大 vs 単調」だけが残る貢献かを見極める——が必須。現時点の証拠は
firm な新規性主張を支持しない。

## 参照（provenance）

全文取得・逐語:
- Schwab & Wisse 2001 = Wisse PhD thesis Ch.2: https://filelist.tudelft.nl/me/Organisatie/Afdelingen/BmechE/DBL/Publications/doc/thesis_wisse.pdf
- Tedrake, Underactuated Robotics Ch.4（rimless wheel = McGeer の最簡モデル）: https://underactuated.mit.edu/simple_legs.html

abstract / 検索要約（全文未取得、要確認）:
- Adamczyk, Collins & Kuo 2006「The advantages of a rolling foot」J Exp Biol 209:3953: https://journals.biologists.com/jeb/article/209/20/3953/16394
- Adamczyk & Kuo 2013: https://pubmed.ncbi.nlm.nih.gov/23580717/
- Kuo 2002「Energetics of actively powered locomotion」: https://pubmed.ncbi.nlm.nih.gov/11871597/
- McGeer 1990a「Passive Dynamic Walking」IJRR 9(2): https://journals.sagepub.com/doi/10.1177/027836499000900206
- Kwan & Hubbard 2007「Optimal foot shape for a passive dynamic biped」JTB 248(2): https://pubmed.ncbi.nlm.nih.gov/17570405/
- Hobbelen & Wisse 2007「Gait Sensitivity Norm」IEEE T-RO 23(6): https://dl.acm.org/doi/10.1109/TRO.2007.904908
- 「Falling dynamics…」Mechanism and Machine Theory 2025/26: https://www.sciencedirect.com/science/article/abs/pii/S0094114X25004306
- Byl & Tedrake 2009「Metastable Walking Machines」IJRR 28(8): https://journals.sagepub.com/doi/10.1177/0278364909340446
- Heim & Spröwitz 2018「Beyond Basins of Attraction」arXiv:1806.08081: https://arxiv.org/pdf/1806.08081
- Full & Koditschek 1999「templates and anchors」: https://dx.doi.org/10.1242/jeb.202.23.3325

**未取得（paywall 403、新規性確定の決定的テスト）:**
- MMT 2025/26「Falling dynamics」全文（S0094114X25004306）
- GSN 2007 全文（本文・図）
- companion: 「Basin boundary metamorphoses…」Nonlinear Dynamics 2026
