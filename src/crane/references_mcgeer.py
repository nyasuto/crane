# src/crane/references_mcgeer.py
"""McGeer (1990a) "Passive Dynamic Walking" の公表数値。Phase 3.5 検証ゲートの照合対象。

対象モデル: 原典 §5「Steady Walking of a General 2D Biped」の STRAIGHT-LEG
ROUND-FOOT（直脚・円弧足）2D biped（Fig. 2）。これは §4 の synthetic wheel を
2 本脚 biped に展開したもので、各脚は剛体（分布質量）、足裏は半径 R の円弧、
hip に点質量。本プロジェクトの rocker-foot compass モデルの直接の出典。

【重要: §6 以降の test machine とは別モデル】
原典には 2 つのモデルがある:
  (a) §5 の直脚・円弧足 biped（解析モデル。本ファイルの対象）
  (b) §6 以降の "test machine"（膝付きの実機。5.5 m ランプ実験 + 解析）
Table 1（§8「S-to-S Modes of the Test Machine」, p71）の固有値
  Speed z = 0.70, Swing z = −0.05, Totter z = −0.83（2.5% 斜面）
および Fig. 7（"trials of the test biped", 2.5% 斜面, 5.5 m ランプ）は
**すべて test machine (b) の値**であり、直脚モデル (a) の固有値ではない。
（プラン bootstrap の「0.70/−0.65/−0.43 on 2.9% slope」は記憶由来の誤り。
原典 Table 1 を直接読むと Totter は −0.83、斜面は 2.5%。確認済み 2026-06-16。）
直脚モデル (a) の固有値は数値表として印刷されていない（複素平面プロット Fig. 9
のみ）。したがって POINCARE_EIGENVALUES = None とし、存在性 + 安定性 + R→0 で
点足 compass に戻ること、を一次ゲートとする（プラン Task 1 の指示どおり）。

座標規約の照合（原典 vs 本プロジェクト）:
    原典 Fig. 2: 脚角は「surface normal（斜面法線）から測る」(p68, eq.19 直前
    "Δθ is the rotation from the surface normal")。本プロジェクトも斜面法線基準の
    絶対角なので角度オフセットは不要（Goswami/Hsu Chen は鉛直基準だったため +γ が
    必要だったが、McGeer は最初から斜面法線基準）。前進 = 斜面下り +x。
    接地時の半脚間角は α_0（原典: "the angle subtended by the feet, at the overall
    mass center", p73 §10.3）。本規約では post-impact stance 角 θ_st⁺ ≈ +α_0、
    swing 角 θ_sw⁺ ≈ −α_0。
    速度はオフセット変換不要（θ̇_ours = θ̇_paper）。

パラメータ変換（原典 §5 + §10 の nominal 値 → 本規約。すべて脚長 l で無次元化）:
    原典 §5（p67）が直脚モデルの可変パラメータを列挙:
        1. Foot radius R
        2. Leg radius of gyration r_gyr
        3. Leg center of mass height c
        4. Fore/aft center of mass offset w
        5. Hip mass fraction m_h
    本規約 (m, m_h, c, rho, L, R, gamma, g) への対応:
        L_LEG   ← l（脚長 = hip→接地点 when upright。原典は l で無次元化。L=1 採用）
        R_FOOT  ← R（足裏円弧半径 / l）
        RHO_GYR ← r_gyr（脚 CoM 回りの回転半径 / l）
        C_HIP_TO_COM ← (l − c)。原典 c は "leg center of mass HEIGHT"（接地点=足から
            上向きの高さ, Fig. 2 で下端 R 付近から CoM までを c とラベル）。本規約の
            C_HIP_TO_COM は hip→CoM 距離なので c_ours = L − c_paper = 1 − 0.63 = 0.37。
            （検算: c が hip に近すぎると swing 脚が間に合わず cycle 消失、c を下げると
            support transfer 非効率, という §10.3 の記述は「c = 床からの高さ」と整合。
            また α_0 を "feet が overall mass center に張る角" と定義する記述も、
            mass center が脚中ほど〜やや下にある配置と整合）。
        M_LEG   ← 1.0 を採用。原典は質量を l, g とともに無次元化の基準にし
            "changing m scales the forces but doesn't change the gait"(p72) と明記。
            質量比のみが gait を決めるので脚質量を 1.0 に正規化。
        M_HIP   ← m_h（hip mass fraction）。原典 §5 の nominal 値は §10 の各
            パラメータ掃引図に明示されず（Fig. 8/9/10 の凡例は α_0, r_gyr, c, R のみ）。
            Fig. 12（§10.4 Hip Mass）が hip mass を 0 から掃引しており、§10.1〜10.3 の
            nominal は m_h ≈ 0（hip 質量なし、脚のみ）と読める。本規約では
            M_HIP = 0.0 を採用し、これは「脚 2 本のみ・hip 点質量なし」の最小構成。
            （プラン: R→0, rho→0, c=b で Goswami 点足 compass に戻る — その compass は
            m_h を持つが、McGeer §5 nominal は m_h=0。質量比は別パラメータとして扱う。）
        gamma   ← Y（slope）。R=0.3 の nominal cycle に対応する斜面を Fig. 8 から
            digitize（後述 GAMMA_GAIT）。
        w（fore/aft offset）: 本規約モデルは左右対称 leg を仮定するので w=0。原典の
            test machine は w ≈ −0.001（Fig. 7 凡例）≈ 0 で、直脚モデル nominal も
            w=0（対称）。本ファイルでは記録のみ（W_OFFSET）。
    単位: 原典は長さを l、時間を √(l/g)、質量を m で無次元化（p72）。本規約は
        SI 量を持つので L_LEG=1 [m], M_LEG=1 [kg] とし g=9.81 で次元を回復。
        gait の照合は無次元量（角度・速度比・斜面）で行うため絶対スケールは任意。

数値の出所について（印刷値 vs 図読み取りの区別）:
    印刷値（原典本文/凡例にテキストとして印字）:
        α_0 = 0.3, r_gyr = 0.32, c = 0.63（Fig. 8/9/10 凡例。p72-73）
        R nominal = 0.3（Fig. 10 凡例の固定値。p73）
        w ≈ −0.001（Fig. 7 凡例 = test machine 値, ≈0）
    図読み取り値（digitized, スキャン図のピクセル校正）:
        gamma（R=0.3 の nominal slope）: Fig. 8（SLOPE γ vs FOOT RADIUS）を
            220→6× レンダリングし、プロット枠（左端 R=0 / 右端 R=1）と y 軸目盛
            （Y=0.05 ライン = 上から 2 本目の major tick, 以下 0.01 刻み）で校正。
            R=0.3 列で SLOPE γ 曲線を読み取り γ ≈ 0.030 rad（±10%）。
            検算: 同列で STEP PERIOD 曲線 τ_0 ≈ 2.4–2.5（凡例 τ_0=2.5 ラインと整合）。
        step period τ_0 ≈ 2.5 √(l/g)（同 Fig. 8 digitize, ±10%）。
            検算: pendulum period 2π ≈ 6.28 との比 τ_0/2π ≈ 0.40 は、原典が §5/§6 で
            述べる「step period ≈ 0.36–0.39 of pendulum period」と整合。
    安定性（Fig. 9 複素平面プロット, p73）: 直脚モデルの S-to-S 固有値 3 個を
        foot radius をパラメータに複素平面へプロット。caption 明記:
        「one eigenvalue (for the 'totter mode') lies outside [unit circle for point
        feet] ... However, with other choices for c and r_gyr, walking is stable
        even on point feet.」すなわち nominal (α_0=0.3, r_gyr=0.32, c=0.63) では
        R=0（点足）で totter 固有値が単位円外 → 不安定だが、R を増やすと内側へ入り
        **円弧足 (R=0.3) では安定**。本ファイル対象は R=0.3 なので PUBLISHED_STABLE=True。
        （R→0 で点足 compass へ退化するとき、この nominal 質量配置だと totter 不安定に
        なるのは原典の重要な定性的結論。R による安定化が rocker-foot の眼目。）
"""

PROVENANCE = (
    "McGeer, T. (1990), 'Passive Dynamic Walking', "
    "International Journal of Robotics Research 9(2):62-82, "
    "doi:10.1177/027836499000900206. Scanned PDF (no text layer), 21 pages, "
    "md5 09cbf2d90d634adfdde0399094c97daf. Public copy: "
    "https://www.cs.cmu.edu/~cga/legs/McGeer1990.pdf (retrieved 2026-06-16). "
    "Values read from page images; §5 + §10 nominal parameters printed in "
    "Fig. 8/9/10 legends, gait slope/period digitized from Fig. 8."
)

# --- §5 直脚・円弧足 biped の nominal パラメータ（本規約へ変換。docstring 参照）---
# すべて脚長 l で無次元化。L_LEG=1, M_LEG=1 を基準にとる。
M_LEG: float = 1.0  # 脚質量。原典: "changing m ... doesn't change the gait"(p72) → 1.0 正規化
M_HIP: float = 0.0  # m_h（hip 点質量）。§5 nominal は hip 質量なし（Fig.12 で別途掃引）
# C_HIP_TO_COM = L − c_paper = 1 − 0.63（原典 c は床からの CoM 高さ。docstring 参照）
C_HIP_TO_COM: float = 0.37  # 無次元（hip→脚 CoM 距離 / l。L_LEG=1 なので [m] と同値）
RHO_GYR: float = 0.32  # r_gyr（脚 CoM 回りの回転半径 / l。Fig. 8/9/10 凡例に印刷）
L_LEG: float = 1.0  # l（脚長 = hip→接地点。無次元化基準）
R_FOOT: float = 0.3  # R（足裏円弧半径 / l。Fig. 10 凡例の nominal 固定値として印刷）
G: float = 9.81  # 重力。原典は g で無次元化（数値印刷なし）。標準値で次元回復

# 半脚間角（接地時、overall mass center に feet が張る角）。Fig. 8/9/10 凡例に印刷
ALPHA_HALF_INTERLEG: float = 0.3  # α_0 [rad]
# fore/aft CoM offset。対称 leg なので 0（test machine 実測 ≈ −0.001 ≈ 0, Fig. 7 凡例）
W_OFFSET: float = 0.0

# 公表 gait の斜面: R=0.3 nominal に対応する slope を Fig. 8（SLOPE γ vs R）から digitize
GAMMA_GAIT: float = 0.030  # [rad]（≈1.72°、勾配 tan≈3.0%。図読み取り ±10%。docstring 参照）

# 断面 y = (θ_st, θ̇_st, θ̇_sw) の初期推定（衝突直後）。出所:
#   McGeer は不動点状態ベクトルを印刷していない（gait は Fig. 8/9 の曲線のみ）。
#   よって以下は Newton shooting の seed であり「公表値ではない」。
#   θ_st⁺ ≈ +α_0 = 0.30（接地時 half-interleg 角。slope-normal 絶対角、前進 +x）
#   θ̇_st⁺ ≈ −1.0（Phase 2 compass 不動点 −1.092 を γ≈0.03 の浅い斜面向けに調整した seed）
#   θ̇_sw⁺ ≈ −0.40（同 compass −0.377 近傍の seed）
# 円弧足 (R=0.3) は接地点が足中心 R 上にあり点足 compass と幾何が異なるため、
# 収束後の値は seed からずれる見込み（あくまで初期推定）。
SECTION_GUESS: tuple[float, float, float] = (0.30, -1.0, -0.40)

# step period（√(l/g) 単位）: Fig. 8 の STEP PERIOD 曲線を R=0.3 で digitize
STEP_PERIOD: float | None = 2.5  # [√(l/g)]（図読み取り ±10%。τ_0/2π≈0.40 と整合）

# 安定性: Fig. 9（複素平面）caption より、nominal (α_0=0.3, r_gyr=0.32, c=0.63) は
# R=0 点足では totter 固有値が単位円外で不安定だが、R を増やすと内側へ入り
# **円弧足 R=0.3 では安定**。本ファイル対象は R=0.3 なので True。
PUBLISHED_STABLE: bool = True

# 直脚モデルの固有値は数値表として印刷されていない（Fig. 9 の複素平面プロットのみ。
# 個々の z の数値ラベルなし）。Table 1 の z=0.70/−0.05/−0.83 は test machine (膝付き)
# の値であり直脚モデルではない（docstring 参照）。よって None。
# 一次ゲートは「cycle の存在 + R=0.3 で安定 + R→0,rho→0 で点足 compass に退化」。
POINCARE_EIGENVALUES: tuple[complex, ...] | None = None

# 照合許容（相対）: GAMMA_GAIT と STEP_PERIOD はスキャン図からの読み取り値
# （枠ピクセル校正 + y 軸目盛 0.01 刻み校正、線幅・校正誤差込み）。±10% とする
# （印刷パラメータ α_0/r_gyr/c/R は exact なので、許容は主に digitized 量に効く）。
GAIT_TOLERANCE: float = 0.10

# --- rimless wheel = 最簡受動歩行モデル（issue #19 / Phase 4a.1 の概念的支柱）---
# McGeer, T. (1990) "Passive Dynamic Walking", IJRR 9(2):62-82。
# 一次資料（本文）取得状況: McGeer1990 の本文は scanned PDF（テキスト層なし、
#   https://www.cs.cmu.edu/~cga/legs/McGeer1990.pdf）でこの Task では WebFetch が
#   timeout し、新規にテキスト抽出できなかった。本ファイル docstring 記載のとおり
#   過去フェーズで page images から §4 "synthetic wheel"（剛体スポーク車輪を斜面で
#   転がす）を直脚 biped の前段モデルとして読んでいるが、本 Task では一次本文を
#   再取得できていない。
# 取得日: 2026-06-16
# 確認した記述（本 Task で実際に読めた二次資料のみ）:
#   - MIT Underactuated Robotics, Ch.4 "Simple Models of Walking and Running"
#     (Russ Tedrake), https://underactuated.mit.edu/simple_legs.html, 取得 2026-06-16。
#     直接引用: "Perhaps the simplest possible model of a legged robot, introduced as
#     the conceptual framework for passive-dynamic walking by McGeer[McGeer90], is the
#     rimless wheel." — McGeer1990 を出典に rimless wheel を最簡受動歩行モデルと位置づけ。
#   - Garcia, Chatterjee, Ruina, Coleman (1998) "The Simplest Walking Model", J.
#     Biomech. Eng. 120(2):281-288（検索要約経由、本文 PDF 直接取得は timeout）:
#     McGeer の 2D rimless wheel は転がる代わりに剛体スポークで pivot/衝突し、
#     非滑り間欠接触で並進する点を歩行と共有する、と記述。
RIMLESS_WHEEL_PROVENANCE: str | None = (
    "rimless（synthetic spoked）wheel を最簡受動歩行モデルとする位置づけは "
    "MIT Underactuated Robotics Ch.4 (Tedrake, https://underactuated.mit.edu/"
    "simple_legs.html, 取得 2026-06-16) で確認: 'Perhaps the simplest possible "
    "model of a legged robot, introduced as the conceptual framework for "
    "passive-dynamic walking by McGeer[McGeer90], is the rimless wheel.' "
    "Garcia et al. 1998 'The Simplest Walking Model' も McGeer の 2D rimless wheel を "
    "前段モデルとして言及。McGeer1990a 一次本文（scanned PDF, テキスト層なし）は "
    "本 Task では WebFetch timeout で再取得できず、§4 synthetic wheel の記述は本ファイル "
    "docstring 記載の過去フェーズ読解に依拠する（一次本文の本 Task 内引用は控える）。"
)
