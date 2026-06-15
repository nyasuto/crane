# src/crane/references_mcgeer_knees.py
"""McGeer (1990b) "Passive Walking with Knees" の公表数値。Phase 3.6 検証ゲートの照合対象。

対象モデル: 原典 Fig. 1 の ROUND-FOOT（円弧足）KNEED walker。各脚が thigh + shank の
2 リンクで膝を持ち、足裏は半径 R の円弧、hip に点質量 m_H。Fig. 2 が nominal cycle を
示し、その斜面・パラメータ凡例（Fig. 2 右）と固有値（Table 1, p1644）を印刷している。
本プロジェクトの rocker-foot kneed モデル（Phase 3 点足 kneed を R 拡張）の直接の出典。
R→0 で Phase 3（Hsu Chen 点足 kneed）に退化する。

【6 ページ会議録・スキャン PDF。数値は Fig. 2 凡例 / Table 1 / 図読み取りに分散】
原典は McGeer 1990a (IJR) の膝版会議録。詳細な定式化は本文に圧縮されており、
数値ゲートに使える印刷値は (a) Fig. 2 の nominal パラメータ凡例、(b) Table 1 の
S-to-S 固有値（"cycle of figure 2" の値と明記）、の 2 か所。不動点状態ベクトルと
step period は数値印刷されておらず、Fig. 2 の時系列から図読み取りする。

座標規約の照合（原典 vs 本プロジェクト）:
    原典 Fig. 1 / Fig. 2 caption: 「Angles ... are relative to the surface normal,
    which is tilted at 4.6% from the vertical」(p1641)。すなわち脚角は斜面法線基準の
    絶対角 — 本プロジェクトと同じ規約なので角度オフセットは不要（Goswami/Hsu Chen は
    鉛直基準で +γ 必要だったが、McGeer は最初から斜面法線基準。1990a と同じ）。
    前進 = 斜面下り +x。
    stride function の状態は (θ_C, Ω_C, Ω_FT, Ω_FS)（p1642 左下）:
        θ_C  = contact（接地時の inter-leg）角。step length を決める指定量。
        Ω_C  = contact 角速度（stance 脚の角速度に相当）
        Ω_FT = forward thigh 角速度（swing 大腿）
        Ω_FS = forward shank 角速度（swing 脛）
    本規約断面 y = (θ_st, θ̇_st, θ̇_sw) に対応:
        θ_st⁺ ≈ +θ_C（接地直後の stance 角）、θ̇_st⁺ ≈ Ω_C、θ̇_sw⁺ ≈ Ω_FT（=Ω_FS, locked）。
    速度はオフセット変換不要（θ̇_ours = θ̇_paper）。

パラメータ変換（原典 Fig. 2 凡例 → 本規約。すべて伸展脚長 l で無次元化, √(l/g) が時間単位）:
    原典 Fig. 1 が thigh/shank の幾何を明示:
        l_T = thigh 全長（hip→knee）, c_T = hip から thigh 質量点 m_T までの距離,
        l_S = shank 全長（knee→foot）, c_S = knee から shank 質量点 m_S までの距離,
        r_gyr = 各リンク CoM 回りの回転半径, w = fore/aft CoM offset, R = 足裏円弧半径。
    重要: McGeer の c は **近位関節から遠位向きに測った CoM 距離**（Fig. 1 で c_T は
    hip から下向き、c_S は knee から下向きの矢印）。これは本規約の b_t(hip→thigh 質量),
    b_s(knee→shank 質量) と **同一定義**。1990a (IJR) では c が「床からの高さ」だったが、
    1990b の膝モデルでは Fig. 1 のとおり関節基準なので変換は直結（高さ反転は不要）:
        L_THIGH (l_t, hip→knee)        ← l_T = 0.46
        L_SHANK (l_s, knee→foot)       ← l_S = 0.54   （l_t + l_s = 1.0 = 伸展脚長 l ✓）
        B_THIGH (b_t, hip→thigh 質量)   ← c_T = 0.20
        B_SHANK (b_s, knee→shank 質量)  ← c_S = 0.24
        M_THIGH ← thigh m  = 0.1
        M_SHANK ← shank m  = 0.062
        R_FOOT  ← R = 0.20
        gamma   ← Y（slope）= 0.0456 rad（凡例に印刷。caption の "4.6% from vertical"
                  = tan(0.0456) ≈ 0.0456 ≈ 4.6% と整合 ✓）
    hip 質量 m_H: Fig. 2 凡例は thigh/shank ブロックのみで hip 質量を列挙しないが、
        Fig. 2 caption が「68% of the overall mass concentrated at the hip」と印刷。
        脚 1 本 = 0.1 + 0.062 = 0.162、2 本 = 0.324。
        m_H /(m_H + 0.324) = 0.68 ⇒ m_H = 0.68/0.32 × 0.324 = 0.6885。
        検算: 総質量 = 0.6885 + 0.324 = 1.0125、hip 比 0.6885/1.0125 = 0.680 ✓。
        （m_H は「印刷された 68% + 印刷された脚質量」からの導出値。直接の数値印刷ではない。）
    r_gyr（thigh 0.135 / shank 0.186）と w（thigh 0, shank 0.01）も凡例に印刷。本実装の
        kneed モデルは各リンクを **点質量**で扱う（Phase 3 Hsu Chen と同じ）ため、リンク
        固有慣性 r_gyr は使わず記録のみ（RGYR_THIGH/RGYR_SHANK）。w≈0（対称）も記録のみ。
        注意: McGeer は分布質量（r_gyr≠0）なので、点質量近似での数値一致は限定的になりうる。
        本 Phase の一次ゲートは「cycle 存在 + R→0 で Phase 3 点足 kneed に退化 + 安定」とし、
        固有値・step period は参考照合（GAIT_TOLERANCE 参照）。
    単位: 原典は長さ l、時間 √(l/g)、質量を無次元化。本規約は SI を持つので
        L_LEG=l=1 [m], g=9.81 で次元回復。gait 照合は無次元量で行うため絶対スケールは任意。

数値の出所について（印刷値 vs 図読み取り vs 不在）:
    印刷値（Fig. 2 凡例 / caption, p1641）:
        thigh: m=0.1, r_gyr=0.135, l=0.46, c=0.20, w=0
        shank: m=0.062, r_gyr=0.186, l=0.54, c=0.24, w=0.01
        R=0.20, Y(slope)=0.0456 rad, ε_K（knee 整列の小角オフセット）≈0.2
        hip 質量比 68%（caption）→ m_H=0.6885（導出）
    印刷値（Table 1 "Step-to-step modes for the cycle of figure 2", p1644）:
        固有値 z: mode1 = −0.001, mode2 = 0.073, mode3,4 = mag 0.447 ∠±0.947 rad。
        4 個すべて |z|<1 ⇒ 安定（PUBLISHED_STABLE=True）。
        複素対 0.447·e^{±0.947i} = 0.447 cos0.947 ± i·0.447 sin0.947
              = 0.2611 ± 0.3628i（位相 0.947 は rad と判断: Table の他 phase 列 ±3.140≈π と
              整合し、phase を rad とみなすと mode3,4 の Δθ_FS 位相 ±3.140=±π = 反位相で自然）。
    図読み取り値（Fig. 2 時系列, p1641 — digitized）:
        half inter-leg（接地時 stance/contact 角）≈ 0.34（stance angle 曲線の t=0 値
            ≈ −0.34、step 末の +0.34 へ向かう。図読み取り ±10%）。
        step period（heel-strike 間隔）: 数値印刷なし。Fig. 2 で knee lock が t≈1.65 √(l/g)、
            stance angle が +0.34 に達する heel-strike は図の右端（t≈2.5）をやや超える
            ⇒ τ ≈ 2.7 √(l/g)（外挿込みの図読み取り ±15%）。
            検算: pendulum period 2π≈6.28 との比 τ/2π≈0.43 は McGeer の "step period ≈ 0.4×
            pendulum period" 級と同オーダー。
    不在（数値印刷なし → None or seed）:
        不動点状態ベクトル (θ_C, Ω_C, Ω_FT, Ω_FS): Table 1 は固有値とモード形 Δθ のみ印刷し、
            不動点そのものは印刷しない。よって SECTION_GUESS は seed（公表値ではない）。
"""

PROVENANCE = (
    "McGeer, T. (1990), 'Passive Walking with Knees', Proc. IEEE Int. Conf. on "
    "Robotics and Automation (ICRA), pp.1640-1645 (CH2876-1/90/0000/1640). "
    "Scanned PDF (no text layer), 6 pages. Public copy: "
    "http://ruina.tam.cornell.edu/research/history/mcgeer_1990_passive_walking_knees.pdf "
    "(retrieved 2026-06-16). Values read from page images: nominal parameters from "
    "the Fig. 2 legend/caption (p1641); S-to-S eigenvalues from Table 1 (p1644, "
    "'cycle of figure 2'); half inter-leg angle and step period digitized from Fig. 2."
)

# --- Fig. 2 round-foot kneed walker の nominal パラメータ（本規約へ変換。docstring 参照）---
# 各リンクは点質量近似（Phase 3 Hsu Chen と同じ）。長さは伸展脚長 l=1 で無次元化。
M_HIP: float = 0.6885  # m_H [kg]（caption「68% at hip」+ 脚質量から導出。直接印刷ではない）
M_THIGH: float = 0.1  # thigh m（Fig. 2 凡例に印刷）
M_SHANK: float = 0.062  # shank m（Fig. 2 凡例に印刷）
L_THIGH: float = 0.46  # l_T = hip→knee（Fig. 2 凡例に印刷）
L_SHANK: float = 0.54  # l_S = knee→foot（Fig. 2 凡例。l_T + l_S = 1.0 = 伸展脚長 ✓）
B_THIGH: float = 0.20  # b_t = c_T = hip→thigh 質量（Fig. 1 で hip 基準, Fig. 2 凡例に印刷）
B_SHANK: float = 0.24  # b_s = c_S = knee→shank 質量（Fig. 1 で knee 基準, Fig. 2 凡例に印刷）
R_FOOT: float = 0.20  # R = 足裏円弧半径 / l（Fig. 2 凡例に印刷）。R→0 で Phase 3 点足 kneed
G: float = 9.81  # 重力。原典は g で無次元化（数値印刷なし）。標準値で次元回復

# 各リンクの回転半径（分布質量パラメータ）。Fig. 2 凡例に印刷。本実装は点質量近似のため
# 記録のみ（数値ゲートには使わない。点質量 vs 分布質量の差は固有値照合の許容に織り込む）。
RGYR_THIGH: float = 0.135  # r_gyr,thigh（印刷）
RGYR_SHANK: float = 0.186  # r_gyr,shank（印刷）
W_THIGH: float = 0.0  # fore/aft CoM offset thigh（印刷, 対称）
W_SHANK: float = 0.01  # fore/aft CoM offset shank（印刷, ≈0）

# 公表 gait の斜面（Fig. 2 凡例 Y=0.0456。caption「4.6% from vertical」と整合）
GAMMA_GAIT: float = 0.0456  # [rad]（≈2.61°。原典が rad で直接印刷）

# 接地時 half inter-leg（contact）角。Fig. 2 の stance angle 曲線から digitize（±10%）
ALPHA_HALF_INTERLEG: float = 0.34  # [rad]（図読み取り）

# 断面 y = (θ_st, θ̇_st, θ̇_sw) の初期推定（衝突直後）。出所:
#   McGeer は不動点状態ベクトル (θ_C, Ω_C, Ω_FT, Ω_FS) を印刷していない（Table 1 は
#   固有値とモード形のみ）。よって以下は Newton shooting の SEED であり「公表値ではない」。
#   plan 指示に従い Phase 3 点足 kneed の不動点 (0.23859, −1.10959, −0.05715) を R 調整の
#   出発点とする。ただし McGeer Fig. 2 の half inter-leg は ≈0.34 と Phase 3 の 0.238 より
#   大きいので、θ_st seed は図読み取り 0.34 寄りに置く（R=0.2 と質量配置差で前進）。
#   θ̇ は McGeer 未印刷のため Phase 3 値を流用（angles は両者とも斜面法線基準でオフセット不要）。
SECTION_GUESS: tuple[float, float, float] = (0.34, -1.10, -0.06)

# 膝屈曲の符号: swing 中は脛が大腿より後方に折れる → −1（Phase 3 references_kneed.py と同じ。
# McGeer Fig. 2 の shank angle が thigh angle より先行して立ち上がり knee lock 時に交差する
# 挙動も同じ屈曲向きと整合）。
KNEE_FLEXION_SIGN: float = -1.0

# 安定性: Table 1 の固有値 4 個（−0.001, 0.073, 0.447∠±0.947）はすべて |z|<1。
# 原典本文も「passive walking requires ... a repetitive cycle, but also stability」の文脈で
# Fig. 2 cycle を安定例として提示。
PUBLISHED_STABLE: bool = True

# 印刷固有値（Table 1, "cycle of figure 2"）。状態 (θ_C, Ω_C, Ω_FT, Ω_FS) の S-to-S map:
#   mode1 z=−0.001, mode2 z=0.073, mode3,4 = 0.447·e^{±0.947i}=0.2611±0.3628i。
# 注意: (a) McGeer は分布質量（r_gyr≠0）、本実装は点質量近似のため数値一致は限定的、
# (b) 断面が 4D（θ_C, Ω_C, Ω_FT, Ω_FS）— 本実装は locked 相で Ω_FS=Ω_FT のため実効 3D に
#   縮約されうる（Phase 3 と同じ rigid-through-collision 採用時）。直接ゲートには使わず参考値。
POINCARE_EIGENVALUES: tuple[complex, ...] | None = (
    -0.001 + 0.0j,
    0.073 + 0.0j,
    0.447 * 0.58412 + 0.447 * 0.81167j,  # 0.447·e^{+0.947i} = 0.2611 + 0.3628i
    0.447 * 0.58412 - 0.447 * 0.81167j,  # 0.447·e^{-0.947i} = 0.2611 - 0.3628i
)

# step period（√(l/g) 単位）: 数値印刷なし。Fig. 2（knee lock t≈1.65, heel-strike が
# stance angle=+0.34 到達 ≈ t 2.5 超）から外挿読み取り（±15%）。τ/2π≈0.43 と整合。
STEP_PERIOD: float | None = 2.7  # [√(l/g)]（図読み取り＋外挿、±15%）

# 照合許容（相対）: パラメータ（l, c, m, R, γ）は Fig. 2 凡例に exact 印刷だが、
# (a) M_HIP は 68% からの導出、(b) ALPHA/STEP_PERIOD はスキャン図読み取り、
# (c) McGeer は分布質量 r_gyr≠0 で本実装の点質量近似と力学が厳密一致しない。
# 以上を見込んで ±15% とする（一次ゲートは cycle 存在 + R→0 退化 + 安定で、許容は
# 主に図読み取り量と点質量近似差に効く）。
GAIT_TOLERANCE: float = 0.15
