# src/crane/references_kneed.py
"""Hsu Chen (2007) MIT 修士論文の公表数値。Phase 3 検証ゲートの照合対象。

点足 kneed walker（点質量 5 個: hip + 両脚の thigh/shank）。原典 Ch.3 が
本プロジェクトと同じ 4 相 hybrid（unlocked 3-link → knee-strike →
locked 2-link → heel-strike）を定義し、パラメータ（Table 3.1）・斜面角・
不動点（Sec. 3.3）・固有値（Sec. 3.3.1）を数値として印刷している。

座標規約の照合（原典 vs 本プロジェクト）:
    原典 Fig. 3-1: 「all angles are defined globally from the vertical axis」。
    q1 = stance 脚, q2 = swing 大腿, q3 = swing 脛（並び順は本プロジェクトと同じ）。
    本プロジェクトは斜面法線基準の絶対角なので、Goswami (Phase 2) と同じく
        θ_ours = θ_paper + γ,  θ̇_ours = θ̇_paper（オフセットは定数）
    検算: heel-strike 面は本規約で θ_st + θ_th = 0 ⇔ 原典で q1 + q2 = −2γ。
    印刷不動点 0.1877 + (−0.2884) = −0.1007 ≈ −2(0.0504) = −0.1008
    （4 桁丸め誤差内で一致）。さらに Sec. 3.3.2 印刷の interleg angle
    0.4761 = q1 − q2 = 0.1877 + 0.2884 とも一致。符号の向き（前進 = 斜面下り
    +x、遠位点が +x に出る向きが正）も同一であることをこの 2 検算が固定する。

パラメータ変換（原典 Table 3.1 → 本規約）:
    原典は各リンク長を質量点で分割して与える: lt = a2 + b2, ls = a1 + b1,
    L = lt + ls。「below/above point mass」の向きは重力ベクトル G (3.2c) で確定:
        V/g = m_H L c1 + m_t(ls+a2) c1 + m_s a1 c1
              + m_t(L c1 − b2 c2) + m_s(L c1 − lt c2 − b1 c3)   (c_i = cos q_i)
        ∂V/∂q1 = −(m_s a1 + m_t(ls+a2) + (m_H+m_t+m_s)L) g sin q1 = 印刷 G1 ✓
        ∂V/∂q2 = (m_t b2 + m_s lt) g sin q2 = 印刷 G2 ✓
        ∂V/∂q3 = m_s b1 g sin q3 = 印刷 G3 ✓
    すなわち swing 大腿質量は hip から b2 下、swing 脛質量は knee から b1 下:
        b_t(ours: hip→大腿質量点) = b2 = 0.325
        b_s(ours: knee→脛質量点) = b1 = 0.125
        l_t(hip→knee) = a2 + b2 = 0.5,  l_s(knee→foot) = a1 + b1 = 0.5
    stance 脚は foot から測って m_s が a1 = ls − b1, m_t が ls + a2 = L − b2 に
    あり、同一脚の同じ質量配置（本実装の stance 直線脚と同じ）。
    慣性行列でも照合: H22 = m_t b2² + m_s lt², H33 = m_s b1²,
    H13 = −m_s b1 L cos(q3−q1) — いずれも上記配置と整合 ✓。
    単位: 時間は秒（Sec. 3.3「time step of 1e-3 s」）、無次元化なし。変換不要。

定式化の照合（cross-check 済み 2026-06-13）:
    knee-strike（Sec. 3.2.1 / Appendix A）: 保存則 2 本 —
        (i) 系全体の角運動量を stance 足回り、(ii) swing 脚（大腿+脛）を hip
        回り。3 速度 → 2 速度、位置不変。本プロジェクト早見表と同一 ✓。
    heel-strike（Sec. 3.2.2, eq. 3.5）: 保存則 2 本 —
        (i) 系全体を衝突足（新接地点）回り、(ii) 衝突後 swing 脚を hip 回り。
        2×2 行列で q̇1⁺, q̇2⁺ を解き、新 swing 脛は q̇3⁺ = q̇2⁺ と
        **運動学的に仮定**（衝突を通して新 swing 脚を剛体扱いし、直後に
        大腿速度のままアンロック）。
        【相違】本プロジェクトは 3 本目に「新 swing 脛の角運動量を knee 回りで
        保存」を立て、θ̇_sh⁺ を独立に解く。保存則 (i)(ii) は一致、3 本目のみ
        モデリング選択が異なる。m_s/m_t = 0.1 なので影響は小さいはずだが
        系統差であり、GAIT_TOLERANCE の根拠に算入（勝手にどちらかへ
        合わせる変更はしない）。
    イベント条件: knee-strike =「脛が大腿と整列した瞬間」、heel-strike =
        locked 相での swing 足接地、膝が曲がったままの接地は転倒扱い
        （原典 Sec. 3.1 / Fig. 3-2）— 本実装の相機械と同じ順序・条件 ✓。

数値の出所について:
    Table 3.1（パラメータ）、γ = 0.0504 rad、heel-strike 直後の不動点
    （q1, q2=q3, q̇1, q̇2=q̇3）、interleg angle 0.4761、固有値 3 個は
    本文に数値として印刷されている。
    step period は印刷されておらず、Fig. 3-4（5 s シミュレーションの力学的
    エネルギー、横軸 time step、dt = 1e-3 s）の heel-strike 鋸歯 8 回から
    digitize した（400% レンダリング + 軸枠ピクセル校正 + PE 曲線の V 字
    ノッチ検出。strike ≈ 643, 1202, 1761, 2321, 2880, 3439, (3998), 4558 —
    7 番目は凡例と重なり検出が乱れたが両隣の間隔和が 2×559.3 で整合。
    等間隔 559.3 step → 0.559 s、読み取り ±2%）。
    g は数値として印刷されていない。Fig. 3-4 の総力学的エネルギー ≈ 12.3
    （全質量 1.6, 平均 CoM 高さ ≈ 0.76 で g·h ≈ 7.7 → g ≈ 9.8–9.81）から
    SI 重力であることは確定（g=1 なら ≈ 1.2 のはず）。9.8 と 9.81 は図上
    区別がつかないため標準値 9.81 を採用（Phase 2 references_goswami と同じ規律）。
"""

PROVENANCE = (
    "Hsu Chen, V. F. (2007), 'Passive Dynamic Walking with Knees: A Point Foot Model', "
    "M.Eng. thesis, MIT EECS (advisor: R. Tedrake). Retrieved 2026-06-13 from "
    "https://dspace.mit.edu/handle/1721.1/41635 (PDF: "
    "https://groups.csail.mit.edu/robotics-center/public_papers/Hsu07.pdf)"
)

# --- パラメータ（Table 3.1 を本規約へ変換。変換根拠は docstring 参照）---
M_HIP: float = 0.5  # m_H [kg]
M_THIGH: float = 0.5  # m_t [kg]（大腿の点質量。m_t : m_s = 10 : 1 が安定性に重要と明記）
M_SHANK: float = 0.05  # m_s [kg]
L_THIGH: float = 0.5  # l_t = a2 + b2 = 0.175 + 0.325 [m]（hip→knee）
L_SHANK: float = 0.5  # l_s = a1 + b1 = 0.375 + 0.125 [m]（knee→foot）
B_THIGH: float = 0.325  # b_t = b2 [m]（hip→大腿質量点。G2 = (m_t b2 + m_s lt) g sin q2 で確定）
B_SHANK: float = 0.125  # b_s = b1 [m]（knee→脛質量点。G3 = m_s b1 g sin q3 で確定）
# 重力: 原典に数値の印刷なし。Fig. 3-4 のエネルギー水準で SI 重力と確定し標準値を採用
G: float = 9.81

# 公表 gait のある斜面（Sec. 3.3:「a ramp with a downward angle of γ = 0.0504 rad」）
GAMMA_GAIT: float = 0.0504  # [rad]（≈ 2.89°。原典が rad で直接印刷）

# --- heel-strike 直後の不動点（Sec. 3.3 に印刷、4 桁）---
# 原典: q1 = 0.1877, q2 = q3 = −0.2884, q̇1 = −1.1014, q̇2 = q̇3 = −0.0399
# 本規約断面 y = (θ_st, θ̇_st, θ̇_sw):
#   θ_st⁺ = q1 + γ = 0.1877 + 0.0504 = 0.2381
#   （検算: θ_th⁺ = q2 + γ = −0.2380 = −θ_st⁺、断面の位置拘束 θ_th = θ_sh = −θ_st
#    が 4 桁丸め内で成立 ✓）
#   速度はオフセット変換で不変: θ̇_st⁺ = −1.1014, θ̇_sw⁺ = −0.0399
# 注意: 原典の heel-strike は q̇3⁺ = q̇2⁺ を仮定（rigid-through-collision モデル）。
# 本プロジェクトも同モデルを採用したため（plan revision 7094d2a）、θ̇_sw⁺ = θ̇_th⁺ = θ̇_sh⁺
# は構造的な等価。断面を 3D に縮約（原典の q̇2=q̇3 拘束が構造化）したためエントリは 3 個。
SECTION_GUESS: tuple[float, float, float] = (0.2381, -1.1014, -0.0399)

# 膝屈曲の符号: swing 中は θ_sh − θ_th < 0（脛が大腿より後方に折れる）→ −1
# 根拠: (1) 原典 Fig. 3-2/3-1 の屈曲膝は脛が後方へ折れ、足先が遅れて追随する
# （前進 = 遠位点が +x = 角度増加の向き、なので遅れる脛は θ_sh < θ_th）。
# (2) knee-strike は「脛が前方へ振れて大腿と整列・伸展した瞬間」(Sec. 3.1.1)、
# すなわち θ_sh − θ_th が負側から 0 へ到達する。(3) 膝が曲がったままの接地
# （θ_sh < θ_th のまま）を転倒扱いとする記述 (Sec. 3.1) とも整合。
KNEE_FLEXION_SIGN: float = -1.0

# 安定性: Sec. 3.3.1 が 3 次元縮約断面（α, q̇1, q̇2=q̇3）の first return map を
# 摂動 25 run の最小二乗で推定し「all the magnitudes are within the unit circle,
# the system is locally stable」と明言
PUBLISHED_STABLE: bool = True
# 印刷固有値（同節）: 0.4053 と複素対 −0.2129 ± 0.3454i（|λ| ≈ 0.405, 0.406, 0.406）。
# 注意: (a) Jacobian の有限差分ではなく摂動軌道の最小二乗推定（精度は粗め）、
# (b) 断面は 3D（θ_st, θ̇_st, θ̇_sw）— rigid-through-collision 採用後は本実装と同次元。
# 直接の数値ゲートには使わず参考値（比較対象の等価 compass は 0.15, 0.13±0.57i と印刷）
POINCARE_EIGENVALUES: tuple[complex, ...] = (
    0.4053 + 0.0j,
    -0.2129 + 0.3454j,
    -0.2129 - 0.3454j,
)

# --- gait 照合量 ---
# step period: 数値印刷なし。Fig. 3-4（dt=1e-3 s）の heel-strike 鋸歯間隔
# 559.3 time step から digitize（ピクセル校正、読み取り ±2%。docstring 参照）
STEP_PERIOD: float | None = 0.559  # [s]（図読み取り、±2%）
# heel-strike 直前の stance 角（本規約）: 脚交換の逆算で θ_st⁻ = θ_sw⁺ = q2⁺ + γ
# = −0.2884 + 0.0504 = −0.2380（印刷不動点からの導出値。SECTION_GUESS と独立でない
# ことに注意。interleg angle 0.4761 の半分 0.23805 とも一致）
THETA_STRIKE: float | None = -0.2380  # [rad]

# 照合許容（相対）: 不動点は 4 桁印刷だが (a) dt=1e-3 の固定刻み前進積分で
# 得た値（イベント検出精度は高いが積分器バイアス残り）、(b) heel-strike の
# 3 本目の扱いが原典（q̇3⁺=q̇2⁺ 仮定）と本実装（脛保存則）で異なる系統差、
# (c) STEP_PERIOD は図読み取り ±2%。以上を見込んで ±5% とする
GAIT_TOLERANCE: float = 0.05
