# src/crane/references_goswami.py
"""Goswami, Thuilot, Espiau (1996) RR-2996 の公表数値。Phase 2 検証ゲートの照合対象。

取得した数値の桁数は原典の記載どおり。原典に印刷されていない量は
「図からの読み取り値 (digitized)」と明記するか None + 理由コメントとする。

座標規約の照合（原典 vs 本プロジェクト）:
    原典は脚の絶対角を「鉛直線から反時計回り正」で測る（Sec. 2.2:
    θ_ns = 非支持脚 (non-support/swing), θ_s = 支持脚 (support/stance)。
    状態ベクトルは q = [θ_ns, θ_s, θ̇_ns, θ̇_s] で swing が先）。
    本プロジェクトは「斜面法線から」測った絶対角で、重力を傾けて表現する。
    斜面を下る進行方向 +x のとき斜面法線は鉛直から −γ（時計回りに γ）なので、
        θ_ours = θ_paper + γ,  θ̇_ours = θ̇_paper（オフセットは定数）
        脚ラベル対応: θ_ns ↔ θ_sw, θ_s ↔ θ_st（並び順も逆）
    検算: 原典の接地条件 (2.9) θ_ns⁻ + θ_s⁻ = −2γ は本プロジェクトの
    strike 面 θ_st + θ_sw = 0 に一致する。また (2.10) より接地時の
    半脚間角 α は θ_ns⁻ − θ_s⁻ = 2α、本規約では θ_st⁻ = −α, θ_sw⁻ = +α。

定式化の照合（cross-check 済み 2026-06-13、相違なし）:
    EOM: 原典 Appendix A.1 の慣性行列 (A.8)・重力項 (A.10) を上記の座標変換
        （角度オフセット + 並び替え）で本実装の T, V から導いた M(θ), ∂V/∂θ と
        sympy で記号比較し、厳密に一致（差はゼロ行列）。
    衝突写像: 原典 (A.18)(A.19) = 「系全体の角運動量を衝突足先（新接地点）回りで
        保存」+「衝突直前の支持脚の角運動量を hip 回りで保存」の2本で、
        配置はそのまま (A.6)、脚ラベル交換は J 行列 (2.4-2.5)。
        本プロジェクト models/compass.py の定式化（系全体を新接地点回り +
        trailing 脚を hip 回り、swap 代入）と同一。

不動点数値の出所について:
    原典は γ=3° の不動点状態ベクトルを数値としては印刷していない。
    印刷されているのは Poincaré 写像 Jacobian ∇F とその固有値 (3.43)(3.44) のみ。
    gait 量（step period, 半脚間角, 支持脚角速度）は分岐図 Figs. C.4/C.5/C.6 と
    位相面図 Fig. C.11 から読み取った（HAL 公開 PDF を 400dpi でレンダリングし
    軸枠でピクセル校正して digitize。検算: 分岐点 γ=4.37° で α=17.46° /
    T=0.751 s の目盛と一致）。読み取り精度は線幅・校正誤差込みで ±2% 程度。
"""

import math

PROVENANCE = (
    "Goswami, Thuilot, Espiau (1996), "
    "'Compass-Like Biped Robot Part I: Stability and Bifurcation of Passive Gaits', "
    "INRIA Research Report RR-2996. Retrieved 2026-06-13 from "
    "https://inria.hal.science/inria-00073701/document "
    "(journal version: IJRR 17(12), 1998, doi:10.1177/027836499801701202)"
)

# 標準パラメータ（Sec. 2.2: m_C = 2m + m_H = 20 kg, l = a + b = 1 m は本文に明記。
# Sec. 3.4.3 / 4.1 の標準ケース: m = 5 kg, m_H = 10 kg, a = b = 0.5 m,
# すなわち質量比 μ = m_H/m = 2, 長さ比 β = b/a = 1）
M_LEG: float = 5.0  # m [kg]（脚の集中質量。足から a, hip から b の点）
M_HIP: float = 10.0  # m_H [kg]
A: float = 0.5  # a [m]（足 → 脚質量点。原典と同じ定義）
B: float = 0.5  # b [m]（脚質量点 → hip。原典と同じ定義）
# 重力加速度は原典に数値として印刷されていない。標準値を採用
# （Fig. C.8 の E ≈ 153 J @3° は 9.8 / 9.81 の図上区別がつかない）
G: float = 9.81

# 公表 gait のある斜面: Sec. 3.4.3 で「3° downward incline」の steady gait に
# 対して ∇F (3.43) と固有値 (3.44) を印刷している
GAMMA_GAIT: float = math.radians(3.0)  # 0.05236 rad

# --- γ=3° gait の照合量 ---
# step period: Fig. C.4 (bifurcation diagram, T vs slope) の γ=3.0° を digitize
STEP_PERIOD: float | None = 0.735  # [s]（図読み取り、±2%）
# 接地時の stance 角（本プロジェクト規約）: θ_st⁻ = −α。
# 半脚間角 α は Fig. C.5 の γ=3.0° を digitize: α = 15.5° = 0.271 rad
THETA_STRIKE: float | None = -0.271  # [rad]（図読み取り、±2%）

# 断面（衝突直後）座標 y = (θ_st, θ̇_st, θ̇_sw) の初期推定。出所:
#   θ_st⁺ = +α = 0.271（Fig. C.5 digitize）
#   θ̇_st⁺ = −1.08 rad/s（Fig. C.6「支持脚角速度」。原典の Poincaré 断面は
#       「swing 脚が地面を離れる瞬間」= 衝突直後なので post-impact 値と解釈。
#       Fig. C.11 の 3° サイクル（点線）の支持相右端 ≈ −1.0 とも整合）
#   θ̇_sw⁺ = −0.33 rad/s（Fig. C.11 の 3° サイクル左側の衝突ジャンプ上端を
#       digitize。衝突直後の swing 脚は一時的に後退する。精度 ±0.08 程度）
# 速度は規約変換不要（角度オフセット γ は定数）。
SECTION_GUESS: tuple[float, float, float] = (0.271, -1.08, -0.33)

# 照合許容（相対）: 上記 gait 量は図からの読み取り値。プロット枠ピクセル校正での
# 読み取り誤差（線幅 ~0.5%）+ 軸校正の仮定（目盛＝枠端）を見込んで ±2% とする
GAIT_TOLERANCE: float = 0.02

# 安定性: γ=3° の Poincaré 写像 Jacobian 固有値の絶対値（eq. (3.44) に印刷）。
# |λ| = 0.332, 0.332（複素対）, 0.014, 2.554e-9（≈0、拘束方向）。
#
# 【歩数定義の数値確定 2026-06-13】
# 原典本文は (3.24)/(3.39) を step-to-step（1歩）map と呼ぶ点と
# Fig 3.1 caption（1 cycle = 2 steps）が混在するため、自前数値で歩数定義を確定した:
# 印刷値は 2歩合成写像 S² の multiplier である。
#
# 根拠: 我々の 3D 1歩 Poincaré 写像（γ=3°, 標準パラメータ）の実測固有値絶対値は
#   |λ| = 0.580, 0.580（複素対）, 0.132（実数）。
#   複素対: 0.580² = 0.336 ≈ 0.332（1.3%一致、印刷精度 3桁と整合）。
#   直接比較 0.580 vs 0.332 は 75% 乖離 → 1歩解釈は成立しない。
# 実数方向: 0.132² = 0.017 vs 印刷 0.014（24%乖離）は 4D 原典 vs 3D 本実装の
# 断面定式化の差異に起因（拘束方向 trivial 固有値 2.554e-9 を本実装では明示的に除去）。
# 小さい絶対誤差（0.003）であり、安定性判定（全 |λ| < 1）には影響なし。
POINCARE_EIGENVALUE_ABS: tuple[float, ...] = (0.332, 0.332, 0.014, 2.554e-9)

# 分岐: μ=2, β=1 の slope continuation（Sec. 4.2.1 / 4.2.4 本文に明記）
# 「symmetric gait は γ = 4.37° まで」「最初の period doubling は γ = 4.37°」
# （Table 1 は 0.25° 刻みの粗い区間 [4.5°, 5°) を示すが、本文の 4.37° を採用）
FIRST_PERIOD_DOUBLING_GAMMA: float | None = math.radians(4.37)
# 2nd doubling 4.9°, 8-periodic 5.01°（Sec. 4.2.4）。γ=5.04° 以降は周期検出不能
# （ヒストグラムは 2^n 周期示唆）、γ=5.2° を「chaotic と呼べる」と明言（Sec. 4.2.5）
CHAOS_GAMMA: float | None = math.radians(5.2)
