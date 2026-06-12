# src/crane/references.py
"""Garcia et al. 1998 の公表数値。Phase 1 検証ゲートの照合対象。

取得した数値の桁数は論文の記載どおり。論文に無い項目は None とし、
その旨をコメントに残す。

座標規約の照合（論文 vs 本プロジェクト）:
    論文の規約は θ = 斜面法線から測った stance 脚角、φ = stance 脚から測った
    swing 脚角（Fig. 1）。本プロジェクトの規約と完全に一致し、符号変換は不要。
    EOM（論文 eq.1-2, β=m/M→0 極限, l=g=1 スケーリング）:
        θ̈ = sin(θ − γ)
        φ̈ = sin(θ − γ) + θ̇² sin φ − cos(θ − γ) sin φ
    heelstrike 条件（eq.3）: φ − 2θ = 0（θ < 0 側。脚平行付近の scuffing は無視）
    heelstrike 写像（eq.4）: θ⁺ = −θ⁻, θ̇⁺ = cos(2θ⁻)·θ̇⁻,
        φ⁺ = −2θ⁻, φ̇⁺ = cos(2θ⁻)(1 − cos(2θ⁻))·θ̇⁻
    いずれも本プロジェクトの実装計画と一致（cross-check 済み 2026-06-12）。

不動点数値の出所について:
    論文は γ=0.009 の不動点を数表としては印刷していない（Fig. 2 はグラフのみ）。
    代わりに O(γ) 解析近似式の係数を6桁で印刷している（Fig. 3 キャプション
    および Appendix A.2 の表）。下記 *_THETA / *_THETA_DOT は、その印刷係数を
    γ=0.009 で評価した値（出典数値からの算術導出。論文 Fig. 2 の図示値
    θ* ≈ 0.2 と整合）。
        θ*  = Θ₀(0)·γ^(1/3) + Θ₁(0)·γ
        θ̇* = α·Θ₀(0)·γ^(1/3) + (α·Θ₁(0) + c₁)·γ
    数値積分による厳密な不動点は論文に数値として印刷されていないため、
    Phase 1 のゲートでは O(γ) 近似であることを考慮した許容誤差を設定すること。
"""

PROVENANCE = (
    "Garcia, Chatterjee, Ruina, Coleman (1998), "
    "'The Simplest Walking Model: Stability, Complexity, and Scaling', "
    "ASME J. Biomech. Eng. 120(2). Retrieved 2026-06-12 from "
    "http://ruina.tam.cornell.edu/research/topics/locomotion_and_robotics/"
    "simplest_walking/simplest_walking.pdf"
)

GAMMA_REF = 0.009

# 論文印刷の O(γ) 解析近似係数（Fig. 3 caption / Appendix A.2 table, 6桁）
# θ* ≈ Θ₀(0)·γ^(1/3) + Θ₁(0)·γ,  θ̇* = α·Θ₀(0)·γ^(1/3) + (α·Θ₁(0)+c₁)·γ
LONG_PERIOD_COEFFS = {
    "tau0": 3.812092,
    "Theta0": 0.970956,
    "Theta1": -0.270837,
    "alpha": -1.045203,
    "c1": 1.062895,
}
SHORT_PERIOD_COEFFS = {
    "tau0": 3.141592653589793,  # 論文では τ₀ = π（厳密値）
    "Theta0": 0.943976,
    "Theta1": -0.264561,
    "alpha": -1.090331,
    "c1": 0.866610,
}

# γ=0.009 の long-period gait 不動点 (heel-strike 直後の断面座標)
# 上記印刷係数を γ=0.009 で評価した O(γ) 近似値（導出値）
LONG_PERIOD_THETA: float = 0.199529  # stance angle θ* [rad]
LONG_PERIOD_THETA_DOT: float = -0.198983  # θ̇*

# γ=0.009 の short-period gait 不動点（同じく O(γ) 近似の導出値）
SHORT_PERIOD_THETA: float | None = 0.193974
SHORT_PERIOD_THETA_DOT: float | None = -0.203696

# 安定性: 論文は γ=0.009 での固有値を数値としては印刷していない
# （Fig. 4 はグラフのみ。解析近似 Jacobian は eq.27-28 に行列として印刷）。
LONG_PERIOD_EIGENVALUE_ABS: tuple[float, ...] | None = None
# Sec. 5.1: "stable period-one gait cycles for slopes of 0 < γ < 0.0151"
STABLE_GAMMA_MAX: float | None = 0.0151
# 参考: period-2 への分岐は固有値が -1 を通過して発生（Sec. 5.2）、
# γ = 0.017→0.019 で period-doubling、γ ≈ 0.019 超で歩行消失（Sec. 5.3）。

# スケーリング則 θ* ∝ γ^(1/3)（論文 eq.6、abstract でも明言）
SCALING_EXPONENT = 1.0 / 3.0
