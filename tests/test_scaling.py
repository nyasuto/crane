"""θ* ∝ γ^(1/3) スケーリング則の検証テスト（Garcia 1998 eq.6）。"""

import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams
from crane.search import find_limit_cycle


def test_stance_angle_scales_as_gamma_one_third():
    """θ* ∝ γ^(1/3) (Garcia 1998)。continuation で γ を変えながら追跡。"""
    gammas = [0.004, 0.006, 0.009, 0.012]
    thetas = []
    y = np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])
    # γ=0.009 から出発し、近い γ へ順に continuation
    for gamma in sorted(gammas, key=lambda g: abs(g - ref.GAMMA_REF)):
        fp = find_limit_cycle(SimplestParams(gamma=gamma), y)
        assert fp.converged, f"no convergence at gamma={gamma}"
        y = fp.y  # 次の γ の初期推定に使う
        thetas.append((gamma, fp.y[0]))

    thetas.sort()
    log_g = np.log([g for g, _ in thetas])
    log_t = np.log([t for _, t in thetas])
    slope = np.polyfit(log_g, log_t, 1)[0]
    assert abs(slope - ref.SCALING_EXPONENT) < 0.05
