"""θ* ∝ γ^(1/3) スケーリング則の検証テスト（Garcia 1998 eq.6）。"""

import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams
from crane.search import find_limit_cycle


def test_stance_angle_scales_as_gamma_one_third():
    """θ* ∝ γ^(1/3) (Garcia 1998)。単調 continuation チェーンで γ を変えながら追跡。

    distance-sorted 順（例: 0.006→0.012）では Newton 探索が非隣接不動点から始まり
    unstable branch（γ=0.012 で max|λ|≈4.41）に着地することがある。
    正しい continuation は seed 点 γ=0.009 から単調に外向きに歩む必要がある:
        upward chain:   0.009 → 0.012
        downward chain: 0.009 → 0.006 → 0.004
    各収束不動点が次の γ の初期推定になる。γ=0.009 は一度だけ計算して両チェーンで共用。
    """
    seed_y = np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT])

    # γ=0.009 の不動点を一度だけ計算（upward/downward 両チェーンの共通起点）
    seed_gamma = ref.GAMMA_REF
    seed_fp = find_limit_cycle(SimplestParams(gamma=seed_gamma), seed_y)
    assert seed_fp.converged, f"no convergence at gamma={seed_gamma}"
    assert seed_fp.eigenvalues is not None
    assert np.max(np.abs(seed_fp.eigenvalues)) < 1.0, (
        f"unstable branch at gamma={seed_gamma} — continuation jumped branches"
    )

    thetas: list[tuple[float, float]] = [(seed_gamma, seed_fp.y[0])]

    # upward chain: 0.009 → 0.012
    y = seed_fp.y.copy()
    for gamma in [0.012]:
        fp = find_limit_cycle(SimplestParams(gamma=gamma), y)
        assert fp.converged, f"no convergence at gamma={gamma}"
        assert fp.eigenvalues is not None
        assert np.max(np.abs(fp.eigenvalues)) < 1.0, (
            f"unstable branch at gamma={gamma} — continuation jumped branches"
        )
        y = fp.y.copy()
        thetas.append((gamma, fp.y[0]))

    # downward chain: 0.009 → 0.006 → 0.004
    y = seed_fp.y.copy()
    for gamma in [0.006, 0.004]:
        fp = find_limit_cycle(SimplestParams(gamma=gamma), y)
        assert fp.converged, f"no convergence at gamma={gamma}"
        assert fp.eigenvalues is not None
        assert np.max(np.abs(fp.eigenvalues)) < 1.0, (
            f"unstable branch at gamma={gamma} — continuation jumped branches"
        )
        y = fp.y.copy()
        thetas.append((gamma, fp.y[0]))

    thetas.sort()
    log_g = np.log([g for g, _ in thetas])
    log_t = np.log([t for _, t in thetas])
    slope = np.polyfit(log_g, log_t, 1)[0]
    assert abs(slope - ref.SCALING_EXPONENT) < 0.05
