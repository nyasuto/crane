import numpy as np

from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.search import find_limit_cycle
from crane import references_mcgeer as ref


def _fp_at_R(R, guess):
    p = RockerCompassParams(
        m=ref.M_LEG,
        m_h=ref.M_HIP,
        c=ref.C_HIP_TO_COM,
        rho=ref.RHO_GYR,
        L=ref.L_LEG,
        R=R,
        gamma=ref.GAMMA_GAIT,
        g=ref.G,
    )
    return find_limit_cycle(make_rocker_compass(p), np.array(guess))


def test_R030_reproduces_phase35_cycle():
    """最強アンカー: R=0.3 で Phase 3.5 検証済みサイクルを再現。"""
    fp = _fp_at_R(0.3, ref.SECTION_GUESS)
    assert fp.converged
    np.testing.assert_allclose(fp.y, [0.30844, -1.26256, -0.87914], atol=1e-3)
    assert np.max(np.abs(fp.eigenvalues)) < 1.0


def test_continuation_connects_neighboring_R():
    """連続性: R=0.3 の解を seed に R=0.25 が近傍へ連続接続する。"""
    fp030 = _fp_at_R(0.3, ref.SECTION_GUESS)
    assert fp030.converged
    fp025 = _fp_at_R(0.25, fp030.y)
    assert fp025.converged
    assert np.linalg.norm(fp025.y - fp030.y) < 0.3
