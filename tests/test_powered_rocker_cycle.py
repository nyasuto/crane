import numpy as np

from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.search import find_limit_cycle

NOMINAL = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)
GUESS = np.array([0.30844, -1.26256, -0.87914])


def test_powered_cycles_exist_and_stable_across_pushoff():
    """γ=0.030 で push_off∈{0,0.04,0.08} の安定サイクルが受動 seed から直接収束。"""
    guess = GUESS
    results = {}
    for po in [0.0, 0.04, 0.08]:
        model = make_powered_rocker_compass(PoweredRockerCompassParams(**NOMINAL, push_off=po))
        fp = find_limit_cycle(model, guess)
        assert fp.converged, f"push_off={po} not converged"
        assert np.max(np.abs(fp.eigenvalues)) < 1.0, f"push_off={po} unstable"
        results[po] = fp.y
        guess = fp.y
    # push-off で不動点が動く（受動と異なる動力サイクル）
    assert not np.allclose(results[0.0], results[0.08], atol=1e-3)
