import numpy as np

from crane import references_goswami as ref
from crane.models.compass import CompassParams, make_compass
from crane.search import find_limit_cycle


def test_stable_family_continues_below_published_gait_slope():
    """公表 slope から下向きの単調 continuation で安定 family が続く。"""
    y = np.array(ref.SECTION_GUESS)
    gammas = np.linspace(ref.GAMMA_GAIT, ref.GAMMA_GAIT * 0.6, 5)
    for gamma in gammas:
        p = CompassParams(m=ref.M_LEG, m_h=ref.M_HIP, a=ref.A, b=ref.B, gamma=float(gamma), g=ref.G)
        fp = find_limit_cycle(make_compass(p), y)
        assert fp.converged, f"no convergence at gamma={gamma}"
        assert np.max(np.abs(fp.eigenvalues)) < 1.0, f"unstable at gamma={gamma}"
        y = fp.y
