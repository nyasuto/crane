import numpy as np

from crane.basin import CONVERGED, basin_slice
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane import references as ref


def test_small_neighborhood_all_converges_simplest():
    """局所漸近安定（max|λ|=0.589<1）なら不動点の微小近傍は全て収束する。

    basin が固有値解析と整合することの内部ゲート。
    """
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    assert fp.converged
    assert np.max(np.abs(fp.eigenvalues)) < 1.0  # 前提（Phase 1 ゲート）

    # 不動点の |y*|·0.5% 程度の微小球
    hw = 0.005 * np.abs(fp.y)
    res = basin_slice(
        make_simplest,
        SimplestParams(gamma=ref.GAMMA_REF),
        fp.y,
        axes=(0, 1),
        half_widths=(float(hw[0]), float(hw[1])),
        resolution=5,
        model_name="simplest",
        n_workers=1,
    )
    assert np.all(res.grid == CONVERGED)
