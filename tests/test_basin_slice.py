import numpy as np

from crane.basin import CONVERGED, BasinResult, basin_slice
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane import references as ref


def _simplest_fp():
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    assert fp.converged
    return fp.y


def test_basin_slice_shape_and_center():
    y_star = _simplest_fp()
    res = basin_slice(
        make_simplest,
        SimplestParams(gamma=ref.GAMMA_REF),
        y_star,
        axes=(0, 1),
        half_widths=(0.002, 0.002),
        resolution=5,
        model_name="simplest",
        n_workers=1,
    )
    assert isinstance(res, BasinResult)
    assert res.grid.shape == (5, 5)
    # 奇数解像度なので中心セルは不動点そのもの → CONVERGED
    assert res.grid[2, 2] == CONVERGED
    assert 0.0 < res.basin_fraction <= 1.0


def test_basin_slice_serial_parallel_agree():
    y_star = _simplest_fp()
    kw = dict(
        axes=(0, 1),
        half_widths=(0.01, 0.01),
        resolution=5,
        model_name="simplest",
    )
    serial = basin_slice(
        make_simplest, SimplestParams(gamma=ref.GAMMA_REF), y_star, n_workers=1, **kw
    )
    parallel = basin_slice(
        make_simplest, SimplestParams(gamma=ref.GAMMA_REF), y_star, n_workers=2, **kw
    )
    assert np.array_equal(serial.grid, parallel.grid)
